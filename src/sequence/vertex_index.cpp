//(c) 2019 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <queue>

#include "vertex_index.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "../common/config.h"


void VertexIndex::countKmers(size_t hardThreshold, int genomeSize)
{
	auto start = std::chrono::high_resolution_clock::now();
	bool option = (Parameters::get().kmerSize > 15);
	if (Parameters::get().kmerSize > 31)
	{
		throw std::runtime_error("Maximum k-mer size is 31");
	}

	Logger::get().debug() << "Hard threshold set to " << hardThreshold;
	if (hardThreshold == 0)
	{
		throw std::runtime_error("Wrong hard threshold value: " + 
								 std::to_string(hardThreshold));
	}
	Logger::get().debug() << "Started k-mer counting";
	genomeSize = (genomeSize > (int)Config::get("big_genome_threshold"))?16:1;
	Logger::get().info() << "GenomeSize " << genomeSize<<"\n";
	//BitMapsHash
	if (option)
		_hbm = new BitMapBig(hardThreshold, genomeSize);
	else
		_hbm = new BitMap(std::min(2,(int)hardThreshold/2), genomeSize);
	Logger::get().info() << "Option " << option << "\n";
	size_t cont = 0;
	/*size_t preCountSize = 1024 * 1024 * 1024;	//1G by default
	if (genomeSize > (int)Config::get("big_genome_threshold"))
	{
		preCountSize *= 4 * 4;					//16G in case of larger genomes
	}
	auto preCounters = new std::atomic<unsigned char>[preCountSize];
	for (size_t i = 0; i < preCountSize; ++i) preCounters[i] = 0;*/

	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		/*if (allReads.size() > 100000)
			break;*/
		allReads.push_back(seq.id);
		cont++;
	}

	//first pass: filling up naive hash counting filter
	if (_outputProgress) Logger::get().info() << "Counting k-mers (1/2):";
	std::function<void(const FastaRecord::Id&)> preCountUpdate =
			[ hardThreshold, this, option](const FastaRecord::Id& readId){
		if (!readId.strand()) return;
		int32_t nextKmerPos = _sampleRate;
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			if (_sampleRate > 1) //subsampling
			{
				if (--nextKmerPos > 0) continue;
				nextKmerPos = _sampleRate + (int32_t)(kmerPos.kmer.hash() % 3) - 1;
			}
			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) -
								   kmerPos.position -
								   Parameters::get().kmerSize;
			}
			(option)?_hbm->insert(kmerPos.kmer,0):_hbm->insert_alter(kmerPos.kmer,0);
		}
	};
	auto start_cont_1 = std::chrono::high_resolution_clock::now();
	processInParallel(allReads, preCountUpdate, 
					  Parameters::get().numThreads, _outputProgress);
	auto finish_cont_1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_cont_1 = finish_cont_1 - start_cont_1;
	Logger::get().info() << "Elapsed time Cont 1/2: " << elapsed_cont_1.count() << " s\n";
	//Pasos intermedios
	auto start_cont_u = std::chrono::high_resolution_clock::now();
	_hbm->upgrade(0);
	auto finish_cont_u = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_cont_u = finish_cont_u - start_cont_u;
	Logger::get().info() << "Upgrade: " << elapsed_cont_u.count() << " s\n";

	//second pass: counting kmers that have passed the filter
	if (_outputProgress) Logger::get().info() << "Counting k-mers (2/2):";
	std::function<void(const FastaRecord::Id&)> countUpdate =
			[hardThreshold, this,option]
					(const FastaRecord::Id& readId)
			{
				if (!readId.strand()) return;
				int32_t nextKmerPos = _sampleRate;
				//IterKmers: dado una lectura permite iterar por kmers.
				for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
				{
					if (_sampleRate > 1) //subsampling
					{
						if (--nextKmerPos > 0) continue;
						nextKmerPos = _sampleRate + (int32_t)(kmerPos.kmer.hash() % 3) - 1;
					}

					bool revCmp = kmerPos.kmer.standardForm();
					/*if (revCmp)
					{
						kmerPos.position = _seqContainer.seqLen(readId) -
										   kmerPos.position -
										   Parameters::get().kmerSize;
					}*/
					(option)?_hbm->insert(kmerPos.kmer,1):_hbm->update_frec(kmerPos.kmer);
				}
			};
	auto start_cont = std::chrono::high_resolution_clock::now();
	processInParallel(allReads, countUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	auto finish_cont = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_cont = finish_cont - start_cont;
	Logger::get().info() << "Elapsed time Cont 2/2: " << elapsed_cont.count() << " s\n";
	_hbm->size();
	for (size_t i = 0; i < _hbm->size(); i++)
	{
		if (!_hbm->what())
			for (size_t j = 0; j < _hbm->size(i);j++) {
				_kmerDistribution[_hbm->get_frec_int(i, j)]++;
			}
		else {
			_kmerDistribution[_hbm->getfrecIn(i)] += 1;
		}
	}
	for (int i = 0; i < 40; i++)
		Logger::get().info() << "KmerDistribution: " << i<<" "<<_kmerDistribution[i];
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	Logger::get().info() << "Elapsed time: " << elapsed.count() << " s\n";
}

namespace
{
	struct KmerFreq
	{
		Kmer kmer;
		size_t position;
		size_t freq;
	};
}

void VertexIndex::buildIndexUnevenCoverage(int minCoverage, float selectRate,
										   int tandemFreq)
{
	//_solidMultiplier = 2;
	_solidMultiplier = 1;

	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}

	//first, count the number of k-mers that will be actually stored in the index
	_kmerIndex.reserve(_kmerCounts.size() / 10);
	if (_outputProgress) Logger::get().info() << "Filling index table (1/2)";
	std::function<void(const FastaRecord::Id&)> initializeIndex = 
	[this, minCoverage, selectRate, tandemFreq] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		thread_local std::unordered_map<Kmer, size_t> localFreq;
		localFreq.clear();
		std::vector<KmerFreq> topKmers;
		topKmers.reserve(_seqContainer.seqLen(readId));

		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			kmerPos.kmer.standardForm();
			size_t freq = 1;
			_kmerCounts.find(kmerPos.kmer, freq);

			topKmers.push_back({kmerPos.kmer, (size_t)kmerPos.position, freq});
			++localFreq[kmerPos.kmer];
		}

		if (topKmers.empty()) return;
		std::sort(topKmers.begin(), topKmers.end(),
				  [](const KmerFreq& k1, const KmerFreq& k2)
				   {return k1.freq > k2.freq;});
		const size_t maxKmers = selectRate * topKmers.size();
		const size_t minFreq = std::max((size_t)minCoverage, topKmers[maxKmers].freq);

		for (auto kmerFreq : topKmers)
		{
			if (kmerFreq.freq < minFreq) break;
			if (kmerFreq.freq > _repetitiveFrequency ||
				localFreq[kmerFreq.kmer] > (size_t)tandemFreq) continue;

			ReadVector defVec((uint32_t)1, (uint32_t)0);
			_kmerIndex.upsert(kmerFreq.kmer, 
							  [](ReadVector& rv){++rv.capacity;}, defVec);
		}
	};
	processInParallel(allReads, initializeIndex, 
					  Parameters::get().numThreads, _outputProgress);
	
	_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
	size_t chunkOffset = 0;
	//Important: since packed structures are apparently not thread-safe,
	//make sure that adjacent k-mer index arrays (that are accessed in parallel)
	//do not overlap within 8-byte window
	const size_t PADDING = 1;
	for (auto& kmer : _kmerIndex.lock_table())
	{
		if (MEM_CHUNK < kmer.second.capacity + PADDING) 
		{
			throw std::runtime_error("k-mer is too frequent");
		}
		if (MEM_CHUNK - chunkOffset < kmer.second.capacity + PADDING)
		{
			_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
			chunkOffset = 0;
		}
		kmer.second.data = _memoryChunks.back() + chunkOffset;
		chunkOffset += kmer.second.capacity + PADDING;
	}

	if (_outputProgress) Logger::get().info() << "Filling index table (2/2)";
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, minCoverage, selectRate, tandemFreq] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		thread_local std::unordered_map<Kmer, size_t> localFreq;
		localFreq.clear();
		std::vector<KmerFreq> topKmers;
		topKmers.reserve(_seqContainer.seqLen(readId));

		for (const auto& kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			auto stdKmer = kmerPos.kmer;
			stdKmer.standardForm();
			size_t freq = 1;
			_kmerCounts.find(stdKmer, freq);

			++localFreq[stdKmer];
			topKmers.push_back({kmerPos.kmer, (size_t)kmerPos.position, freq});
		}

		if (topKmers.empty()) return;
		std::sort(topKmers.begin(), topKmers.end(),
				  [](const KmerFreq& k1, const KmerFreq& k2)
				   {return k1.freq > k2.freq;});
		const size_t maxKmers = selectRate * topKmers.size();
		const size_t minFreq = std::max((size_t)minCoverage, topKmers[maxKmers].freq);

		for (auto kmerFreq : topKmers)
		{
			if (kmerFreq.freq < minFreq) break;

			KmerPosition kmerPos(kmerFreq.kmer, kmerFreq.position);
			FastaRecord::Id targetRead = readId;
			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
				targetRead = targetRead.rc();
			}

			if (kmerFreq.freq > _repetitiveFrequency ||
				localFreq[kmerPos.kmer] > (size_t)tandemFreq) continue;

			//in case kmer not in index yet, creates a new vector
			//with a single element in it
			_kmerIndex.update_fn(kmerPos.kmer, 
				[targetRead, &kmerPos, this](ReadVector& rv)
				{
					if (rv.size == rv.capacity) 
					{
						Logger::get().warning() << "Index size mismatch " << rv.capacity;
						return;
					}
					size_t globPos = _seqContainer
							.globalPosition(targetRead, kmerPos.position);
					rv.data[rv.size].set(globPos);
					++rv.size;
				});
		}
	};
	processInParallel(allReads, indexUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	_kmerCounts.clear();
	_kmerCounts.reserve(0);

	Logger::get().debug() << "Sorting k-mer index";
	for (const auto& kmerVec : _kmerIndex.lock_table())
	{
		std::sort(kmerVec.second.data, kmerVec.second.data + kmerVec.second.size,
				  [](const IndexChunk& p1, const IndexChunk& p2)
				  	{return p1.get() < p2.get();});
	}
	
	size_t totalEntries = 0;
	for (const auto& kmerRec : _kmerIndex.lock_table())
	{
		totalEntries += kmerRec.second.size;
	}
	Logger::get().debug() << "Selected k-mers: " << _kmerIndex.size();
	Logger::get().debug() << "Index size: " << totalEntries;
}

namespace
{
	template <class T>
	size_t getFreq(T& histIter)
		{return histIter->first;};

	template <class T>
	size_t getCount(T& histIter)
		{return histIter->second;};

}

void VertexIndex::setRepeatCutoff(int minCoverage)
{
	size_t totalKmers = 0;
	size_t uniqueKmers = 0;
	for (auto mapPair = this->getKmerHist().begin();
		 mapPair != this->getKmerHist().end(); ++mapPair)
	{
		if (minCoverage <= (int)getFreq(mapPair))
		{
			totalKmers += getCount(mapPair) * getFreq(mapPair);
			uniqueKmers += getCount(mapPair);
		}
	}
	float meanFrequency = (float)totalKmers / (uniqueKmers + 1);
	_repetitiveFrequency = (float)Config::get("repeat_kmer_rate") * meanFrequency;
	
	size_t repetitiveKmers = 0;
	for (auto mapPair = this->getKmerHist().rbegin();
		 mapPair != this->getKmerHist().rend(); ++mapPair)
	{
		if (getFreq(mapPair) > _repetitiveFrequency)
		{
			repetitiveKmers += getCount(mapPair);
		}
	}
	float filteredRate = (float)repetitiveKmers / uniqueKmers;
	Logger::get().debug() << "Repetitive k-mer frequency: " 
						  << _repetitiveFrequency;
	Logger::get().debug() << "Filtered " << repetitiveKmers 
						  << " repetitive k-mers (" <<
						  filteredRate << ")";
}

void VertexIndex::buildIndex(int minCoverage)
{
	auto start = std::chrono::high_resolution_clock::now();
	if (_outputProgress) Logger::get().info() << "Filling index table";
	_solidMultiplier = 1;
	
	//"Replacing" k-mer couns with k-mer index. We need multiple passes
	//to avoid peaks in memory usage during the hash table extensions +
	//prevent memory fragmentation
	size_t kmerEntries = 0,kmerEntries2 = 0;
	size_t solidKmers = 0, solidKmers2 = 0;
	_hbm->update(true);
	Logger::get().info() << "MinCoverage: "<<minCoverage<< " Repetitive freq: "<<_repetitiveFrequency
				<<" Size: "<<(_hbm->size())<<"\n";

	size_t chunkOffset = 0, wasted = 0, chunk = 0;
	if (_hbm->what())
		_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
    const size_t PADDING = 1;
	for (size_t i = 0; i < _hbm->size(); i++) {
		for (size_t j = 0; j < _hbm->size(i); j++) {
			size_t frec = (_hbm->what()) ? _hbm->getfrecIn(i) : _hbm->get_frec_int(i, j);
			if ((size_t) minCoverage <= frec &&
				frec <= _repetitiveFrequency) {
				Kmer kmer((_hbm->what())?_hbm->return_kmer(i+1):
						  _hbm->return_kmer(i, j));
				kmerEntries2 += frec;
				++solidKmers2;++solidKmers;
				if (_hbm->what())
				{
					if (MEM_CHUNK < frec + PADDING) {
						throw std::runtime_error("k-mer is too frequent");
					}
					if (MEM_CHUNK - chunkOffset < frec + PADDING) {
						wasted += MEM_CHUNK - chunkOffset;
						chunk++;
						_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
						chunkOffset = 0;
					}
				}
				_hbm->mod_frec_pos(kmer, nullptr);
                chunkOffset += frec + PADDING;
			}
			if (frec > _repetitiveFrequency)
			{
				Kmer kmer((_hbm->what())?_hbm->return_kmer(i+1):
						  _hbm->return_kmer(i, j));
				_repetitiveKmers.insert(kmer, true);
			}
		}
	}

	Logger::get().info() << "Sampling rate: " << _sampleRate;
	Logger::get().info() << "Solid k-mers: " << solidKmers;
	Logger::get().debug() << "K-mer index size: " << kmerEntries;
	Logger::get().debug() << "Mean k-mer frequency: " 
		<< (float)kmerEntries / solidKmers;
	//Segunda actualización
	_hbm->update(false);
	Logger::get().info() << "Size: "<<_hbm->size()<<"\n";
	//Segunda pasada sobre kmers para llenar los huecos
	chunkOffset = 0;solidKmers = 0;
	size_t ind = 0;
	for (size_t i = 0; i < _hbm->size(); i++) {
		for (size_t j=0; j < _hbm->size(i);j++) {
			size_t frec = (_hbm->what()) ? _hbm->getfrecIn(i) : _hbm->get_frec_int(i, j);
			if ((size_t) minCoverage <= frec &&
				frec <= (size_t)_repetitiveFrequency) {
				Kmer kmer((_hbm->what())?_hbm->return_kmer(i+1):
						  _hbm->return_kmer(i, j));
                if (MEM_CHUNK - chunkOffset < frec + PADDING) {
					ind++;
					chunkOffset = 0;
				}
				_hbm->mod_frec_pos(kmer,(_hbm->what())?_memoryChunks[ind]+chunkOffset:new IndexChunk[0]);
                chunkOffset += frec + PADDING;
			}
		}
	}
	//Tercera Pasada
	if (!_hbm->what()) {
		_hbm->upgrade(true);
		ind = 0;chunkOffset = 0;
		for (size_t i = 0; i < _hbm->size(); i++) {
			for (size_t j = 0; j < _hbm->size(i); j++) {
				size_t frec = (_hbm->what()) ? _hbm->getfrecIn(i) : _hbm->get_frec_int(i, j);
				if ((size_t) minCoverage <= frec &&
					frec <= (size_t) _repetitiveFrequency) {
					Kmer kmer((_hbm->what()) ? _hbm->return_kmer(++solidKmers) :
							  _hbm->return_kmer(i, j));
					_hbm->FinalStep(kmer, frec);
				}
			}
		}
	}
	//Eliminamos sobrante
	_hbm->clear();
	/*
	 * Last
	 */
	if (!_hbm->what()) {
		_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
		ind = 0;chunkOffset = 0;
		for (size_t i = 0; i < _hbm->size_total(); i++)
		{
			size_t frec = _hbm->get_frec_solid(i);
			if (frec > 0) {
				if (MEM_CHUNK < frec + PADDING) {
					throw std::runtime_error("k-mer is too frequent");
				}
				if (MEM_CHUNK - chunkOffset < frec + PADDING) {
					wasted += MEM_CHUNK - chunkOffset;
					chunk++;
					_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
					chunkOffset = 0;
					ind++;
				}
				_hbm->is_unique(i, _memoryChunks[ind] + chunkOffset);
				chunkOffset += frec + PADDING;
			}else{
				size_t num_chunks = _hbm->size_total_chunk(i);
				for (size_t j = 0; j < num_chunks; j++)
				{
					frec = _hbm->get_frec_solid(i,j);
					if (MEM_CHUNK < frec + PADDING) {
						throw std::runtime_error("k-mer is too frequent");
					}
					if (MEM_CHUNK - chunkOffset < frec + PADDING) {
						wasted += MEM_CHUNK - chunkOffset;
						chunk++;
						_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
						chunkOffset = 0;
						ind++;
					}
					_hbm->is_unique(i,_memoryChunks[ind]+chunkOffset,j);
					chunkOffset += frec + PADDING;
				}
			}

		}
	}
	//Checking
	Logger::get().info() << "Solid kmers: " << solidKmers2;
	Logger::get().info() << "Kmer index size: " << kmerEntries2;
	Logger::get().info() << "Num chunks: "<<chunk<<" "<<wasted<<"\n";

	uint32_t conter = 0;
	std::function<void(const FastaRecord::Id&)> indexUpdate =
			[this, minCoverage, &conter] (const FastaRecord::Id& readId)
			{
				//std::cout << "ReadId: "<<readId.getId()<<std::endl;
				if (!readId.strand()) return;
				int32_t nextKmerPos = _sampleRate;
				//kmerPos: kmerIterator: seq_read,position y kmer (atributos)
				for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
				{
					if (_sampleRate > 1) //subsampling
					{
						if (--nextKmerPos > 0) continue;
						nextKmerPos = _sampleRate + (int32_t)(kmerPos.kmer.hash() % 3) - 1;
					}

					FastaRecord::Id targetRead = readId;
					bool revCmp = kmerPos.kmer.standardForm();

					int pos = _hbm->check(kmerPos.kmer, true);
					if (pos != -1) {
						conter++;
						size_t offset = kmerPos.position, rc_offset;
						if (revCmp)
						{
							kmerPos.position = _seqContainer.seqLen(readId) -
											   kmerPos.position -
											   Parameters::get().kmerSize;
							targetRead = targetRead.rc();
							rc_offset = kmerPos.position;
						}else
							rc_offset = _seqContainer.seqLen(readId)
										- kmerPos.position
										- Parameters::get().kmerSize;
						_read_positions_vector[readId.getId()][offset] = 1;
						if (offset != 0 || revCmp) {
                            /*if (rc_offset < 0 || rc_offset > _read_positions_vector[readId.getId()+1].size()) {
                                std::cout << rc_offset << " " << _read_positions_vector[readId.getId() + 1].size()
                                          << std::endl;
                                exit(1);
                            }*/
							_num_kmers[readId.getId()+1]++;
							_read_positions_vector[readId.getId()+1][rc_offset] = 1;
						}
						_num_kmers[readId.getId()]++;
						size_t globPos = _seqContainer.globalPosition(targetRead, kmerPos.position);
						_hbm->update_read_pos(kmerPos.kmer,targetRead,globPos,pos);
					}
				};
			};
	size_t cont = 0;
	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		while (cont < seq.id.getId()) {
			update_pos_bit(0);
			cont++;
		}
		update_pos_bit(_seqContainer.seqLen(seq.id));
		cont++;
		/*if (allReads.size() > 100000)
			break;*/
		allReads.push_back(seq.id);
	}
	auto start_cont = std::chrono::high_resolution_clock::now();
	processInParallel(allReads, indexUpdate, 
					  Parameters::get().numThreads, _outputProgress);
	auto finish_cont = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_cont = finish_cont - start_cont;
	Logger::get().info() << "Elapsed time Cont Indexing: " << elapsed_cont.count() << " s\n";

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	Logger::get().info() << "Elapsed time: " << elapsed.count() << " s\n";

	Logger::get().debug() << "Sorting k-mer index";
	_hbm->sortIndexChunks();
}


void VertexIndex::clear()
{
	for (auto& chunk : _memoryChunks) delete[] chunk;
	_memoryChunks.clear();

	_kmerIndex.clear();
	_kmerIndex.reserve(0);

	_kmerCounts.clear();
	_kmerCounts.reserve(0);
}

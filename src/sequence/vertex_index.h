//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>
#include <iostream>
#include <cstring>

#include <cuckoohash_map.hh>

#include "kmer.h"
#include "sequence_container.h"
#include "../common/config.h"
#include "../common/logger.h"
#include "Optionk30.h"

class VertexIndex
{
public:
	~VertexIndex()
	{
		this->clear();
	}
	VertexIndex(const SequenceContainer& seqContainer, int sampleRate):
		_seqContainer(seqContainer), _outputProgress(false), 
		_sampleRate(sampleRate), _repetitiveFrequency(0),
		_solidMultiplier(1)
		//_flankRepeatSize(flankRepeatSize)
	{}

	VertexIndex(const VertexIndex&) = delete;
	void operator=(const VertexIndex&) = delete;

private:

	//static const size_t MAX_INDEX = 1ULL << (sizeof(IndexChunk) * 8);

	struct ReadPosition
	{
		ReadPosition(FastaRecord::Id readId = FastaRecord::ID_NONE, 
					 int32_t position = 0):
			readId(readId), position(position) {}
		FastaRecord::Id readId;
		int32_t position;
	};

	struct ReadVector
	{
		ReadVector(uint32_t capacity = 0, uint32_t size = 0, IndexChunk * data = nullptr):
			capacity(capacity), size(size), data(data) {}
		uint32_t capacity;
		uint32_t size;
		IndexChunk* data;
	};

public:
	typedef std::map<size_t, size_t> KmerDistribution;

	class KmerPosIterator
	{
	public:
		KmerPosIterator(ReadVector rv, size_t index, bool revComp, 
						const SequenceContainer& seqContainer, bool iterhelper = false):
			rv(rv), index(index), revComp(revComp), 
			seqContainer(seqContainer), kmerSize(Parameters::get().kmerSize), iterhelper(iterhelper)
		{}

		bool operator==(const KmerPosIterator& other) const
		{
			return index == other.index && rv.data == other.rv.data;
		}
		bool operator!=(const KmerPosIterator& other) const
		{
			return !(*this == other);
		}

		//__attribute__((always_inline))
		ReadPosition operator*() const
		{
			size_t globPos = rv.data[index].get();
			FastaRecord::Id seqId;
			int32_t position;
			int32_t seqLen;
			seqContainer.seqPosition(globPos, seqId, position, seqLen);
			if (!revComp)
			{
				return ReadPosition(seqId, position);
			}
			else
			{
				return ReadPosition(seqId.rc(), seqLen - position - kmerSize);
			}
		}

		KmerPosIterator& operator++()
		{
			++index;
			return *this;
		}
	
	private:
		ReadVector rv;
		size_t index;
		bool   revComp,iterhelper;
		const  SequenceContainer& seqContainer;
		size_t kmerSize;
	};

	class IterHelper
	{
	public:
		IterHelper(ReadVector rv, bool revComp, 
				   const SequenceContainer& seqContainer): 
			rv(rv), revComp(revComp), seqContainer(seqContainer)
			{}

		KmerPosIterator begin()
		{
			return KmerPosIterator(rv, 0, revComp, seqContainer, true);
		}

		KmerPosIterator end()
		{
			return KmerPosIterator(rv, rv.size, revComp, seqContainer, true);
		}

	private:
		ReadVector rv;
		bool revComp;
		const SequenceContainer& seqContainer;
	};

	void countKmers(size_t hardThreshold, int genomeSize);
	void setRepeatCutoff(int minCoverage);
	void buildIndex(int minCoverage);
	void buildIndexUnevenCoverage(int minCoverage, float selectRate, 
								  int tandemFreq);
	std::pair<size_t,select_support_mcl<1>> get_select_structure(FastaRecord::Id readId) const
	{
		return std::pair<size_t, select_support_mcl<1>>(_num_kmers[readId.getId()]
				,select_support_mcl<1>(&(_read_positions_vector[readId.getId()])));
	}
	void clear();

	size_t check(Kmer kmer) const {
		return _hbm->check(kmer,true);
	}

	Read_Pos * get_pos(Kmer kmer,size_t pos) {
		return _hbm->get_pos(kmer, pos);
	}

	IterHelper iterKmerPos(Kmer kmer, size_t pos) const
	{
		bool revComp = kmer.standardForm();
		std::pair<size_t,IndexChunk *> pair = _hbm->get_kmer_reads(kmer,pos);
		ReadVector rv(pair.first,pair.first,pair.second);
		return IterHelper(rv, revComp,
						  _seqContainer);
	}

	//__attribute__((always_inline))
	/*bool isSolid(Kmer kmer) const
	{
		kmer.standardForm();
		return _kmerIndex.contains(kmer);
	}*/

	bool isRepetitive(Kmer kmer) const
	{
		kmer.standardForm();
		return _repetitiveKmers.contains(kmer);
	}
	
	size_t kmerFreq(Kmer kmer) const
	{
		kmer.standardForm();
		size_t freq = _hbm->getfrec(kmer, _hbm->check(kmer, true));
		/*ReadVector rv;
		_kmerIndex.find(kmer, rv);*/
		return freq;
	}

	void outputProgress(bool set) 
	{
		_outputProgress = set;
	}

	const KmerDistribution& getKmerHist() const
	{
		return _kmerDistribution;
	}

	int getSampleRate() const {return _sampleRate * _solidMultiplier;}

	//New Version

	void update_pos_bit(size_t length){
		_read_positions_vector.push_back(bit_vector(length,0));
		_num_kmers.push_back(0);
	}

	void _construct_select(FastaRecord::Id readId){
		_read_positions_select.push_back(
				select_support_mcl<1>(&_read_positions_vector[readId.getId()]));
	}

private:
	void addFastaSequence(const FastaRecord& fastaRecord);

	const SequenceContainer& _seqContainer;
	KmerDistribution 		 _kmerDistribution;
	bool    _outputProgress;
	int32_t _sampleRate;
	size_t  _repetitiveFrequency;
	int32_t _solidMultiplier;

	const size_t MEM_CHUNK = 32 * 1024 * 1024 / sizeof(IndexChunk);
	std::vector<IndexChunk*> _memoryChunks;

	cuckoohash_map<Kmer, ReadVector> _kmerIndex;
	cuckoohash_map<Kmer, size_t> 	 _kmerCounts;
	cuckoohash_map<Kmer, char> 	 	 _repetitiveKmers;

	//Try mix
	std::vector<bit_vector> _read_positions_vector;
	std::vector<size_t> _num_kmers;
	std::vector<select_support_mcl<1>> _read_positions_select;

	//Parte de ellos
	Control * _hbm;
};

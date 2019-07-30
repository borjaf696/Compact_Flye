#include "Optionk30.h"
#include "../common/config.h"
#include "math.h"
#include "unistd.h"
#include <chrono>

BitMap::BitMap(int hardThreshold, size_t genome):frec(){
    genome_size*=genome;
    hardThreshold = (hardThreshold==0)?1:hardThreshold;
    for (int i = 0; i < hardThreshold; i++)
        b_2.push_back(bit_vector(genome_size));
}

void BitMap::clear(){
    this->frec.clear();
    this->frec.reserve(0);
    (&this->b)->~bit_vector();
    b_rank = rank_support_v<1>(&_b_final);
    (&b_rank_tmp)->~rank_support_v<1>();
}

void BitMap::mod_frec_pos(Kmer kmer_pos, IndexChunk * read_pos) {
    if (read_pos == nullptr) {
        _b_final[kmer_pos.getRepresentation()] = 1;
        this->read_pos.push_back(Read_Pos());
    }else
        this->read_pos[b_rank_tmp(kmer_pos.getRepresentation())] = Read_Pos(read_pos);
}

void BitMap::update_read_pos(Kmer kmer_pos, FastaRecord::Id
        ,size_t pos, size_t) {
    uint32_t place = b_rank(kmer_pos.getRepresentation());
    read_pos[place].data[read_pos[place].size++].set(pos);
}

void BitMap::set_to_zero(Kmer) {}

int BitMap::check(Kmer kmer, bool standar) {
    if (standar)
        kmer.standardForm();
    if (kmer.getRepresentation() >= _b_final.size())
        return -1;
    return (_b_final[kmer.getRepresentation()]==1)?1:-1;
}

void BitMap::upgrade(bool) {
    for (size_t i = 0; i < (b_2.size()-1);i++)
        (&this->b_2[i])->~bit_vector();
    b = b_2[b_2.size()-1];
    b_2.clear();
    b_rank = rank_support_v<1>(&b);
}

void BitMap::update_frec(Kmer kmer) {
    size_t kmer_repr = kmer.getRepresentation();
    if (b[kmer_repr])
        frec[b_rank(kmer_repr)]++;
}

uint32_t BitMap::getfrec(Kmer kmer, size_t){
    return read_pos[b_rank(kmer.getRepresentation())].size;
}

uint32_t BitMap::getfrecIn(size_t pos){
    return frec[pos];
}

unsigned int BitMap::get_sumTotal(){
    return _sumTotal;
}

unsigned int BitMap::insert(Kmer) {
    return 0;
}

void BitMap::insert_alter(Kmer kmer, bool type_insert) {
    _insert_alter(kmer,type_insert);
}

int BitMap::locate(Kmer) {
    return 0;
}

void BitMap::setChunckSize(){}

size_t BitMap::size(){return frec.size();}

void BitMap::update(bool choose){
    if (choose) {
        b_select = select_support_mcl<1>(&this->b);
        _b_final = bit_vector(b_select(b_rank(genome_size))+1);
    }else
        b_rank_tmp = rank_support_v<1>(&this->_b_final);
}

Kmer BitMap::return_kmer(int kmer_){
    //Falta select
    size_t kmer_int = this->b_select((size_t)kmer_);
    return Kmer(kmer_int);
}

Read_Pos * BitMap::get_pos(Kmer kmer, size_t) {
    return &(read_pos[kmer.getRepresentation()]);
}

std::pair<size_t,IndexChunk *> BitMap::get_kmer_reads(Kmer kmer, size_t){
    size_t place = b_rank(kmer.getRepresentation());
    return std::pair<size_t,IndexChunk*>(read_pos[place].size,
                                           read_pos[place].data);
}

//BitMap::BitMapBig class
BitMapBig::BitMapBig(int hardThreshold, size_t genome){
    std::cout<<"HardThreshold: "<<hardThreshold<<std::endl;
    _h_threshold = hardThreshold;
    genome_size = genome_size*genome;
    _big_genome = (genome > 1);
    bloom = int_vector<NUM_BITS>(genome_size,0);
    /*for (int i = 0; i < hardThreshold; i++)
        b_2.push_back(bit_vector(genome_size,0));*/
    b = bit_vector(genome_size,0);
}

void BitMapBig::mod_frec_pos(Kmer kmer_local, IndexChunk * rp) {
    size_t place = kmer_local.hash() % (size_t) genome_size;
    if (rp == nullptr){
        if (_b_final[place] == 0) {
            _b_final[place] = 1;
            _cont_local++;
        }
    }else{
        place = b_rank_tmp(place);
        //Nueva Etapa
        if (_uniques[place]==1){
            if (_not_unique[place] == 0) {
                _not_unique[place] = 1;
                _keys_not_unique.push_back((size_t *)
                                                   malloc(_block * sizeof(size_t)));
                _not_unique_read_pos.push_back((Read_Pos *)
                                                       malloc(_block * sizeof(Read_Pos)));
                _info_new_app.push_back({_block, 0});
            }
        }else
            _uniques[place] = 1;
    }
}

void BitMapBig::FinalStep(Kmer kmer_local, size_t frec) {
    size_t place = b_rank_tmp(kmer_local.hash() % (size_t) genome_size);
    if (_not_unique[place] == 0) {
        place -= b_rank_uniques(place);
        _uniques_keys[place] = kmer_local.getRepresentation();
        _unique_read_pos[place] = Read_Pos(frec);
    } else {
        place = b_rank_uniques(place);
        _keys_not_unique[place][_info_new_app[place].second]
                = kmer_local.getRepresentation();
        _not_unique_read_pos[place][_info_new_app[place].second++]
                = Read_Pos(frec);
        if (_info_new_app[place].first == _info_new_app[place].second) {
            _info_new_app[place].first += _block;
            _keys_not_unique[place] = (size_t *) realloc(_keys_not_unique[place],
                                                         _info_new_app[place].first * sizeof(size_t));
            _not_unique_read_pos[place] = (Read_Pos *) realloc(_not_unique_read_pos[place],
                                                               _info_new_app[place].first * sizeof(Read_Pos));
        }
    }
}

void BitMapBig::update_read_pos(Kmer kmer_local, FastaRecord::Id
        ,size_t position, size_t pos) {
    size_t place = b_rank_tmp(kmer_local.hash() % (size_t) genome_size);
    size_t pos_r = b_rank_uniques(place);
    if (_not_unique[place] == 0) {
        _unique_read_pos[place - pos_r].data[_unique_read_pos[place - pos_r].size++].set(position);
    }else {
        _not_unique_read_pos[pos_r][pos].data[_not_unique_read_pos[pos_r][pos].size++].set(position);
    }
}

void BitMapBig::erase_frec(uint32_t i,std::vector<uint32_t> vect_j){
    uint8_t cont = 1;
    for (int i = vect_j[0]; i < _info_frecs[i].second-1;i++) {
        while (i+cont == vect_j[cont])
            cont++;
        _frecs[i] = _frecs[i + cont];
    }
    _info_frecs[i].second -= vect_j.size();

}

void BitMapBig::clear(){
    sdsl::util::clear(b_rank);
    //Remove defa
    //(&b)->~bit_vector();
    size_t cont = 0;
    for (auto f: _frecs){
        if (f != NULL)
            cont++;
        free(f);
    }
    for (auto k: _keys)
        free(k);
    std::cout<<"Clean Done: "<<cont<<std::endl;
    size_t * ej = (size_t*) malloc(cont*sizeof(size_t));
    free(ej);
    sdsl::util::clear(_frec_unique);
    _frec_unique.clear();_frec_unique.reserve(0);_frec_unique.shrink_to_fit();
    _frecs.clear();_frecs.reserve(0);_frecs.shrink_to_fit();
    _info_frecs.clear();_info_frecs.reserve(0);_info_frecs.shrink_to_fit();
    _keys_unique.clear();_keys_unique.reserve(0);_keys_unique.shrink_to_fit();
    _keys.clear();_keys.reserve(0);_keys.shrink_to_fit();
}

void BitMapBig::insert(Kmer kmer, bool type_insert) {
    if (type_insert)
        _insert_2(kmer);
    else {
        _insert_1_2(kmer);
    }
}

void BitMapBig::upgrade(bool type_update) {
    if (type_update) {
        sdsl::util::clear(_uniques);
        b_rank_uniques = rank_support_v<1>(&_not_unique);
        size_t size = _cont_local-b_rank_uniques(_cont_local);
        _uniques_keys = std::vector<size_t>(size,0);
        _unique_read_pos = std::vector<Read_Pos>(size);
    } else {
        /*for (size_t i = 0; i < (b_2.size());i++)
            (&this->b_2[i])->~bit_vector();
        b_2.clear();
        b_2.reserve(0);
        b_2.shrink_to_fit();*/
        for (size_t i = 0; i < bloom.size(); ++i)
        {
            if (bloom[i] >= _h_threshold) {
                b[i] = 1;
                _sumsup++;_sum_unique++;
            }
        }
        sdsl::util::clear(bloom);
        b_rank = rank_support_v<1>(&b);
        size_t max = b_rank(genome_size);
        std::cout << "Num holes: "<<max<<" "<<_sumsup<<std::endl;
        _uniques_frec = bit_vector(max,0);
        _non_uniques_frec = bit_vector(max,0);
        _frecs = std::vector<uint16_t *>(max,NULL);
        _frec_unique = std::vector<uint16_t>(max);
        _keys = std::vector<size_t*>(max,NULL);
        _keys_unique = std::vector<size_t>(max);
        _info_frecs = std::vector<std::pair<uint8_t,uint8_t>>(max,{0,0});
    }
}

void BitMapBig::update(Kmer){}

unsigned int BitMapBig::get_sumTotal(){
    return _sumsup;
}

size_t BitMapBig::size(){
    return _sumsup;
}

size_t BitMapBig::size(int i){
    //return (size_t)frecs[i].size();
    if (_non_uniques_frec[i] == 0)
        return 1;
    else
        return (size_t)_info_frecs[i].second;
}

//Iterators in stead
uint32_t BitMapBig::get_frec_int(int i,int j){
    if (_non_uniques_frec[i] == 0)
        return _frec_unique[i];
    return _frecs[i][j];
}

uint32_t BitMapBig::getfrec(Kmer kmer, size_t pos){
    size_t place = b_rank_tmp(kmer.hash() % (size_t) genome_size);
    size_t pos_r = b_rank_uniques(place);
    if (_not_unique[place] == 0)
        return _unique_read_pos[place - pos_r].size;
    else
        return _not_unique_read_pos[pos_r][pos].size;
}

void BitMapBig::_insert_2(Kmer kmer) {
    size_t kmer_repr = kmer.getRepresentation();
    size_t bucket = kmer.hash() % (size_t) genome_size;
    if (b[bucket] == 1) {
        size_t place = b_rank(bucket);
        if (_uniques_frec[place] == 0 && _non_uniques_frec[place] == 0){
            _keys_unique[place] = kmer_repr;
            _frec_unique[place] = 1;
            _uniques_frec[place] = 1;
        }else if (_uniques_frec[place] == 1 && kmer_repr == _keys_unique[place]) {
            if (_frec_unique[place] < (pow(2,NUM_BITS_FRECS)-1))
                _frec_unique[place]++;
        }else if (_uniques_frec[place] == 0){
            std::function<int(size_t)> find = [this, place] (size_t kmer_repr){
                for (int i = 0; i < _info_frecs[place].second; i++)
                    if (_keys[place][i] == kmer_repr) return i;
                return (int)_info_frecs[place].second;
            };
            int _pos = find((size_t)kmer_repr);
            if (_pos != _info_frecs[place].second) {
                if (_frecs[place][_pos] < (pow(2,NUM_BITS_FRECS)-1))
                    _frecs[place][_pos]++;
            }else {
                if (_info_frecs[place].first == _info_frecs[place].second){
                    _info_frecs[place].first += _block;
                    _frecs[place] = (uint16_t *)realloc(_frecs[place]
                            ,_info_frecs[place].first*sizeof(uint32_t));
                    _keys[place] = (size_t*)realloc(_keys[place]
                            ,_info_frecs[place].first*sizeof(size_t));
                    _initialize_frecs(place,1);
                }
                _keys[place][_pos] = kmer_repr;
                _info_frecs[place].second++;
            }
        }else{
            _info_frecs[place].first += 2;
            _frecs[place] =  (uint16_t *)realloc(_frecs[place]
                    ,_info_frecs[place].first*sizeof(uint32_t));
            _keys[place] = (size_t*)realloc(_keys[place]
                    ,_info_frecs[place].first*sizeof(size_t));
            _initialize_frecs(place,1);
            _frecs[place][0] = _frec_unique[place];
            _keys[place][0] = _keys_unique[place];
            _keys[place][1] = kmer_repr;
            _uniques_frec[place] = 0;_non_uniques_frec[place] = 1;
            _info_frecs[place].second += 2;
        }
    }
}

Kmer BitMapBig::return_kmer(int rank, int pos){
    if (_non_uniques_frec[rank] == 0)
        return Kmer(_keys_unique[rank]);
    return Kmer(_keys[rank][pos]);
}

void BitMapBig::update(bool choose){
    if (choose) {
        //(&_uniques_frec)->~bit_vector();
        sdsl::util::clear(_uniques_frec);
        sdsl::util::clear(b);
        _b_final = bit_vector(genome_size);
    }else{
        //(&this->b)->~bit_vector();
        b_rank_tmp = rank_support_v<1>(&_b_final);
        _uniques = bit_vector(_cont_local);
        _not_unique = bit_vector(_cont_local);
        /*_keys_not_unique.reserve(_cont_local);
        _not_unique_read_pos.reserve(_cont_local);*/
    }
}

int BitMapBig::check(Kmer kmer,bool standar) {
    if (standar)
        kmer.standardForm();
    size_t kmer_repr = kmer.getRepresentation();
    size_t bucket = kmer.hash() % (size_t)(genome_size);
    if (bucket <= _b_final.size())
        if (_b_final[bucket] == 1) {
            size_t place = b_rank_tmp(bucket);
            if (_not_unique[place] == 0) {
                size_t pos = place-b_rank_uniques(place);
                if (_uniques_keys[pos] == kmer_repr)
                    return pos;
                return -1;
            }
            size_t place2 = b_rank_uniques(place);
            std::function<int(size_t)> _find = [this, place2] (size_t kmer_repr){
                for (int i = 0; i < _info_new_app[place2].second; i++)
                    if (_keys_not_unique[place2][i] == kmer_repr) return i;
                return (int)_info_new_app[place2].second;
            };
            int pos = _find(kmer_repr);
            if (pos != _info_new_app[place2].second)
                return pos;
            else
                return -1;
        }
    return -1;
}

Read_Pos * BitMapBig::get_pos(Kmer, size_t) {
    printf("En algun momento entro aqui\n");
    return nullptr;
}

std::pair<size_t,IndexChunk *> BitMapBig::get_kmer_reads(Kmer kmer, size_t pos){
    size_t place = b_rank_tmp(kmer.hash() % (size_t) (genome_size));
    if (_not_unique[place] == 0) {
        size_t pos_r = b_rank_uniques(place);
        return std::pair<size_t, IndexChunk *>(_unique_read_pos[place - pos_r].size,
                                               _unique_read_pos[place - pos_r].data);
    }else {
        size_t position = b_rank_uniques(place);
        return std::pair<size_t, IndexChunk *>(_not_unique_read_pos[position][pos].size,
                                                 _not_unique_read_pos[position][pos].data);
    }
}
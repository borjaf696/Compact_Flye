#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <vector>
#include "unistd.h"
#include "kmer.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <tuple>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/wt_int.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <cuckoohash_map.hh>
#define NUM_BITS 3
#define NUM_BITS_FRECS 16
#define POSITION_WRONG 99999
using namespace sdsl;

struct IndexChunk
{
    IndexChunk():
            hi(0), low(0) {}
    IndexChunk(const IndexChunk& other):
            hi(other.hi), low(other.low) {}
    IndexChunk(IndexChunk&& other):
            hi(other.hi), low(other.low) {}
    IndexChunk& operator=(const IndexChunk& other)
    {
        hi = other.hi;
        low = other.low;
        return *this;
    }

    size_t get() const
    {
        return ((size_t)hi << 32) + (size_t)low;
    }
    void set(size_t val)
    {
        low = val & ((1ULL << 32) - 1);
        hi = val >> 32;
    }

    uint8_t hi;
    uint32_t low;
} __attribute__((packed));
static_assert(sizeof(IndexChunk) == 5,
              "Unexpected size of IndexChunk structure");

struct ReadPosition
{
    ReadPosition(FastaRecord::Id readId = FastaRecord::ID_NONE,
                 int32_t position = 0):
            readId(readId), position(position) {}
    ReadPosition& operator=(const ReadPosition& rp){
        this->readId = rp.readId;
        this->position = rp.position;
        return *this;
    }
    FastaRecord::Id readId;
    int32_t position;
};

struct Read_Pos
{
    Read_Pos():size(0),data(nullptr){}
    Read_Pos(IndexChunk * rp):size(0),data(rp),frec(0){}
    Read_Pos(IndexChunk * rp,size_t frec):size(0),data(rp), frec(frec){}
    Read_Pos(size_t frec):size(0), frec(frec){}
    uint32_t size, frec;
    IndexChunk* data;
};

class Control{
public:
    //Little
    virtual void insert_alter(Kmer,bool){};
    virtual void update_frec(Kmer) {};
    virtual size_t size() {return 0;};
    virtual size_t size(int) {return 1;};
    /*
     * New
     */
    virtual size_t size_total() {return 1;};
    virtual size_t size_total_chunk(size_t) {return 1;};
    virtual void internal_count(Kmer) {};
    virtual size_t is_unique(size_t, IndexChunk*, size_t pos = 0) {return 1;};
    virtual size_t get_frec_solid(size_t, size_t pos = 99999){return 1;};

    virtual unsigned int get_sumTotal(){return 0;};
    virtual uint32_t getfrec( Kmer, size_t){return -1;};
    virtual uint32_t getfrecIn(size_t){return -1;};
    virtual unsigned int insert(Kmer) {return 0;};
    virtual int locate(Kmer) {return 0;};
    virtual int check(Kmer,bool) {return 0;};
    virtual void setChunckSize(){};
    virtual void mod_frec_pos(Kmer,IndexChunk*){};
    virtual void update_read_pos(Kmer, FastaRecord::Id,size_t, size_t){};
    virtual Read_Pos * get_pos(Kmer,size_t){return nullptr;};
    virtual void set_to_zero(Kmer){};
    //Big
    virtual void insert(Kmer,bool) {};
    virtual void upgrade(bool) {};
    virtual void update(Kmer) {};
    virtual uint32_t get_frec_int(int,int){return 0;}
    virtual Kmer return_kmer(int){return Kmer(0);};
    virtual Kmer return_kmer(int,int){return Kmer(0);};
    virtual void update(bool){};
    virtual bool what(){return false;};
    virtual void erase_frec(uint32_t,std::vector<uint32_t>){};
    virtual void FinalStep(Kmer,size_t){};
    virtual void sortIndexChunks(){};

    virtual std::pair<size_t,IndexChunk *> get_kmer_reads(Kmer,size_t){
        return std::pair<size_t,IndexChunk *>(0, nullptr);
    };
    //Destroyer
    virtual void clear(){};
};

class BitMap: public Control
{
public:
    BitMap();
    BitMap(int,size_t);
    ~BitMap()
    {}
    void mod_frec_pos(Kmer,IndexChunk*);
    void update_read_pos(Kmer,FastaRecord::Id,size_t,size_t);
    void insert_alter(Kmer,bool);
    void update_frec(Kmer);
    void upgrade(bool);
    size_t size();
    Read_Pos * get_pos(Kmer,size_t);
    unsigned int get_sumTotal();
    uint32_t getfrec( Kmer, size_t);
    uint32_t getfrecIn(size_t);
    unsigned int insert(Kmer kmer);
    int locate(Kmer kmer);
    std::pair<size_t,IndexChunk *> get_kmer_reads(Kmer,size_t);
    int check(Kmer,bool);
    void setChunckSize();
    Kmer return_kmer(int);
    void update(bool);
    bool what(){
        return true;
    }
    void set_to_zero(Kmer);
    void sortIndexChunks()
    {
        for (auto rp : read_pos)
        {
            std::sort(rp.data, rp.data + rp.size,
                      [](const IndexChunk& p1, const IndexChunk& p2)
                      {return p1.get() < p2.get();});
        }
    }

    void clear();
private:
    unsigned int _sumTotal = 0;
    void _insert_alter(Kmer kmer, bool){
        size_t kmer_repr = kmer.getRepresentation();
        for (size_t i = 0; i < b_2.size();i++)
            if (!b_2[i][kmer_repr]) {
                b_2[i][kmer_repr] = 1;
                if (i == (b_2.size()-1)) {
                    _sumTotal++;
                    frec.push_back(0);
                }
                return;
            }
        return;
    }
    bit_vector b,_b_final;
    rank_support_v<1> b_rank,b_rank_tmp;
    select_support_mcl<1> b_select;
    std::vector<uint16_t> frec;
    std::vector<Read_Pos> read_pos;
    std::vector<bit_vector> b_2;
    int count = 0;
    size_t genome_size = pow(2,30);
};

class BitMapBig: public Control
{
public:
    BitMapBig(int,size_t);
    void mod_frec_pos(Kmer,IndexChunk*);
    void update_read_pos(Kmer,FastaRecord::Id,size_t,size_t);
    size_t size();
    size_t size(int);
    size_t size_total() {
        return b_rank_tmp(genome_size);
    }
    size_t size_total_chunk(size_t chunck)
    {
        chunck = b_rank_uniques(chunck);
        return _info_new_app[chunck].second;
    }
    void internal_count(Kmer kmer)
    {
        size_t place = b_rank_tmp(kmer.hash() % (size_t) genome_size);
        if (_not_unique[place] == 1)
        {
            _info_new_app[place].first += _block;
            _keys_not_unique[place] = (size_t*)realloc(_keys_not_unique[place]
                    ,_info_new_app[place].first*sizeof(size_t));
            _not_unique_read_pos[place] = (Read_Pos*)realloc(_not_unique_read_pos[place]
                    ,_info_new_app[place].first*sizeof(Read_Pos));
        }
    }
    size_t is_unique(size_t place, IndexChunk* rp, size_t pos)
    {
        if (_not_unique[place] == 0) {
            place -= b_rank_uniques(place);
            _unique_read_pos[place].data = rp;
            return 0;
        }
        place = b_rank_uniques(place);
        _not_unique_read_pos[place][pos].data = rp;
        return place;
    }
    size_t get_frec_solid(size_t place, size_t pos)
    {
        if (_not_unique[place] == 0)
        {
            place -= b_rank_uniques(place);
            return _unique_read_pos[place].frec;
        }
        if (pos != 99999) {
            place = b_rank_uniques(place);
            return _not_unique_read_pos[place][pos].frec;
        }
        return 0;
    }
    unsigned int get_sumTotal();
    void insert(Kmer,bool);
    void upgrade(bool);
    void update(Kmer);
    Read_Pos * get_pos(Kmer,size_t);
    void _insert_2(Kmer);
    int check(Kmer,bool);
    uint32_t get_frec_int(int, int);
    uint32_t getfrec(Kmer,size_t);
    std::pair<size_t,IndexChunk *> get_kmer_reads(Kmer,size_t);
    Kmer return_kmer(int,int);
    void update(bool);
    bool what(){return false;}
    void clear();
    void erase_frec(uint32_t,std::vector<uint32_t>);
    void FinalStep(Kmer,size_t);
    void sortIndexChunks()
    {
        for (auto rp : _unique_read_pos)
        {
            std::sort(rp.data, rp.data + rp.size,
                      [](const IndexChunk& p1, const IndexChunk& p2)
                      {return p1.get() < p2.get();});
        }
        for (size_t i = 0; i < _not_unique_read_pos.size(); ++i)
        {
            for (size_t j = 0; j < _info_new_app[j].second; ++j)
                std::sort(_not_unique_read_pos[i][j].data, _not_unique_read_pos[i][j].data+_not_unique_read_pos[i][j].size,
                          [](const IndexChunk& p1, const IndexChunk& p2)
                          {return p1.get() < p2.get();});
        }
    }
private:
    void _initialize_frecs(size_t place, uint32_t frec){
        size_t limit_l = _info_frecs[place].second, limit_h = _info_frecs[place].first;
        for (size_t i = limit_l; i < limit_h; i++)
            _frecs[place][i] = frec;
    }
    void _insert_1_2(Kmer kmer){
        size_t place = kmer.hash() % (size_t) genome_size;
        if (bloom[place] == (pow(2,NUM_BITS)-1))
            return;
        bloom[place]++;
        /*if (_h_threshold > 1)
            if (b[place] == 1)
                return;
        if (bloom[place] < _h_threshold)
            bloom[place]++;
        else
        {
            _sumsup++;
            _sum_unique++;
            b[place] = 1;
        }*/
        /*for (size_t i = 0; i < _h_threshold; i++) {
            if (!b_2[i][place]) {
                b_2[i][place] = 1;
                if (i == (b_2.size() - 1)) {
                    _sumsup++;
                    _sum_unique++;
                    b[place] = 1;
                }
                return;
            }
        }*/
        return;
    }
    void _update_2(){
        printf("Not supported yet\n");
    }
    unsigned int _sumsup = 0,_sum_unique = 0, _cont_local = 0, _cuenta_interna,_cuckoos_count = 0;
    size_t _h_threshold;
    bool _big_genome = false;
    rank_support_v<1> b_rank,b_rank_tmp,b_rank_uniques;
    bit_vector b,_b_final;
    //New Approach
    int_vector<NUM_BITS> bloom;
    bit_vector _uniques,_not_unique, _uniques_frec, _non_uniques_frec;
    std::vector<size_t> _uniques_keys;
    std::vector<Read_Pos> _unique_read_pos;
    std::vector<size_t *> _keys_not_unique;
    std::vector<Read_Pos *> _not_unique_read_pos;
    std::vector<bit_vector> b_2;
    size_t genome_size = pow(2,30);

    //Test arrays
    std::vector<size_t> _keys_unique;
    std::vector<uint16_t> _frec_unique;
    std::vector<uint16_t *> _frecs;
    std::vector<size_t *> _keys;
    std::vector<std::pair<uint8_t,uint8_t>> _info_frecs, _info_new_app;
    uint8_t _block = 1;

    //Cuckoos incorporation
    cuckoohash_map<Kmer, size_t> _big_case;
};
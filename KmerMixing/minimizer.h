#ifndef MINIMIZER_H
#define MINIMIZER_H
#include <assert.h>
#include <string>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <iomanip>
#include <bitset>
#include <iostream>
#include <fstream>
#include "robin_hood.h"


namespace mzr
{
const uint8_t POSITIVE = 0;
const uint8_t NEGATIVE = 1;

uint8_t char_to_num(char c)
{
    switch (c)
    {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    default:
        return 0;
    }
}
uint8_t char_to_num_complement(char c)
{
    switch (c)
    {
    case 'A':
        return 3;
    case 'C':
        return 2;
    case 'G':
        return 1;
    case 'T':
        return 0;
    default:
        return 0;
    }
}

char num_to_char(uint8_t b)
{
    switch (b)
    {
    case 0:
        return 'A';
    case 1:
        return 'C';
    case 2:
        return 'G';
    case 3:
        return 'T';
    default:
        return 'N';
    }
}


/* void my_assert(bool b, std::string m = "")
{
    if (!b)
    {
        logger.Error("condition failed with message: " + m);
        exit(-1);
    }
} */

const char int_to_char_map[] = {'A','C','G','T'};



std::string unpack_kmer(uint64_t packed_kmer, uint32_t kmer_len);


/* inline char num_to_char(uint8_t b)
{
    if (b < 4)
        return int_to_char_map[b];
    return 'N';
} */

struct Kmer
{
    uint64_t seq;
    uint32_t pos;
    uint8_t sign;
    Kmer(uint64_t s, uint32_t p, uint8_t sg) : seq(s), pos(p), sign(sg) {}
    Kmer(const Kmer& k)
    {
        seq = k.seq;
        pos = k.pos;
        sign = k.sign;
    }
    Kmer() : seq(0), pos(0), sign(POSITIVE){};

    std::string as_string(int8_t k) const
    {        
        return unpack_kmer(seq, k);
    }    
};

struct KmerGenerator
{
    KmerGenerator(std::string const &s, uint8_t k, uint64_t mask, bool lex_low = false) : seq(s), k(k), mask(mask), lxl(lex_low)
    {
        seq_it = seq.begin() + k;
        kmer = Kmer(pack(seq.begin(), seq_it), 0, POSITIVE);

        // replace the kmer with the lexicographically lowest version between it
        // and it's reverse complement
        if (lex_low)
        {
            // std::string rev_seq(s);
            // std::reverse(rev_seq.begin(), rev_seq.end());
            // auto rev_seq_it = rev_seq.end() - k -  1;
            // rev_kmer = Kmer(pack(rev_seq_it + 1, rev_seq.end(), true), 0, NEGATIVE);

            rev_kmer = Kmer(pack_reverse(seq.begin(), seq_it, true), 0, NEGATIVE);
        }
        init = true;
    }

    bool empty() const { return seq_it >= seq.end(); }

/*     std::string kmer_as_string()
    {
        return seq.substr(kmer.pos, k);
    } */

    Kmer get_kmer()
    {
        
        if (!init) // it's an intermediate call, and we need to pack a new base
        {            
            uint32_t pos = seq_it - k - seq.begin() + 1;
            uint64_t n = char_to_num(*seq_it);

            kmer = Kmer(mask & (kmer.seq << 2 | n), pos, POSITIVE);            

            if (lxl)
            {
                uint64_t nc = char_to_num_complement(*seq_it);
                rev_kmer = Kmer(mask & (rev_kmer.seq >> 2 | (nc << (k * 2 - 2))), pos, NEGATIVE);
            }
            seq_it++;
            return (rev_kmer.seq < kmer.seq && lxl) ? rev_kmer : kmer;
        }
        else // first call
        {
            init = false;
            return (rev_kmer.seq < kmer.seq && lxl) ? rev_kmer : kmer;
        }
    }

    Kmer min_kmer_in_window(uint8_t w, std::unordered_set<uint64_t> const &hfk) const
    {
        std::string::const_iterator f = seq_it - w + 1; // move f to the start of the forward window
        uint64_t _kmer = pack(f - k, f);
        uint64_t _rev_kmer = pack_reverse(f - k, f, true);

        // select lowest kmer
        uint64_t min_k = default_kmer.seq;
        uint32_t min_pos = f - seq.begin() - k;
        uint8_t sign = POSITIVE;

        // point to the next unprocessed characters
       
        while (f <= seq_it)
        {
            if (_kmer < min_k && hfk.find(_kmer) == hfk.end())
            {
                min_k = _kmer;
                min_pos = f - seq.begin() - k;
                sign = POSITIVE;
            }
            if (_rev_kmer < min_k && hfk.find(_rev_kmer) == hfk.end())
            {
                min_k = _rev_kmer;
                min_pos = f - seq.begin() - k;
                sign = NEGATIVE;
            }            
            
            _kmer = mask & (_kmer << 2 | char_to_num(*f));
            _rev_kmer = mask & (_rev_kmer >> 2 | char_to_num_complement(*f) << (k * 2 - 2));
            f++;
        }
        
        return {min_k, min_pos, sign};
    }

    uint64_t pack(std::string::const_iterator begin, std::string::const_iterator end, bool complement = false) const
    {
        uint64_t _kmer = 0;
        uint64_t base;
        while (begin < end)
        {
            base = char_to_num(*begin);
            _kmer = (_kmer << 2) | base;
            begin++;
        }
        return mask & (!complement ? _kmer : ~_kmer);
    }
    uint64_t pack_reverse(std::string::const_iterator begin, std::string::const_iterator end, bool complement = false) const
    {
        uint64_t _kmer = 0;
        uint64_t base;
        uint64_t i = 0;
        while (begin < end)
        {
            base = char_to_num(*begin);
            _kmer = _kmer | base << 2 * i;
            i++;
            begin++;
        }
        return mask & (!complement ? _kmer : ~_kmer);
    }

    std::string seq;
    std::string::iterator seq_it;
    uint8_t k;
    uint64_t mask;
    bool lxl, init;
    Kmer kmer, rev_kmer, default_kmer{0xFFFFFFFFFFFFFFFF, 0, POSITIVE};
    
};
}

#endif
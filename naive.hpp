//
// Created by nicco on 7/31/24.
//

//Naive implementations for rank and select as ground-truth for validation

#pragma once

#include "bv.hpp"

class naive : public Bitvector {
public:
    void build_auxiliaries() override {
        // No auxiliary structures needed for naive implementation
    }

    uint64_t space_in_bits() const override {
        return 64 * bits.size();
    }

    bool access(const uint64_t i) const override {
        uint64_t W = bits[i/64];
        return (W >> (i%64)) & 1ULL;
    }

    //Rank: Does a raw linear scan in O(n)
    uint64_t rank_1(uint64_t i) const override {
        uint64_t rank=0;
        uint64_t w = (i/64);
        //Scan linearly over words
        for(uint64_t v=0; v < w; v++) {
            rank += std::popcount(bits[v]);
        }
        //Scan inside the word
        uint64_t W = bits[w];
        uint64_t i1 = i%64;
        for(uint64_t j=0; j<i1; j++) {
            rank += (W >> j) & 1ULL;
        }
        return rank;
    }

    uint64_t rank_0(const uint64_t i) const override {
        return i - rank_1(i);
    }

    /*
     * Select: Does a raw linear scan in O(n).
     * Type=0/1
     */
    static uint64_t select(const int type, const uint64_t k) {
        uint64_t i=0;
        uint64_t c=0;
        uint64_t v=0;
        //Scan linearly over words
        while(c<k) {
            int pop = std::popcount(bits[v]);
            type==1 ? c+=pop : c+=64-pop;
            v++;
            i+=64;
        }
        v--;
        i-=64;
        int pop = std::popcount(bits[v]);
        type==1 ? c-=pop : c-=64-pop;
        uint64_t W = bits[v];
        //Scan inside the word
        for(int j=0; j<WORDSIZE; j++) {
            int bit = (W>>j) & 1ULL;
            type==1 ? c+=bit : c+=1-bit;
            if(c==k) {
                break;
            }
            i++;
        }
        return i;
    }

    uint64_t select_0(int64_t k) const override {
        return select(0,k);
    }

    uint64_t select_1(int64_t k) const override {
        return select(1,k);
    }

}; //class naive
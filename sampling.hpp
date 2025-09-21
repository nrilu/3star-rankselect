//
// Created by nicco on 8/10/24.
//

#pragma once
#include "bv.hpp"

/*
    * Sample every k-th 0 and 1 and store the index of its superblock in S0 or S1_
    * Default is to sample every 16384-th digit or each type
    */
   inline void buildSelect(std::vector<uint32_t> &S0, std::vector<uint32_t> &S1, bool storeIndex) {
       uint64_t s0N = N_zeros / SAMPLEDIST +2;
       uint64_t s1N = N_ones / SAMPLEDIST +2;
       S0 = std::vector<uint32_t>(s0N, 0);
       S1 = std::vector<uint32_t>(s1N, 0);
       int c0 = 0;
       int c1 = 0;
       int s0 = 0;
       int s1 = 0;
       bool first0 = true;
       bool first1 = true;
       int words_in_header = 2;
       int superblocks_in_header = 1;
       int hdsize = 0;
       storeIndex? hdsize = words_in_header : hdsize = superblocks_in_header;
       //Go through whole bitvector
       for(uint64_t w=0; w<N_words; w++) {
           //Save index of the h1-word of the next superblock
           int pop = std::popcount(bits[w]);
           c0 += 64 - pop;
           c1 += pop;
           //Edge-Case that first superblocks are completely void of one of the two digits
           if(first0 && c0>0) {
               S0[s0++] = hdsize*(w/WORDS_IN_L1) + hdsize; // = index of next header = h1
               first0 = false;
           }
           if(first1 && c1>0) {
               S1[s1++] = hdsize*(w/WORDS_IN_L1) + hdsize;
               first1 = false;
           }
           //General case after a digit has been seen the first time
           if(c0>= SAMPLEDIST) {
               S0[s0++] = hdsize*(w/WORDS_IN_L1) + hdsize;
               c0 -= SAMPLEDIST;
           }
           if(c1>= SAMPLEDIST) {
               S1[s1++] = hdsize*(w/WORDS_IN_L1) + hdsize;
               c1 -= SAMPLEDIST;
           }
       }
       //sentinels
       S0[s0] = S0[s0-1];
       S1[s1] = S1[s1-1];
   }

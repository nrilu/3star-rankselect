/*
 *  Converts raw bitvector bits into our uint64_t representation
 *  Takes input either from a textfile or generates it synthetically during runtime
 *
 */

#pragma once

#include "bv.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream> //on some machines needed for isstringstream

namespace file_io {

    using std::cout;
    using std::endl;

inline void checkAlignement(const std::vector<uint64_t> &v) {
    constexpr size_t CACHE_LINE_SIZE = 64;
    uintptr_t address = reinterpret_cast<uintptr_t>(v.data());
    if (address % CACHE_LINE_SIZE == 0) {
        std::cout << "Vector is cache line aligned." << std::endl;
    } else {
        std::cout << "Vector is not cache line aligned." << std::endl;
    }

}


    /*Read the bitvector form file, store in a continous <uint64> vector*/
    inline void readBits(const std::string &path) {

        //Reset internal storate
        bits = std::vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>>(0);
        N_bits=0; N_queries=0; N_ones=0; N_zeros=0; N_words=0; dRead = 0;

        //Read number of queries
        std::ifstream file(path);
        if(!file) {cout << "Error: Could not read file or did not find file" << endl;}
        file >> N_queries;
        file.close();

        //Open file again, now reading bits in buffered chunks (Ca. 2x faster than Character by Character)
        constexpr std::size_t bufferSize = 1 << 14;
        std::vector<char> buffer(bufferSize);
        std::ifstream inFile(path, std::ios::binary);
        if (!inFile) {std::cerr << "Error opening file for reading." << std::endl;}

        //skip first line, have already read it
        inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        //prepare reading bits
        uint64_t n=0; //total bit count
        uint64_t W=0; //current word to fill with bits
        int i=0;      //bit count in current word
        bool EOBITS = false;
        //start reading bits
        while (!EOBITS && (inFile.read(buffer.data(), buffer.size()) || inFile.gcount() > 0)) {
            std::size_t bytesRead = inFile.gcount();
            std::string chunk(buffer.data(), bytesRead);
            for(int j=0; j<bytesRead; j++) {
                if(chunk[j]=='\n') {
                    EOBITS = true;
                    break;
                }
                if(chunk[j]=='1') {
                    W |= 1ULL << i;
                }
                i++;
                n++;
                if(i==64) {
                    //Save the accumulated 64 bits
                    bits.push_back(W);
                    int pop = std::popcount(W);
                    N_ones += pop;
                    N_zeros += 64 - pop;
                    W=0;
                    i=0;
                }
            }
        }
        if(i!=0) {
            //Save the last partial bits that didnt fill a full word
            bits.push_back(W);
            int pop = std::popcount(W);
            N_ones += pop;
            N_zeros += i-pop;
            N_words = 1;
        }
        N_bits = n;
        N_words += N_bits/WORDSIZE;
        //bits.shrink_to_fit();
        inFile.close();
    }

/* Reads the Queries and store them in a vector, to later iterate over the vector
 *
 * Each query is stored as a 64 Bit Integer,
 * with the query type in the highest 3 bits and the query index in the lower 61 bits
 Query type codes:
     access   = 0
     rank 0   = 1
     rank 1   = 2
     select 0 = 3
     select 1 = 4
*/

inline void readQueries(const std::string &path) {
    //Reset storage
    queries = std::vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>>(0);
    //Open file
    std::ifstream inFile(path);
    if(!inFile) {cout << "Error: Could not read file or did not find file" << endl;}
    //Skip query number and bitvector
    inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    //Codes for query types
    uint64_t ACCESS_CODE = 0;
    uint64_t RANK_CODE = 1; //+1 when adressing 1s
    uint64_t SELECT_CODE = 3; //+1 when adressing 1s

    //Read queries
    std::string qString;
    while (std::getline(inFile, qString)) {
        uint64_t qCode = 99;
        uint64_t offset= 0;
        //is access-query
        if(qString[0]=='a') {
            qString.erase(0,6);
            qCode = ACCESS_CODE;
        }
        else {
            //is rank-query
            if(qString[0]=='r') {
                qCode = RANK_CODE;
                offset = 5;
            }
            //is select-query
            else if(qString[0]=='s') {
                qCode = SELECT_CODE;
                offset=7;
            }
            //check whether query is about 0s or 1s
            char digit = qString[offset];
            if(digit=='1') {
                qCode+=1;
            }
            //leave only the index information
            qString.erase(0,offset+2);
        }
        if(qCode!=99) {
            //Read Index
            uint64_t qIndex;
            std::istringstream iss(qString);
            iss >> qIndex;
            //Combine Query Type and Index into single 64-Bit word
            uint64_t query_coded=0;
            query_coded |= qCode << 61;
            query_coded |= qIndex;
            //Save Query
            queries.push_back(query_coded);
        }
    }
    if(queries.size() != N_queries) {cout << "Found different number of queries than stated at beginning" << endl;}
    inFile.close();
}

inline void readFile(const std::string& path) {
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    readBits(path);
    readQueries(path);

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    dRead = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
}

    /* * Write the query results to a file*/
    inline void write_to_file(const std::vector<uint64_t> &results, const std::string &out_path) {
    std::ofstream file;
    file.open(out_path);
    uint64_t k = 0;
    for(auto r : results) {
        file << r;
        k++;
        if(k<results.size()) {
            file << endl;
        }
    }
    file.close();
}

    inline void printWord(int w) {
    uint64_t W = bits[w];
    for(int i=0; i<64; i++) {
        std::cout << ((W>>i) & 1ULL);
    }
}

    inline void printInitialBits(int numRows) {
        for(int w=0; w<numRows; w++) {
            printWord(w);
            std::cout << std::endl;
        }
    }

}

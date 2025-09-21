//
// Created by nicco on 8/11/24.
//

#if LOAD_PASTA

#pragma once


	#include "bv.hpp"
	#include "pasta/bit_vector/bit_vector.hpp"
	#include "pasta/bit_vector/support/flat_rank_select.hpp"
	#include "pasta/bit_vector/support/rank.hpp"


	class pastaflat : public Bitvector {

		pasta::BitVector bv;
		pasta::FlatRankSelect<> rs;
		//The SELECT_SAMPLE_RATE of pasta in set in the file "flat_rank.hpp"

	public:

		void build_auxiliaries() override {
			//std::cout << "pasta FlatRankSelectConfig::SELECT_SAMPLE_RATE = " << pasta::FlatRankSelectConfig::SELECT_SAMPLE_RATE << std::endl;

			bv.resize(N_bits);
			uint64_t i = 0;
			for(uint64_t W : bits) {
				for(int j=0; j<64; j++) {
					bv[i++] = (W >> j) & 1ULL;
				}
			}
			tCopyFinished = std::chrono::steady_clock::now();

			rs = pasta::FlatRankSelect<>(bv);
		}

		uint64_t space_in_bits() const override {
			uint64_t total = ((double)N_bits) * 1.0358;
			return total;
		}

		bool access(const uint64_t i) const override {
			return 0; //not supported
		}

		uint64_t rank_0(const uint64_t i) const override {
			return rs.rank0(i);
		}
		uint64_t rank_1(const uint64_t i) const override {
			return rs.rank1(i);
		}

		uint64_t select_0(int64_t k) const override {
			return rs.select0(k);
		}

		uint64_t select_1(int64_t k) const override {
			return rs.select1(k);
		}


	}; //class pastaflat
#endif
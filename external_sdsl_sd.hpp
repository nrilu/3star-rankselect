//
// Created by nicco on 12/8/24.
//
#if LOAD_SDSL

#pragma once
#include "bv.hpp"
#include <sdsl/bit_vectors.hpp>

class external_sdsl_sd : public Bitvector {

	sdsl::sd_vector<> sd_vec;
	sdsl::sd_vector<>::select_1_type sd_select_1;
	sdsl::sd_vector<>::select_0_type sd_select_0;
	sdsl::sd_vector<>::rank_1_type   sd_rank_1;
	sdsl::sd_vector<>::rank_0_type   sd_rank_0;
public:

	void build_auxiliaries() override {
		cout << "Copying data into sdsl_sd_vector" << endl;

		// auto bv = sdsl::bit_vector(N_bits);
		// uint64_t i = 0;
		// for(uint64_t W : bits) {
		// 	for(int j=0; j<64; j++) {
		// 		bv[i++] = (W >> j) & 1ULL;
		// 	}
		// }

		sdsl::sd_vector_builder builder(N_bits, N_ones);
		uint64_t i = 0;
		for (uint64_t word : bits) {
			for (int j = 0; j < 64; ++j, ++i) {
				if ((word >> j) & 1ULL) {
					builder.set(i);
				}
			}
		}


		tCopyFinished = std::chrono::steady_clock::now();
		// sd_vec = sdsl::sd_vector<>(bv);
		sd_vec = sdsl::sd_vector<>(builder);
		sd_select_1 = sdsl::sd_vector<>::select_1_type(&sd_vec);
		sd_select_0 = sdsl::sd_vector<>::select_0_type(&sd_vec);
		sd_rank_1   = sdsl::sd_vector<>::rank_1_type(&sd_vec);
		sd_rank_0   = sdsl::sd_vector<>::rank_0_type(&sd_vec);

	}

	uint64_t select_1(int64_t k) const override {
		return sd_select_1(k);
	}

	uint64_t select_0(int64_t k) const override {
		return sd_select_0(k);
	}

	uint64_t rank_1(const uint64_t i) const override {
		return sd_rank_1.rank(i);
	}

	uint64_t rank_0(const uint64_t i) const override {
		return sd_rank_0.rank(i);
	}

	bool access(const uint64_t i) const override {return 0;}


	uint64_t space_in_bits() const override {
		int64_t sdvec_bits = 8*sdsl::size_in_bytes(sd_vec);
		int64_t s1_bits    = 8*sdsl::size_in_bytes(sd_select_1);
		int64_t s0_bits    = 8*sdsl::size_in_bytes(sd_select_0);
		int64_t r_bits     = 8*sdsl::size_in_bytes(sd_rank_1);
		int64_t r0_bits    = 8*sdsl::size_in_bytes(sd_rank_0);
		int64_t total	   = sdvec_bits + s1_bits + s0_bits + r_bits + r0_bits;
		if(show_overhead) {
			cout << "sdsl_sd space" << endl;
			print_overhead_line("input", N_bits);
			print_separator();
			print_overhead_line("sd_vec", sdvec_bits);
			print_overhead_line("s1", s1_bits);
			print_overhead_line("s0", s0_bits);
			print_overhead_line("r ", r_bits);
			print_overhead_line("r0", r0_bits);
			print_separator();
			print_overhead_line("total", total);
			print_overhead_line("ovhd", total - N_bits);
			cout << endl;
		}
		return total;
	}

}; //class external_sdsl_sd


#endif

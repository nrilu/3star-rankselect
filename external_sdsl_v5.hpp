//
// Created by nicco on 12/8/24.
//
#if LOAD_SDSL

#pragma once
#include "bv.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v5.hpp>

class external_sdsl_v5 : public Bitvector {

	sdsl::bit_vector bv;;
	sdsl::rank_support_v5<> rank_v5;

public:

	void build_auxiliaries() override {
		bv.bit_resize(N_bits);
		uint64_t i = 0;
		for(uint64_t W : bits) {
			for(int j=0; j<64; j++) {
				bv[i++] = (W >> j) & 1ULL;
			}
		}
		tCopyFinished = std::chrono::steady_clock::now();
		rank_v5 = sdsl::rank_support_v5<>(&bv);
	}

	uint64_t select_1(int64_t k) const override {
		cout << "select_1: not supported by sdsl_rank_v5" << endl;
		return 0;
	}

	uint64_t select_0(int64_t k) const override {
		cout << "select_0: not supported by sdsl_rank_v5" << endl;
		return 0;
	}

	uint64_t rank_1(const uint64_t i) const override {
		return rank_v5.rank(i);
	}

	uint64_t rank_0(const uint64_t i) const override {
		return i - rank_1(i);
	}

	bool access(const uint64_t i) const override {return 0;}

	uint64_t space_in_bits() const override {
		int64_t bv_bits = 8*sdsl::size_in_bytes(bv);
		int64_t r_bits  = 8*sdsl::size_in_bytes(rank_v5);
		int64_t total   = bv_bits + r_bits;
		if(show_overhead) {
			cout << "sdsl_v5 space" << endl;
			print_overhead_line("input:", N_bits);
			print_separator();
			print_overhead_line("bv", bv_bits);
			print_overhead_line("rank_v5", r_bits);
			print_separator();
			print_overhead_line("total", total);
			print_overhead_line("ovhd", total-N_bits);
			cout << endl;
			//string t = get_unixtimestamp();
			//sdsl::write_structure<sdsl::HTML_FORMAT>(bvs1, "plotting/visuals/sdsl_htmls/sdsl_default_bvs1"+t+".html");
			//sdsl::write_structure<sdsl::HTML_FORMAT>(bvs0, "plotting/visuals/sdsl_htmls/sdsl_default_bvs0"+t+".html");
		}
		return total;
	}

}; //class external_sdsl_v5

#endif

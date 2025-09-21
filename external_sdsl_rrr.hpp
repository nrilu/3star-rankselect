//
// Created by nicco on 12/8/24.
//
#if LOAD_SDSL

#pragma once
#include "bv.hpp"
#include <sdsl/bit_vectors.hpp>

class external_sdsl_rrr : public Bitvector {

	sdsl::rrr_vector<RRR_BLOCKSIZE> rrrb;
	sdsl::rrr_vector<RRR_BLOCKSIZE>::select_1_type rrr_select_1;
	sdsl::rrr_vector<RRR_BLOCKSIZE>::select_0_type rrr_select_0;
	sdsl::rrr_vector<RRR_BLOCKSIZE>::rank_1_type rrr_rank_1;

public:

	void build_auxiliaries() override {
		auto bv = sdsl::bit_vector(N_bits);
		uint64_t i = 0;
		for(uint64_t W : bits) {
			for(int j=0; j<64; j++) {
				bv[i++] = (W >> j) & 1ULL;
			}
		}
		tCopyFinished = std::chrono::steady_clock::now();

		rrrb = sdsl::rrr_vector<RRR_BLOCKSIZE>(bv);
		rrr_select_1 = sdsl::rrr_vector<RRR_BLOCKSIZE>::select_1_type(&rrrb);
		rrr_select_0 = sdsl::rrr_vector<RRR_BLOCKSIZE>::select_0_type(&rrrb);
		rrr_rank_1   = sdsl::rrr_vector<RRR_BLOCKSIZE>::rank_1_type(&rrrb);

	}

	uint64_t select_1(int64_t k) const override {
		return rrr_select_1(k);
	}

	uint64_t select_0(int64_t k) const override {
		return rrr_select_0(k);
	}

	uint64_t rank_1(const uint64_t i) const override {
		return rrr_rank_1.rank(i);
	}

	uint64_t rank_0(const uint64_t i) const override {
		return i-rank_1(i);
	}

	bool access(const uint64_t i) const override {return 0;}



	uint64_t space_in_bits() const override {
		int64_t rrr_bits = 8*sdsl::size_in_bytes(rrrb);
		int64_t s1_bits  = 8*sdsl::size_in_bytes(rrr_select_1);
		int64_t s0_bits  = 8*sdsl::size_in_bytes(rrr_select_0);
		int64_t r_bits   = 8*sdsl::size_in_bytes(rrr_rank_1);
		int64_t total    = rrr_bits + s1_bits + s0_bits + r_bits;

		if(show_overhead) {
			print_overhead_line("input",N_bits);
			print_separator();
			print_overhead_line("rrr-bv", rrr_bits);
			print_overhead_line("s1", s1_bits);
			print_overhead_line("s0", s0_bits);
			print_overhead_line("r ", r_bits);
			print_separator();
			print_overhead_line("total", total);
			print_overhead_line("ovhd", total - N_bits);
			cout << endl;
			//string t = get_unixtimestamp();
			//sdsl::write_structure<sdsl::HTML_FORMAT>(rrrb, "plotting/visuals/sdsl_htmls/sdsl_rrrb"+t+".html");
		}
		return total;
	}

}; //class external_sdsl_rrr

#endif

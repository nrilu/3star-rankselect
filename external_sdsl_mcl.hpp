//
// Created by nicco on 12/8/24.
//
#if LOAD_SDSL

#pragma once
#include "bv.hpp"
#include <sdsl/bit_vectors.hpp>

class external_sdsl_mcl : public Bitvector {

public:
	sdsl::bit_vector bv;
	sdsl::bit_vector::select_1_type mcl_select_1;
	sdsl::bit_vector::select_0_type mcl_select_0;
	//namespace corresponds to Munroe/Clarks mcl algorithm

	void build_auxiliaries() override {
		bv.bit_resize(N_bits);
		uint64_t i = 0;
		for(uint64_t W : bits) {
			for(int j=0; j<64; j++) {
				bv[i++] = (W >> j) & 1ULL;
			}
		}
		tCopyFinished = std::chrono::steady_clock::now();

		mcl_select_1 = sdsl::bit_vector::select_1_type(&bv);
		mcl_select_0 = sdsl::bit_vector::select_0_type(&bv);
	}

	uint64_t select_1(int64_t k) const override {
		return mcl_select_1(k);
	}

	uint64_t select_0(int64_t k) const override {
		return mcl_select_0(k);
	}

	uint64_t rank_1(const uint64_t i) const override {
		cout << "rank_1: not supported by sdsl_mcl" << endl;
		return 0;
	}

	uint64_t rank_0(const uint64_t i) const override {
		cout << "rank_0: not supported by sdsl_mcl" << endl;
		return 0;
	}

	bool access(const uint64_t i) const override {return 0;}



	uint64_t space_in_bits() const override {
		int64_t bv_bits = 8*sdsl::size_in_bytes(bv);
		int64_t s1_bits = 8*sdsl::size_in_bytes(mcl_select_1);
		int64_t s0_bits = 8*sdsl::size_in_bytes(mcl_select_0);
		int64_t total   = bv_bits + s1_bits + s0_bits;
		if(show_overhead) {
			cout << "   sdsl_mcl space" << endl;
			print_overhead_line("input",   N_bits);
			print_separator();
			print_overhead_line("bv",       bv_bits);
			print_overhead_line("s1", s1_bits);
			print_overhead_line("s0", s0_bits);
			print_separator();
			print_overhead_line("total", total);
			print_overhead_line("ovhd", total - N_bits);
			cout << endl;
			//string t = get_unixtimestamp();
			//sdsl::write_structure<sdsl::HTML_FORMAT>(bvs1, "plotting/visuals/sdsl_htmls/sdsl_default_bvs1"+t+".html");
			//sdsl::write_structure<sdsl::HTML_FORMAT>(bvs0, "plotting/visuals/sdsl_htmls/sdsl_default_bvs0"+t+".html");
		}
		return total;
	}

}; //class external_sdsl_mcl


#endif

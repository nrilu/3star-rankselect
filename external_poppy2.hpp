//
// Created by nicco on 12/8/24.
//
#if LOAD_POPPY

#pragma once
#include "bv.hpp"
#include <rankselect2/bitmap.h>
#include <rankselect2/shared.h>

class external_poppy2 : public Bitvector {

	std::unique_ptr<BitmapPoppy> rs_ptr;

public:

	void build_auxiliaries() override {
		cout << "Starting Poppy2 Bit reversing" << endl;
		synth::create_64bit_reversed_bitvector();
		tCopyFinished = std::chrono::steady_clock::now();
		rs_ptr = std::make_unique<BitmapPoppy>(bits_64bit_reversed.data(), N_bits);
	}

	uint64_t select_1(int64_t k) const override {
		return rs_ptr->select(k);
	}

	uint64_t select_0(int64_t k) const override {
		cout << "select_0: not supported by poppy2" << endl;
		return 0;
	}

	uint64_t rank_1(const uint64_t i) const override {
		return rs_ptr->rank(i);
	}

	uint64_t rank_0(const uint64_t i) const override {
		return i-rank_1(i);
	}

	bool access(const uint64_t i) const override {return 0;}


	uint64_t space_in_bits() const override {
		int64_t bv_bits   = N_bits;				//from paper & their code
		int64_t L1L2_bits = 0.03125 * N_bits;   //from paper & their code
		int64_t s1_bits   = 32 * N_ones / 8192; //from paper & their code
		int64_t total = bv_bits + L1L2_bits + s1_bits;
		if(show_overhead) {
			cout << "   poppy2 space" << endl;
			print_overhead_line("input",N_bits);
			print_separator();
			print_overhead_line("bv",	bv_bits);
			print_overhead_line("L1+L2",L1L2_bits);
			print_overhead_line("s1",s1_bits);
			print_separator();
			print_overhead_line("total", total);
			print_overhead_line("ovhd", total - N_bits);
		}
		return total;
	}

}; //class external_poppy2


#endif

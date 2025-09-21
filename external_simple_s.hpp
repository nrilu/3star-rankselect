//
// Created by nicco on 12/8/24.
//
#if LOAD_SUX

#pragma once
#include "bv.hpp"
#include <sux/bits/SimpleSelect.hpp>

template<int SIMPLE_SELECT_PARAM>
class external_simple_select : public Bitvector {

	// Store pointer to avoid copy assignment issues
	std::unique_ptr<sux::bits::SimpleSelect<>> ss_ptr;

public:

	void build_auxiliaries() override {
		tCopyFinished = std::chrono::steady_clock::now();
		ss_ptr = std::make_unique<sux::bits::SimpleSelect<>>(bits.data(), N_bits, SIMPLE_SELECT_PARAM);
	}

	uint64_t select_1(int64_t k) const override {
		return ss_ptr->select(k-1);
	}

	uint64_t select_0(int64_t k) const override {
		cout << "select_0: not supported by simple_select" << endl;
		return 0;
	}

	uint64_t rank_1(const uint64_t i) const override {
		cout << "rank_1: not supported by simple_select" << endl;
		return 0;
	}

	uint64_t rank_0(const uint64_t i) const override {
		cout << "rank_0: not supported by simple_select" << endl;
		return 0;
	}

	bool access(const uint64_t i) const override {return 0;}


	uint64_t space_in_bits() const override {
		int64_t bv_bits = N_bits;
		int64_t ss_bits = ss_ptr->bitCount();
		int64_t total = bv_bits + ss_bits;
		if(show_overhead) {
			cout << "   simple_select space"<<endl;
			print_overhead_line("input: ", N_bits);
			print_overhead_line("bv: ", bv_bits);
			print_overhead_line("ss: ", ss_bits);
			print_overhead_line("ovhd sum:", total - N_bits);
		}
		return total;
	}

}; //class external_simple_select

#endif

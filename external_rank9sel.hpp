//
// Created by nicco on 12/8/24.
//
#if LOAD_SUX

#pragma once
#include "bv.hpp"
#include <sux/bits/Rank9Sel.hpp>

class external_rank9sel : public Bitvector {

	// Store as pointer to avoid copy assignment issues
	std::unique_ptr<sux::bits::Rank9Sel<>> rs_ptr;

public:

	void build_auxiliaries() override {
		tCopyFinished = std::chrono::steady_clock::now();
		rs_ptr = std::make_unique<sux::bits::Rank9Sel<>>(bits.data(), N_bits);
	}

	uint64_t select_1(int64_t k) const override {
		return rs_ptr->select(k-1);
	}

	uint64_t select_0(int64_t k) const override {
		cout << "select_0: not supported by rank9sel" << endl;
		return 0;
	}

	uint64_t rank_1(const uint64_t i) const override {
		return rs_ptr->rank(i);
	}

	uint64_t rank_0(const uint64_t i) const override {
		return i - rank_1(i);
	}

	bool access(const uint64_t i) const override {return 0;}

	uint64_t space_in_bits() const override {
		int64_t bv_bits    = N_bits;
		int64_t rank9_bits = rs_ptr->bitCount();
		int64_t total 	   = bv_bits + rank9_bits;
		if(show_overhead) {
			cout << "Rank9Sel space"<<endl;
			print_overhead_line("input: ", N_bits);
			print_overhead_line("bv: ", N_bits);
			print_overhead_line("s0,s1,r:", rank9_bits);
			print_overhead_line("ovhd sum:", total - N_bits);
		}
		return total;
	}

}; //class external_rank9sel

#endif

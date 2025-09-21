//// Created by nicco on 3/16/25.
//

/**
 *  Our 3* rank-select implementation.
 *  For building the sample-tree, see also Tree3.hpp.
 */

//
#pragma once
#include "bv.hpp"
#include "alignedallocator.hpp"
#include "Tree3.hpp"
#include <cstring> //for memcpy

#define DEBUGINFO_L2_SUMMARY    0
#define DEBUGINFO_L2_SUMM_BUILD 0
#define DEBUGINFO_SB_SELECT		0
#define DEBUGINFO_SB_BUILDING   0
#define DEBUGINFO_T3_WALK       0
#define DEBUGINFO_T3_SCANLENGTH 0
#define DOUBLECHECK_T3_VANILLA  0
#define ALERT_T3_ALPHA_BREAKING 0
#define DEBUGINFO_RANK          0

//#define COUNT_HIT_TYPES 0

/**
 * The Superblock Structure
 *		[L1    sparse_w0    dense_blocks   d0...d7    c0...c7  b0...b7   A2...A7    alldense  A1    minflag  A0    sparse_w0]
 *       43    21           32             8*12       8*12     8*12      5*16       1         15    1        15    16
 *
 *    L1: Number of ones before the superblock
 *    sparse_w0: Start of the sparse representation of this superblock.
 *    dense_blocks: Number of dense blocks before this superblock.
 *    A,b,c,d: L0-values
 *    A: Cumulative count of min-bits in block nr. 4k+3    ( 3, 7,11,15,19,23,27),
 *    b: Relative count of the min-bits in block nr. 4k    ( 0, 4, 8,12,16,20,24,28)
 *    c: Relative count of the min-bits in block nr. 4k+1  ( 1, 5, 9,13,17,21,25,29)
 *    d: Relative count of the min-bits in block nr. 4k+2  ( 2, 6,10,14,18,22,26,30)
 *
 */

using namespace std;



class m3 final : public Bitvector {

	static constexpr bool show_packing_situation = false;
	static constexpr bool show_array_packing = false;

	#if CUSTOM_L0==2048 || CUSTOM_L0==0
	static constexpr uint64_t L0_SIZEz		= 2048;
	static constexpr uint64_t Z_VALUE_SIZE  = 11;  //
	static constexpr uint64_t MAX_SPARSE    = 186; //	2048/11
	#elif CUSTOM_L0==1024
	static constexpr uint64_t L0_SIZEz		= 1024;
	static constexpr uint64_t Z_VALUE_SIZE  = 10;
	static constexpr uint64_t MAX_SPARSE    = 102; //	1024/10
	#elif CUSTOM_L0==512
	static constexpr uint64_t L0_SIZEz		= 512;
	static constexpr uint64_t Z_VALUE_SIZE  = 9;
	static constexpr uint64_t MAX_SPARSE    = 56; //	512/9
	#else
	cout << "Error: m3 only supports L0= 512, 1024 and 2048 (default 0 maps to 2048)"
	#endif


	//L1 Superblock
	static constexpr uint64_t L0_IN_L1z	= 32;
	static constexpr uint64_t L1_SIZEz	= L0_IN_L1z * L0_SIZEz;		//65536
	static constexpr uint64_t W_IN_L0z	= L0_SIZEz / WORDSIZE;		//   32
	static constexpr uint64_t W_IN_L1z	= W_IN_L0z * L0_IN_L1z; // 1024
	static constexpr uint64_t WORDS_IN_L1HEADER  = 8;
	static constexpr uint64_t NUM_SENTINEL_SUPERBLOCKS= 5; //such that we also cover more than one extra summary-level

	static constexpr uint64_t L1_POS	= 7;
	static constexpr uint64_t L1_BITSIZE= 43;
	static constexpr uint64_t L1_SHIFT	= 64 - L1_BITSIZE;

	static constexpr uint64_t S_W0_LOWSIZE	= 16;
	static constexpr uint64_t S_W0_LOWMASK	= (1ULL << S_W0_LOWSIZE) - 1;
	static constexpr uint64_t S_W0_HIGHSIZE = 64 - L1_BITSIZE; //21
	static constexpr uint64_t S_W0_TOTALSIZE = S_W0_LOWSIZE + S_W0_HIGHSIZE;
	static constexpr uint64_t S_W0_TOTALMASK = (1ULL << S_W0_TOTALSIZE) - 1;
	static constexpr uint64_t S_W0_LOW_POS = 0; //low part is in word 0 (next to the A's)
	static constexpr uint64_t S_W0_HIGH_POS = 7; //hight part is in word 7 (next to L1)

	static constexpr uint64_t DENSE_BLOCK_COUNT_POS = 6;
	static constexpr uint64_t DENSE_BLOCK_COUNT_BITSIZE = 32;
	static constexpr uint64_t DENSE_BLOCK_COUNT_SHIFT = 64 - DENSE_BLOCK_COUNT_BITSIZE;

	static constexpr uint64_t L0d_POS = 5;
	static constexpr uint64_t L0c_POS = 3;
	static constexpr uint64_t L0b_POS = 2;
	static constexpr uint64_t L0A_POS = 0;

	static constexpr uint64_t MINFLAG_POS	= 0;
	static constexpr uint64_t MINFLAG_SHIFT	= 31;
	// [										   ...  minflag(1)   A(15)  w0(16)]		//16+15=31

	static constexpr uint64_t ALLDENSE_W  = 0;
	static constexpr uint64_t ALLDENSE_SHIFT = 47;
	// [					  ... dense_flag(1)  A(15)  minflag(1)   A(15)  w0(16)]		// 16 + 15 + 1 + 15 = 47

	static constexpr uint64_t ALLSPARSE_W  = 0;
	static constexpr uint64_t ALLSPARSE_SHIFT = 63;
	// [  ... sparseflag(1) A(15) dense_flag(1)  A(15)  minflag(1)   A(15)  w0(16)]		// 16 + 15 + 1 + 15 + 1 + 15 = 63

	static constexpr uint64_t LAST_L0_IS_SPARSE_W  = 1;
	static constexpr uint64_t LAST_L0_IS_SPARSE_SHIFT = 15;

	//L2 Summary level
	static constexpr uint64_t L2_BITSIZE = 47;
	static constexpr uint64_t L2_SHIFT   = 64 - L2_BITSIZE;
	static constexpr uint64_t L1_IN_L2   = 16;
	static constexpr uint64_t WORDS_IN_L2_SYMMLEVEL_HEADER = 8;
	static constexpr uint64_t L2_POS     = 0;
	static constexpr uint64_t L2_SIZE    = L1_IN_L2 * L1_SIZEz;
	static constexpr uint64_t W_IN_L2    = L1_IN_L2 * W_IN_L1z;
	// static constexpr uint64_t L1B_SPLIT_BITS = 12;
	// static constexpr uint64_t L1B_SHIFT  = 32 - L1B_SPLIT_BITS; //20
	// static constexpr uint64_t L1B_POS = WORDS_IN_L2_SYMMLEVEL_HEADER - 2;
	static constexpr uint64_t NUM_SENTINEL_SUMMARY_BLOCKS= 2;


	typedef std::vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>> aligned_vector;
	aligned_vector L1L0;		//The Superblock Headers (containing 1xL1, 31xL0 and some more meta-data per superblock)
	aligned_vector dense;	    //Exact copy of the input bitvector, reduced to the parts where L0-blocks are dense
	aligned_vector sparse;		//A delta-encoded representation of the remaining sparse L0-blocks bits
	aligned_vector S0, S1;		//naive sampling points when not using Tree
	aligned_vector L2L1;		//Optional summary level above the superblocks

	//runtime variables tracking sparse/dense situation
	uint64_t running_dense_block_count = 0; //count of sparse l0 blocks
	uint64_t running_sparse_w0 = 0;	   //next free word in the s_row for the next superblock
	uint64_t running_zcount_of_block = 0;  //count of the z-entries in s_row for the current superblock

	uint64_t superblocks_counting_0s=0;
	uint64_t superblocks_counting_1s=0;

	//Kick w0 from lowest slot when loading A. Also clip A to only 15 bits each
	__m128i const remove_clutter_around_A =_mm_set_epi16(0x7FFF,0x7FFF,0x7FFF,0x7FFF,0x7FFF,0x7FFF,0x7FFF,    0);

	//Unpacking packed 12-bit entries into 16bit slots         x  x        x  x        x  x        x  x
	__m128i const unpack12to16_shuff_low = _mm_set_epi8(-1,-1,10, 9,-1,-1, 7, 6,-1,-1, 4, 3,-1,-1, 1, 0);
	//                                                  15,14,13,12,11,10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
	//                                                   x  x        x  x        x  x        x  x
	__m128i const unpack12to16_shuff_hi  = _mm_set_epi8(11,10,-1,-1, 8, 7,-1,-1, 5, 4,-1,-1, 2, 1,-1,-1);
	//                                                  15,14,13,12,11,10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
	__m128i const keepLower12_epi16  = _mm_set1_epi16(0x0FFF);

	//Complementing to the opposite bit type
	__m128i const MAX_A = _mm_set_epi16(28*L0_SIZEz, 24*L0_SIZEz, 20*L0_SIZEz, 16*L0_SIZEz, 12*L0_SIZEz, 8*L0_SIZEz, 4*L0_SIZEz,		  0);
	__m128i const MAX_B = _mm_set_epi16(29*L0_SIZEz, 25*L0_SIZEz, 21*L0_SIZEz, 17*L0_SIZEz, 13*L0_SIZEz, 9*L0_SIZEz, 5*L0_SIZEz, 1*L0_SIZEz);
	__m128i const MAX_C = _mm_set_epi16(30*L0_SIZEz, 26*L0_SIZEz, 22*L0_SIZEz, 18*L0_SIZEz, 14*L0_SIZEz,10*L0_SIZEz, 6*L0_SIZEz, 2*L0_SIZEz);
	__m128i const MAX_D = _mm_set_epi16(31*L0_SIZEz, 27*L0_SIZEz, 23*L0_SIZEz, 19*L0_SIZEz, 15*L0_SIZEz,11*L0_SIZEz, 7*L0_SIZEz, 3*L0_SIZEz);

	//Maximum number of bits such that a pair is considered sparse
	__m256i const range_shift_32768_epi16= _mm256_set1_epi16(32768);
	__m128i const sparse_limit_si128 = _mm_set1_epi16(MAX_SPARSE);
	__m256i const sparse_limit_si256 = _mm256_set1_epi16(MAX_SPARSE);
	__m256i const sparse_limit_plus1_si256 = _mm256_set1_epi16(MAX_SPARSE+1); //need +1 for greater-than comparison

	//The order in which blocks end up in packus_epi16
	__m256i const block_order_in_packus = _mm256_set_epi8(31, 29, 27, 25, 23, 21, 19, 17,
														  30, 28, 26, 24, 22, 20, 18, 16,
														  15, 13, 11,  9,  7,  5,  3,  1,
														  14, 12, 10,  8,  6,  4,  2,  0);


	//Keep only pairs for s_count summation that lie befor our pair
	__m128i const keep_evnpairs_indices = _mm_set_epi16(14,12,10,8,6,4,2,0);
	__m128i const keep_oddpairs_indices = _mm_set_epi16(15,13,11,9,7,5,3,1);
	__m128i const mm_L0size_epi16 = _mm_set1_epi16(L0_SIZEz);


	__m256i const summary_mask_hide_L2 = _mm256_set_epi64x(-1,-1,-1,0x1FFFF); //hide L2 in the higher 47 bits of the last word, only keep the lower 17 bits


	//Tree structure
	static constexpr uint64_t SAMPLEDIST_EXP_LVL0 = TREE_EXP_LVL0;
	static constexpr uint64_t SAMPLEDIST_EXP_LVL1 = SAMPLEDIST_EXP_LVL0 - 1*TREE_DIVISION_EXP;
	static constexpr uint64_t SAMPLEDIST_EXP_LVL2 = SAMPLEDIST_EXP_LVL0 - 2*TREE_DIVISION_EXP;
	static constexpr uint64_t SAMPLEDIST_EXP_LVL3 = SAMPLEDIST_EXP_LVL0 - 3*TREE_DIVISION_EXP;
	static constexpr uint64_t SAMPLEDIST_EXP_LVL4 = SAMPLEDIST_EXP_LVL0 - 4*TREE_DIVISION_EXP;
	// static constexpr uint64_t SAMPLEDIST_EXP_LVL4 = 1;

	static constexpr uint64_t MODULO_SAMPLEDIST_LVL0 = (1ULL << SAMPLEDIST_EXP_LVL0)-1;
	static constexpr uint64_t MODULO_SAMPLEDIST_LVL1 = (1ULL << SAMPLEDIST_EXP_LVL1)-1;
	static constexpr uint64_t MODULO_SAMPLEDIST_LVL2 = (1ULL << SAMPLEDIST_EXP_LVL2)-1;
	static constexpr uint64_t MODULO_SAMPLEDIST_LVL3 = (1ULL << SAMPLEDIST_EXP_LVL3)-1;
	static constexpr uint64_t MODULO_SAMPLEDIST_LVL4 = (1ULL << SAMPLEDIST_EXP_LVL4)-1;

	static constexpr uint64_t FLAG_SIZE = 2; //tree type flags need 2 bits of space:
	static constexpr uint64_t SENTINEL_FLAG = 0; // 00
	static constexpr uint64_t L1_FLAG		= 1; // 01
	static constexpr uint64_t EXPLICIT_FLAG = 2; // 10
	static constexpr uint64_t NEXTLEVEL_FLAG= 3; // 11
	static constexpr uint64_t TREE_MASKVALUE = (1ULL << 30)-1; //read entry from 32-bit Tree, masking the initial 2 flag bits
	static constexpr uint64_t TREE_FLAGSHIFT = 30;

	static constexpr uint64_t SAMPLEDIST_LVL0 = 1ULL << SAMPLEDIST_EXP_LVL0;
	static constexpr uint64_t SAMPLEDIST_LVL4 = 1ULL << SAMPLEDIST_EXP_LVL4;
	static constexpr uint64_t TREE_DIVISION_COUNT = 1ULL << TREE_DIVISION_EXP;
	//known common factors for each entry type in the Tree. Can be divided out for storing, reducing space, and then multiplied again during decoding
	static constexpr uint64_t COMPRESSION_FACTOR[4] = {0, WORDS_IN_L1HEADER, SAMPLEDIST_LVL4, TREE_DIVISION_COUNT};

	static constexpr uint64_t MAX_LEVEL = 4; //Levels 0-4

	//Tree metadata
	typedef std::vector<uint32_t, AlignedAllocator<uint32_t,ALIGNMENT>> aligned_32vector;
	aligned_32vector B0, B1, B2, B3, B4; //tracking 1s
	aligned_32vector C0, C1, C2, C3, C4; //tracking 0s
	std::array<aligned_32vector*, 6> B, C; //Combined access to B0..B4 and C0..C4 to pass them to Tree-building

	aligned_vector EXPL;
	std::array<aligned_vector*, 1> EXPL_wrapper; //pointer-like workaround to pass EXPL reference

	static constexpr string spacer="				 ";

	Tree3 T3;



public:
	[[nodiscard]] uint64_t rank_0(const uint64_t i) const override {return i - rank_1(i);}
	[[nodiscard]] bool access(const uint64_t i) const override {return 0;}

	static constexpr uint64_t MAX_SPARSE_PUBLIC = MAX_SPARSE;

	uint64_t space_in_bits() const override {
		int64_t summary_L2_bits    = 64*L2L1.size();
		int64_t summary_L1_bits    = 64*L1L0.size();
		int64_t dense_bits		   = 64*dense.size();
		int64_t sparse_bits		   = 64*sparse.size();
		#if TREE_TYPE==3
		int64_t tree1_bits = T3.get_packed_size_in_bits();
		int64_t tree0_bits = 0; //potential extension: tree also for 0-bits ?
		int64_t expl_bits  = 0;
		#elif TREE_TYPE==5
		int64_t tree1_bits = 32*(B0.size() + B1.size() + B2.size() + B3.size() + B4.size());
		int64_t tree0_bits = 32*(C0.size() + C1.size() + C2.size() + C3.size() + C4.size());
		int64_t expl_bits  = 64*EXPL.size();
		#else
		cout << "Warning, m3 only supports TREE_TYPE={3,5}" << endl;
		#endif

		int64_t total = summary_L1_bits;
				total += summary_L2_bits;
				total += dense_bits;
				total += sparse_bits;
				total += tree1_bits;
				total += tree0_bits;
				total += expl_bits;

		if(show_overhead) {
			cout << "    m3 space" << endl;
			print_overhead_line("input", N_bits);
			print_separator();
			print_overhead_line("1-Tree", tree1_bits);
			print_overhead_line("0-Tree", tree0_bits);
			print_overhead_line("EXPL", expl_bits);
			print_overhead_line("L2+L1", summary_L2_bits);
			print_overhead_line("L1+L0", summary_L1_bits);
			print_overhead_line("dense arr.", dense_bits);
			print_overhead_line("sparse arr.", sparse_bits);
			print_separator();
			print_overhead_line("total", total);
			print_overhead_line("ovhd", total - N_bits);
			cout << endl;
		}
		return total;
	}
	void build_auxiliaries() override {
		tCopyFinished = std::chrono::steady_clock::now(); //no initial copying needed for our own algos (needed for sdsl-vectors)

		#if TREE_TYPE==3

		uint64_t L_exp;
		if (SUMMARY_LEVELS==2) {
			L_exp = _tzcnt_u64(L2_SIZE);
		}
		else {
			L_exp = _tzcnt_u64(L1_SIZEz);
		}
        T3.build_Tree3(L_exp);

		#else
		cout << "Error: m3 needs TREE_TYPE={3}" << endl;
		#endif

		uint64_t n_superblocks = build_L1_L0_array();

		if (SUMMARY_LEVELS==2) {
			build_L2_L1_summary_level(n_superblocks);
		}
		double ovhd = T3.get_packed_size_in_bits()/(double)N_bits;
		printf("													m3: summ-levels=%i, ALPHA=%i, tree-strategy=%i, a_exp=%lu  b_exp=%lu ovhd=%.10f rel_to_theo=%.10f size=%.10f MB,  bits=%lu  \n",
			SUMMARY_LEVELS, ALPHA,TREE3_STRATEGY, T3.a_exp, T3.b_exp, ovhd ,  ovhd/treestats.theo_worstcase_size_factor, T3.get_overhead_in_MB(), T3.get_packed_size_in_bits());
	}

	static bool get_SB_minority_type(uint64_t const l1n) {
		int c=0;
		uint64_t w_start = l1n*W_IN_L1z;
		//only read the first 31 of the 32 L0-blocks because we dont represent the last 32th block
		uint64_t w_end	 = w_start + 31*W_IN_L0z;
		for(uint64_t w = w_start; w < w_end && w<N_words; w++) {
			c += std::popcount(bits[w]);
		}
		//Technically we dont measure exactly the minority of the full superblock,
		//but rather whether there are less than 1024 1-bits in the first 31 blocks, which is already sufficient for our representation
		//(and which biases the flag for random 50:50 distribution to be almost always 1, helping the branche predictor)
		return (c < L1_SIZEz/2);
	}


	//Store L1 count of superblock
	void store_L1(uint64_t const L1, uint64_t const l1n) {
		uint64_t l1 = l1n * WORDS_IN_L1HEADER;
		L1L0[l1 + L1_POS] |= L1 << L1_SHIFT;
	}

	//store how many sparse pairs have come before this superblock
	void store_dense_blocks_before(uint64_t const l1n) {
		uint64_t l1 = l1n * WORDS_IN_L1HEADER;
		L1L0[l1 + DENSE_BLOCK_COUNT_POS] |= running_dense_block_count << DENSE_BLOCK_COUNT_SHIFT;
	}

	//store the word in the s_row where the sparse entries for this superblock start
	void store_sparse_w0(uint64_t l1n) {
		uint64_t l1 = l1n * WORDS_IN_L1HEADER;
		uint64_t low_w0  = running_sparse_w0 & S_W0_LOWMASK;
		uint64_t high_w0 = running_sparse_w0 >> S_W0_LOWSIZE;
		L1L0[l1 + S_W0_LOW_POS] |= low_w0;
		L1L0[l1 + S_W0_HIGH_POS]|= high_w0;
	}

	//Store whether the 1s or the 0s are in the minority in this superblock
	void store_minflag(uint64_t const m, uint64_t const l1n) {
		uint64_t l1 = l1n * WORDS_IN_L1HEADER;
		L1L0[l1 + MINFLAG_POS] |= m << MINFLAG_SHIFT;
	}

	void set_alldense_flag(uint64_t flag, uint64_t const l1n) {
		uint64_t l1 = l1n * WORDS_IN_L1HEADER;
		if (flag==1) {
			L1L0[l1 + ALLDENSE_W] |= 1ULL << ALLDENSE_SHIFT; //set bit
		} else {
			L1L0[l1 + ALLDENSE_W] &= ~(1ULL << ALLDENSE_SHIFT);  // Clear bit
		}
	}

	void set_allsparse_flag(uint64_t flag, uint64_t const l1n) {
		uint64_t l1 = l1n * WORDS_IN_L1HEADER;
		if (flag==1) {
			L1L0[l1 + ALLSPARSE_W] |= 1ULL << ALLSPARSE_SHIFT; //set bit
		} else {
			L1L0[l1 + ALLSPARSE_W] &= ~(1ULL << ALLSPARSE_SHIFT);  // Clear bit
		}
	}

	void store_last_L0_is_sparse(uint64_t const l1n) {
		uint64_t l1 = l1n * WORDS_IN_L1HEADER;
		L1L0[l1 + LAST_L0_IS_SPARSE_W] |= 1ULL << LAST_L0_IS_SPARSE_SHIFT; //set bit
	}

	//Store the seven absolute L0 counts
	void store_L0_Absolute(uint64_t L0A, uint64_t l0, uint64_t l1n) {
		//Group A represent the l0-indices: 3,7,11,15,19,23,27
		static constexpr uint64_t L0A_size	  = 16; //actual value is only 15 bits, but spaced by 16 bits
		static constexpr uint64_t L0A_PER_WORD = 4; //four L0A fit in 64 bits
		uint64_t l1			= l1n * WORDS_IN_L1HEADER;
		uint64_t index		= 1 + l0/4;  //1..7 (start at the second 16-bit slot, leave the first one free for the lower s_w0)
		uint64_t shift		= L0A_size * (index % L0A_PER_WORD);
		uint64_t word_offs  = index / L0A_PER_WORD;
		L1L0[l1 + L0A_POS + word_offs] |= L0A << shift;
	}

	static void insert_12bit_value(uint64_t* vec, uint64_t value, uint64_t index) {
		uint64_t bit_pos = index * 12;		// Position in bits
		uint64_t word_index = bit_pos / 64; // Which 64-bit word to modify
		uint64_t bit_offset = bit_pos % 64; // Offset within that word
		vec[word_index] |= value << bit_offset;
		// If value overflows into the next 64-bit word, handle carry-over
		if (bit_offset > 52) {  // 64 - 12 = 52, overflow happens if bit_offset > 52
			uint64_t overflow_bits = bit_offset - 52;
			vec[word_index + 1] |= value >> (12 - overflow_bits);
		}
	}

	static void insert_value_at_index(uint64_t* vec, uint64_t value, uint64_t bit_width, uint64_t index) {
		uint64_t bit_pos    = index * bit_width;
		uint64_t word_index = bit_pos / 64;
		uint64_t bit_offset = bit_pos % 64;
		if (value >= (1ULL<<bit_width)) cout << "Warning! Inserting value " << value << "larger than it's bit_width" << bit_width << endl;
		vec[word_index] |= value << bit_offset;
		// Handle overflow if bits spill into the next 64-bit word
		if (bit_offset + bit_width > 64) {
			uint64_t overflow_bits = (bit_offset + bit_width) - 64;
			vec[word_index + 1] |= value >> (bit_width - overflow_bits);
		}
	}



	static inline uint64_t extract_z_value(const uint64_t *vec, uint64_t k) __attribute__((always_inline)) {
		constexpr uint64_t MASK = (1ULL << Z_VALUE_SIZE) - 1;
		const uint64_t bit_offset  = k * Z_VALUE_SIZE;
		const uint64_t byte_offset = bit_offset / 8;
		const uint64_t bit_shift   = bit_offset % 8;
		const auto *vec_as_byte_ptr = reinterpret_cast<const uint8_t*>(vec);
		uint64_t word;
		std::memcpy(&word, vec_as_byte_ptr + byte_offset, sizeof(word));//Safe unaligned load via memcopy.
		return (word >> bit_shift) & MASK;
	}

    //https://godbolt.org/z/faKTjMv4r
	//Read the next 'size' bits from vec, starting at index i
	static uint64_t extract_precise_value(const uint64_t i, const uint64_t size, const aligned_vector &vec) {
		const uint64_t byte_nr         = i / 8;
		const uint64_t offset_in_byte  = i % 8;
		// #if PREFETCH_T0
			// const uint64_t w = i/64;
			// _mm_prefetch((char*)(&vec[w]), _MM_HINT_T0); //We want to keep the tree in the cache
		// #endif
		const auto *vec_as_byte_ptr = reinterpret_cast<const uint8_t*>(vec.data());
		uint64_t extracted_word;
		//Safe unaligned load via memcopy.
		std::memcpy(&extracted_word, vec_as_byte_ptr + byte_nr, sizeof(extracted_word));
		// const uint64_t MASK = (1ULL << size) - 1;
		return _bzhi_u64(extracted_word >> offset_in_byte, size);

	}

	//Store the L0 count of a single block, which is 12 bits (representing 0...2048)
	void store_L0_delta(uint64_t delta_L0, uint64_t group_nr, uint64_t l0, uint64_t l1n) {
		//The delta L0 values are divided in three groups (b,c,d), each representing eight indices
		//b = 0, 4, 8,12,16,20,24,28 (group_nr 0)
		//c = 1, 5, 9,13,17,21,25,29 (group_nr 1)
		//d = 2, 6,10,14,18,22,26,30 (group_nr 2)
		//A = 3, 7,11,15,18,23,27,31 (shown here for completeness)
		//Each group is 96 consecutive bits long (8 values à 12 bits)
		uint64_t l1					 = l1n * WORDS_IN_L1HEADER;
		uint64_t index_in_group		 = l0/4; //0..7
		uint64_t prev_groups_entries = group_nr * 8;
		uint64_t index				 = index_in_group + prev_groups_entries;
		insert_12bit_value(&L1L0[l1+L0b_POS],delta_L0,index);
	}


	uint64_t summary_largest_L1 = 0;
	void build_L2_L1_summary_level(uint64_t n_superblocks) {
		//Layout: [L2  A |  B   C  |  E   F |  G  H | ...]
		//         47 17 |  32  32 | 32  32 | 32 32 |
		//One Summary block represents 16 L1-superblocks
		//One L2 value counts the number of 1-bits before the summary block
		//16 L1 values count the number of 1-bits since L2, up to and including its respective superblock
		//Technically only 20 bits would be needed per L1 entry (log2(16) * log2(2^16))
		//but for simplicity and faster SIMD processing, each gets 32 bits
		// Total: 512 bits

		uint64_t summary_L2_blocks= (n_superblocks / L1_IN_L2) +  NUM_SENTINEL_SUMMARY_BLOCKS ;
		L2L1 = aligned_vector(WORDS_IN_L2_SYMMLEVEL_HEADER*summary_L2_blocks,0);
		summary_largest_L1 = 0;

		// cout << "summary L2 blocks: " << summary_L2_blocks << endl;
		//Build on top of already built L1-Superblocks
		for (uint64_t l2n = 0; l2n < summary_L2_blocks; l2n++) {
			build_new_L2_block(l2n);
		}
	}

	void build_new_L2_block(uint64_t l2n) {
		//Get the L1 values for this L2-block.
		std::vector<uint64_t> L1s = std::vector<uint64_t>(L1_IN_L2);
		uint64_t l1n = l2n * L1_IN_L2;
		uint64_t l1 = l1n * WORDS_IN_L1HEADER;
		for (int i=0; i<L1_IN_L2; i++) {
			uint64_t L1;
			if (l1 + L1_POS < L1L0.size()) {
				L1 = L1L0[l1 + L1_POS] >> L1_SHIFT;
			} else {
				L1 = summary_largest_L1; //reached end of bitvector, fill with sentinels
			}
			summary_largest_L1 = L1;
			L1s[i] = L1;
			l1 += WORDS_IN_L1HEADER;
			if (i > 0 && L1 < L1s[i-1]) {
				cerr << "Error, next L1 ("<< L1<<")smaller than the previous" << endl;
			}
		}
		//The first L1 becomes the L2 value
		uint64_t L2 = L1s[0];
		#if DEBUGINFO_L2_SUMM_BUILD
		cout << "l2n " << l2n << ", L2: " << L2 << endl;
		#endif

		//Convert the remaining L1's to relative counts
		for (int i=0; i<L1_IN_L2; i++) {
			L1s[i] = L1s[i] - L2;
			// cout << "dL1 " << L1s[i] << endl;
			if (L1s[i]>16*L1_SIZEz) {
				cerr << "Error, delta L1 ("<< L1s[i]<<") too large" << endl;
			}
		}
		//Get the first L1
		uint64_t L1_A = L1s[1];
		if (L1_A > L1_SIZEz) {
			std::cerr << "ERROR, L1_A too large:" << L1_A << endl;
		}
		//L2 and first L1 share the first word
		uint64_t l2 = l2n * WORDS_IN_L2_SYMMLEVEL_HEADER;
		L2L1[l2 + L2_POS] |= (L2 << L2_SHIFT);
		L2L1[l2 + L2_POS] |= L1_A;

		//The next 14 L1 values represent the remaining 14 superblocks. We dont need to represent the 16-th L1 value, this can be reconstructed at runtime via L2' - L2
		for (int i=2; i<L1_IN_L2; i++) {
			int word  = i/2;
			int shift = i%2==0 ? 32 : 0;
			uint64_t L1 = L1s[i];
			L2L1[l2 + word] |= L1 << shift;
		}

		//put the second L1 split in the empty spaces before C and D
		// [L1_B_hi  C    L1_B_lo   D]
		//      12  20         12  20
		// uint64_t L1_B = L1s[2];
		// uint64_t L1_B_lo = _bzhi_u64(L1_B, L1B_SPLIT_BITS);
		// uint64_t L1_B_hi = L1_B >> L1B_SPLIT_BITS;
		// L2L1[l2 + L1B_POS] |= L1_B_lo << L1B_SHIFT;
		// L2L1[l2 + L1B_POS] |= L1_B_hi << L1B_SHIFT;
	}


	uint64_t build_L1_L0_array() {
		uint64_t words_to_build = N_words + NUM_SENTINEL_SUPERBLOCKS*W_IN_L1z; //some extra sentinel superblocks, for e.g. scan=2 stepping or alpha=4 explicit loads
		uint64_t n_superblocks = (words_to_build / W_IN_L1z);
		L1L0		= aligned_vector(WORDS_IN_L1HEADER*n_superblocks,0);
		dense		= aligned_vector(0);
		sparse		= aligned_vector(0);
		superblocks_counting_0s = 0;
		superblocks_counting_1s = 0;
		//We now walk through the full bitvector and sort it's information into the dense and sparse arrays
		uint64_t L1 = 0;	//Cumulative count of 1s
		uint64_t L0_cumul = 0;	//Cumulative count of minority-bits in the superblock
		uint64_t L0_block = 0;		//count of minority-bits in current L0 block
		// uint64_t delta_L0_prev = 0; //count of the minority bits in the previous L0 block
		uint64_t minflag = 1; //the minority bit type of the current superblock
		for (uint64_t w = 0; w < words_to_build; w++) {
			if(w % W_IN_L1z == 0) {
				//new L1 superblock starts
				uint64_t l1n = w / W_IN_L1z;
				minflag = get_SB_minority_type(l1n);
				minflag==1 ? superblocks_counting_1s++ : superblocks_counting_0s++;
				running_sparse_w0 += (running_zcount_of_block * Z_VALUE_SIZE + 63) / 64; //the next free slot in sparse
				//In a superblock header we store:
				//1. L1					Number of 1s before this SB
				//2. 31 x L0			The 31 L0 values. 7 of them cumulate (A), the other 24 only relative (b,c,d)
				//2. sparse_w0			The word where the sparse representation of this SB starts
				//3. dense_blocks_before The number of dense blocks that come before this SB
				//4. minflag			The minority bit type in this SB
				store_L1(L1, l1n);
				store_sparse_w0(l1n);
				store_dense_blocks_before(l1n);
				store_minflag(minflag, l1n);
				//initially set both all_x flags to true, and update them to false in case we find a counterexample. So we dont need to manually track their status
				set_alldense_flag(true, l1n);
				set_allsparse_flag(true, l1n);

				#if DEBUGINFO_SB_BUILDING
					cout << "SB "<< (w / W_IN_L1z) << endl;
					cout << "  L1 " << L1 << endl;
					cout << "  minflag=" << minflag << endl;
				#endif
				L0_cumul = 0;
				L0_block = 0;
				running_zcount_of_block = 0;
			}
			uint64_t pop = (w < N_words) ? std::popcount(bits[w]): 0;
			uint64_t mincount = minflag==1 ? pop : WORDSIZE-pop;
			L1	+= pop; //Counts 1-bits
			L0_cumul  += mincount; //Counts the minbit
			L0_block  += mincount;
			if((w+1) % W_IN_L0z == 0) {
				//a block is finished
				uint64_t l1n = w / W_IN_L1z;
				uint64_t l0 = (w % W_IN_L1z) / W_IN_L0z; // l0=0..31
				//We now store the L0 values

				#if DEBUGINFO_SB_BUILDING
					cout << "     l0: " << l0 << endl;
					cout << "    L0c: " << L0_cumul << endl;
					cout << "             L0b: " << L0_block << endl;
					cout << "                        L1curr: " << L1 << endl;
				#endif
				//Every fourth L0 value is stored cumulatve, the others only the relative count of the single block
				uint64_t group_nr = l0%4;
				if(group_nr==3 && l0!=31) {
					store_L0_Absolute(L0_cumul, l0, l1n);	 //Cumulative store (15 bit)
				} else if (l0!=31) {
					store_L0_delta(L0_block, group_nr, l0, l1n);//Relative store (12 bit)
				}
				store_block_sparse_or_dense(L0_block, minflag, l0, l1n);
				if (L0_block<=MAX_SPARSE) {
					set_alldense_flag(false, l1n); //saw a counterexample to all-dense, deactivate it for this superblock (idempotent, we might deactivate it multiple times)
				} else {
					set_allsparse_flag(false, l1n); //saw a counterexample to all-sparse, deactivate it for this superblock
				}
				if (l0==31 && L0_block <= MAX_SPARSE) {
					store_last_L0_is_sparse(l1n); //without this flag, we would need to read L1 from the next header cacheline to determine whether the last L0 block is sparse
				}
				L0_block = 0;
			}
		}
		uint64_t relevant_superblocks = n_superblocks - (NUM_SENTINEL_SUPERBLOCKS-1);
		return relevant_superblocks;
	}


	void store_block_sparse_or_dense(uint64_t L0_block, uint64_t minflag, uint64_t l0, uint64_t l1n) {
		if (L0_block <= MAX_SPARSE && BV_COMPRESSION) {
			#if DEBUGINFO_SB_BUILDING
			cout << "                           l0 = " << l0 << " sparse " << endl;
			#endif
			store_sparse_block(l0, minflag, l1n);
		} else {
			#if DEBUGINFO_SB_BUILDING
			cout << "                           l0 = " << l0 << " dense" << endl;
			#endif
			store_dense_block(l0, l1n);
		}

		#if DEBUGINFO_SB_BUILDING
		cout << "                         dense block count = " << running_dense_block_count << endl;
		#endif
	}


	//store the 2048 bits of a L0-block in an encoded compressed form
	void store_sparse_block(uint64_t l0, uint64_t minflag, uint64_t l1n) {
		uint64_t const w_block_start = l1n * W_IN_L1z + l0 * W_IN_L0z;
		uint64_t const w_block_end   = w_block_start + W_IN_L0z;

		//Go through block, collect the sparse bits.
		//Count how many dense bits we have already seen in z_prev
		uint64_t z_prev = 0;
		for(uint64_t w = w_block_start; w < w_block_end; w++) {
			uint64_t WORD = bits[w];
			//Inverting helps to not have to bother in the following lines over the mintype
			if(minflag == 0) {
				WORD = ~WORD;
			}
			uint64_t sparse_bits_in_word = std::popcount(WORD);
			//iterate over all sparse bits per word
			for(uint64_t k=0; k<sparse_bits_in_word; k++) {
				uint64_t i = Tree3::get_kth_1_position(k+1,WORD);
				uint64_t z = z_prev + i - k;
				store_z(z);
			}
			z_prev += 64 - sparse_bits_in_word;
		}
	}

	//Append a z-value to the sparse array.
	void store_z(uint64_t z) {
		//keep the array large enough
		uint64_t rough_new_index = running_sparse_w0 + (running_zcount_of_block*Z_VALUE_SIZE)/64;
		while(sparse.size() < rough_new_index + 5) {
			sparse.push_back(0);
		}
		//store
		//insert_12bit_value(&sparse[running_sparse_w0],z,running_zcount_of_block);
		insert_value_at_index(&sparse[running_sparse_w0], z , Z_VALUE_SIZE, running_zcount_of_block);
		running_zcount_of_block += 1;
	}

	//just copy and store the 2048 bits of this L0-block exactly as they are
	void store_dense_block(uint64_t l0, uint64_t l1n) {
		uint64_t w_read_start = l1n * W_IN_L1z + l0 * W_IN_L0z;
		uint64_t w_read_end   = w_read_start + W_IN_L0z;
		for(uint64_t w = w_read_start; w < w_read_end; w++) {
			dense.push_back(bits[w]);
			#if DEBUGINFO_RANK
			if (l1n==0) {
				// cout << "dense[" << dense.size()-1 << "]="<< endl;
				cout << "l0="<< l0 << endl;
				cout << "w=" << w << endl;
				printWord(dense.back()," dense index="+to_string(dense.size()-1)) ;
			}
			#endif
		}
		running_dense_block_count++;
	}


	void static print256_epi16(__m256i vec, const string &name, int mult=1, int offs=0) {
		alignas(32) uint16_t tmp[16];
		_mm256_store_si256((__m256i*)tmp, vec);
		cout << name << endl;
		for (int i = 0; i < 16; i++) {
			cout << "    " << (i*mult + offs) << " : " << tmp[i] << endl;
		}
	}


	void static print256_epi32(__m256i vec, const string &name, int mult=1, int offs=0) {
		alignas(32) uint32_t tmp[8];
		_mm256_store_si256((__m256i*)tmp, vec);
		cout << name << endl;
		for (int i = 0; i < 8; i++) {
			cout << "    " << (i*mult + offs) << " : " << tmp[i] << endl;
		}
	}

	void static print128_epi16(__m128i vec, const string &name, bool sign=false) {
		alignas(32) uint16_t tmp[8];
		_mm_store_si128((__m128i*)tmp, vec);
		cout << name << endl;
		for (int i = 0; i < 8; i++) {
			if (sign) {
				cout << "    " <<  i << " : " <<((int16_t)tmp[i]) << endl;
			}
			else {
				cout << "    " <<  i << " : " << tmp[i] << endl;
			}

		}
	}

	static inline void printWord(uint64_t W, string name="") {
		if (name!="") cout << name << endl;
		for(int i=63; i>=0; i--) {
			std::cout << ((W>>i) & 1ULL);
			if (i==32) cout << ".";
		}
		cout << endl;
	}

	[[nodiscard]] uint64_t get_vanilla_sample_tree_result(uint64_t k) const {
		k--; //go zero-indexed
		uint64_t t = k>> T3.a_exp;
		auto TG      = T3.tmp_TOP[t];
		auto TG_next = T3.tmp_TOP[t+1];
		uint64_t l1n = TG.l1; //we store l1n in T but for brevity call it l1 there
		uint64_t l1n_next = TG_next.l1;
		uint64_t r = l1n_next - l1n;

		#if DEBUGINFO_T3_WALK
			printf("vanilla_walk_T3:\n");
			printf("k        %lu\n", k+1); //one-indexed
			printf("l1n      %lu\n", l1n);
			printf("l1n_next %lu\n", l1n_next);
			printf("r        %lu\n", r);
		#endif

		if (r < ALPHA) {
			return l1n;
		}


		k = k & T3.mod_a;
		uint64_t midgroup_nr = k >> T3.b_exp;
		uint64_t gamma = Tree3::rlog2(r);
		uint64_t s = gamma + TG.kappa;
		auto MG = T3.tmp_MID[s][TG.j];
		uint64_t dl      = MG.DL[midgroup_nr];
		uint64_t dl_next = MG.DL[midgroup_nr+1]; //need sentinel midgroup!
		uint64_t r_prime = dl_next - dl;

		#if DEBUGINFO_T3_WALK
			printf("--\n");
			printf("k     %lu\n", k);
			printf("gamma %lu\n", gamma);
			printf("kappa %lu\n", TG.kappa);
			printf("s     %lu\n", s);
			printf("midblock_nr %lu\n", TG.j);
			printf("midgroup_nr %lu\n", midgroup_nr);
			printf("dl      %lu\n", dl);
			printf("dl_next %lu\n", dl_next);
			printf("r_prime %lu\n", r_prime);
		#endif

		if (r_prime < ALPHA) {
			return l1n + dl;
		}

		uint64_t rho = Tree3::rlog2(r_prime);
		uint64_t c  = MG.c_local[midgroup_nr];

		uint64_t Cs_in_tuple = MG.rhomax - T3.rho_min + 1;
		uint64_t i_tuple = MG.h * Cs_in_tuple;
		uint64_t d_K = rho - T3.rho_min;
		uint64_t C_rho = T3.tmp_KK[MG.rhomax][i_tuple + d_K];
		// uint64_t C_rho = MG.C_global[rho - T3.rho_min];
		uint64_t c_hat = c + C_rho;
		uint64_t i_Bsection = T3.b * c_hat;
		uint64_t d_iB = k & T3.mod_b;
		uint64_t ddl = T3.tmp_BOT[rho][i_Bsection + d_iB];

		#if DEBUGINFO_T3_WALK
			printf("--\n");
			printf("rho %lu \n", rho);
			printf("c   %lu\n", c);
			// printf("h   %lu\n", MG.h);
			// printf("C_nr %lu \n", i_Ksection + d_iK);
			printf("C    %lu \n", C_rho);
			printf("c_hat %lu\n", c_hat);
			printf("ddl_nr %lu\n", i_Bsection + d_iB);
			printf("i_Bsection %lu \n", i_Bsection);
			printf("ddl %lu \n", ddl);
		#endif
		return l1n + dl + ddl;
	}





	[[nodiscard]] uint64_t get_packed_sample_tree_result(uint64_t k) const __attribute__((always_inline)) {
		//Calculate the position of the topgroup [l1,j,kappa] in T
		const uint64_t a_exp		= T3.a_exp;
		const uint64_t topgroup_bits= T3.glob_topgroup_bits;
		const uint64_t j_bits		= T3.glob_j_bits;
		const uint64_t kappa_bits	= T3.glob_kappa_bits;
		const auto &TOP = T3.pack_TOP;
		k--; //work zero-indexed
		const uint64_t topgroup_nr	= k >> a_exp;
		const uint64_t i_top		= topgroup_nr * topgroup_bits;

		//Read the top-group T
		const uint64_t topgroup		= extract_precise_value(i_top                , topgroup_bits, TOP);
		const uint64_t topgroup_next= extract_precise_value(i_top + topgroup_bits, topgroup_bits, TOP);
		const uint64_t l1n			= topgroup     >> (j_bits + kappa_bits);
		const uint64_t l1n_next		= topgroup_next>> (j_bits + kappa_bits);
		const uint64_t r			= l1n_next - l1n;

		#if DEBUGINFO_T3_WALK
			printf("\n\nwalk_T3_packed:\n");
			printf("tg_nr    %lu\n", topgroup_nr);
			printf("i_top    %lu\n", i_top);
			printf("topgroup %lu\n", topgroup);
			printf("l1n      %lu\n", l1n);
			printf("l1n_next %lu\n", l1n_next);
			printf("r        %lu\n", r);
		#endif

		if (r < ALPHA) {
			return l1n;
		}

		//Calculate the position of the midgroup [dl, c] in M
		const uint64_t mod_a		= T3.mod_a;
		const uint64_t MIDENTRIES_wS= T3.MIDENTRIES_w_SENTINEL;
		const uint64_t rhomax_bits	= T3.glob_rhomax_bits;
		const uint64_t h_bits		= T3.glob_h_bits;
		const uint64_t b_exp		= T3.b_exp;
		k = k & mod_a;
		const uint64_t gamma		= Tree3::quick_rlog2(r);
		const uint64_t kappa		= _bzhi_u64(topgroup, kappa_bits);
		const uint64_t s			= gamma + kappa;
		const auto &MIDs = T3.pack_MID[s];
		const uint64_t midgroup_size= MIDENTRIES_wS * s + rhomax_bits + h_bits;
		const uint64_t midgroup_nr	= _bzhi_u64(topgroup >> kappa_bits, j_bits);
		const uint64_t midentry_nr	= k >> b_exp;
		const uint64_t i_midgroup	= midgroup_nr * midgroup_size;
		const uint64_t i_midentry	= i_midgroup + midentry_nr * s;
		//Read the [dl, c], [dl', c'] midgroups from M_s
		const uint64_t midentry		= extract_precise_value(i_midentry    , s, MIDs);
		const uint64_t midentry_next= extract_precise_value(i_midentry + s, s, MIDs);
		const uint64_t dl			= midentry     >> kappa;
		const uint64_t dl_next		= midentry_next>> kappa;
		const uint64_t r_prime		= dl_next - dl;

		#if DEBUGINFO_T3_WALK
			printf("--\n");
			printf("k       %lu\n", k);
			printf("gamma   %lu\n", gamma);
			printf("kappa   %lu\n", kappa);
			printf("s       %lu\n", s);
			printf("midblock_nr %lu\n", midblock_nr);
			printf("midgroup_nr %lu\n", midgroup_nr);
			printf("dl      %lu\n", dl);
			printf("dl_next %lu\n", dl_next);
			printf("r_prime %lu\n", r_prime);
		#endif

		if (r_prime < ALPHA) {
			return l1n + dl;
		}

		//read the pointer-info [rhomax, h] from the midblock
		const uint64_t i_rhomax_h= i_midgroup + midgroup_size - rhomax_bits - h_bits;
		const uint64_t rhomax_h	= extract_precise_value(i_rhomax_h, rhomax_bits + h_bits, MIDs);
		const uint64_t rhomax	= rhomax_h >> h_bits;
		const uint64_t h		= _bzhi_u64(rhomax_h, h_bits);

		//calculate the position of the C-counter in K_rhomax
		const uint64_t rho_min	= T3.rho_min;
		const uint64_t C_bits	= T3.glob_C_bits;
		const uint64_t Cs_in_tuple = rhomax - rho_min + 1;
		const uint64_t rho		= Tree3::quick_rlog2(r_prime);
		const uint64_t C_nr		= h * Cs_in_tuple + (rho - rho_min);
		const uint64_t i_C		= C_nr * C_bits;

		//read the C-counter in K_rhomax
		const uint64_t C		= extract_precise_value(i_C, C_bits, T3.pack_KK[rhomax]);

		//calculate the position of the bot-sample ddl in B
		const uint64_t b		= T3.b;
		const uint64_t mod_b	= T3.mod_b;
		const uint64_t c		= _bzhi_u64(midentry, kappa);
		const uint64_t ddl_nr	= (C + c ) * b + (k & mod_b);
		const uint64_t i_ddl	= ddl_nr * rho;

		//Read ddl from the bot-level B_rho
		const uint64_t ddl		= extract_precise_value(i_ddl, rho, T3.pack_BOT[rho]);

		#if DEBUGINFO_T3_WALK
			printf("C      %lu \n", C);
			printf("c      %lu \n", c);
			printf("C+c    %lu \n", C+c);
			printf("ddl_nr %lu \n", ddl_nr);
			printf("i_ddl  %lu \n", i_ddl);
			printf("ddl    %lu \n", ddl);
		#endif

		return l1n + dl + ddl;
	}





	[[nodiscard]] uint64_t select_1(int64_t k0) const override {
		uint64_t k = k0;


		#if TREE_TYPE==3
		#if DOUBLECHECK_T3_VANILLA
		uint64_t vanilla_tree_result = get_vanilla_sample_tree_result(k);
		#endif
		uint64_t tree_result =  get_packed_sample_tree_result(k);

		#if DOUBLECHECK_T3_VANILLA
		if (vanilla_tree_result != tree_result) {
			cout << "Warning, T3 vanilla and T3 packed differ! vanilla " << vanilla_tree_result << "  packed " << tree_result << " k0 =" << k0 <<  endl;
		}
		#endif




		#elif TREE_TYPE==5
		uint64_t const km1 = k-1;
		uint64_t const k_mod_lvl0 = km1 & MODULO_SAMPLEDIST_LVL0;
		uint64_t const k_mod_lvl1 = km1 & MODULO_SAMPLEDIST_LVL1;
		uint64_t const k_mod_lvl2 = km1 & MODULO_SAMPLEDIST_LVL2;
		uint64_t const k_mod_lvl3 = km1 & MODULO_SAMPLEDIST_LVL3;
		uint64_t const k_mod_lvl4 = km1 & MODULO_SAMPLEDIST_LVL4;

		uint64_t const k_offset_lvl1 = k_mod_lvl0 >> SAMPLEDIST_EXP_LVL1;
		uint64_t const k_offset_lvl2 = k_mod_lvl1 >> SAMPLEDIST_EXP_LVL2;
		uint64_t const k_offset_lvl3 = k_mod_lvl2 >> SAMPLEDIST_EXP_LVL3;
		uint64_t const k_offset_lvl4 = k_mod_lvl3 >> SAMPLEDIST_EXP_LVL4;


		uint64_t lvl = 0;
		uint64_t i_expl = -1;
		uint64_t const b0	 = km1 >> SAMPLEDIST_EXP_LVL0;
		#if PREFETCH_T0
		_mm_prefetch((char*)(&B0[b0]), _MM_HINT_T0); //We want to keep the highest tree level in the cache
		#endif
		uint32_t const BB0   = B0[b0];
		uint32_t const flag0 = BB0 >> TREE_FLAGSHIFT;
		uint32_t const i0	 = BB0 & TREE_MASKVALUE;


		// cout << "flag0 : " << flag0 << endl;
		// cout << "i0    : " << i0 << endl;

		if     (flag0==EXPLICIT_FLAG) {i_expl =(i0 << SAMPLEDIST_EXP_LVL4) + k_mod_lvl0;}
		else if(flag0==L1_FLAG)	{l1 = (i0*WORDS_IN_HEADER)+ WORDS_IN_HEADER + L1_POS;
			#if L1MINISKIP
			l1 += WORDS_IN_HEADER*(k_mod_lvl0/L1_SIZEz);
			#endif
		}
		else {
			lvl = 1;
			uint32_t const b1	 = i0 + k_offset_lvl1;
			uint32_t const BB1   = B1[b1];
			uint32_t const flag1 = BB1 >> TREE_FLAGSHIFT;
			uint32_t const i1    = BB1 & TREE_MASKVALUE;
			if     (flag1==EXPLICIT_FLAG) {i_expl = (i1 << SAMPLEDIST_EXP_LVL4) + k_mod_lvl1;}
			else if(flag1==L1_FLAG) {l1 = (i1*WORDS_IN_HEADER)+ WORDS_IN_HEADER + L1_POS;}
			else {
				lvl = 2;
				uint32_t const b2	 = i1 + k_offset_lvl2;
				uint32_t const BB2   = B2[b2];
				uint32_t const flag2 = BB2 >> TREE_FLAGSHIFT;
				uint32_t const i2    = BB2 & TREE_MASKVALUE;
				if     (flag2==EXPLICIT_FLAG) {i_expl = (i2 << SAMPLEDIST_EXP_LVL4) + k_mod_lvl2;}
				else if(flag2==L1_FLAG) {l1 = (i2*WORDS_IN_HEADER)+ WORDS_IN_HEADER + L1_POS; }
				else {
					lvl = 3;
					uint32_t const b3	 = i2 + k_offset_lvl3;
					uint32_t const BB3   = B3[b3];
					uint32_t const flag3 = BB3 >> TREE_FLAGSHIFT;
					uint32_t const i3    = BB3 & TREE_MASKVALUE;
					if     (flag3==EXPLICIT_FLAG) {i_expl = (i3 << SAMPLEDIST_EXP_LVL4) + k_mod_lvl3;}
					else if(flag3==L1_FLAG)		  {l1 = (i3*WORDS_IN_HEADER)+ WORDS_IN_HEADER + L1_POS;}
					else {
						lvl = 4;
						uint32_t const b4	 = i3 + k_offset_lvl4;
						uint32_t const BB4   = B4[b4];
						uint32_t const flag4 = BB4 >> TREE_FLAGSHIFT;
						uint32_t const i4    = BB4 & TREE_MASKVALUE;
						if (flag4==EXPLICIT_FLAG) {i_expl = (i4 << SAMPLEDIST_EXP_LVL4) + k_mod_lvl4;}
						else					  {l1 = (i4*WORDS_IN_HEADER)+ WORDS_IN_HEADER + L1_POS;}
					}
				}
			}
		}

		#if DIAGNOSE_SB_SELECT
			printf("tree lvl %lu\n", lvl);
			printf("i_expl = %ld\n", i_expl);
		#endif
		if(i_expl != -1) {
			if(show_tree_delta_l1n){cout << "ex" << lvl << spacer;}
			return EXPL[i_expl];
		}
		uint64_t l1n_start = l1/WORDS_IN_HEADER;
		//at this point we either already returned the answer (if it was stored explicitly in the tree)
		//or l1 is now sufficiently close in front of the correct superblock such that the remaining scanning is bounded in O(1)
		#else
		cout << "m3 only supports TREE_TYPE = {3,5}" << endl;
		#endif

		uint64_t l1;

		if (SUMMARY_LEVELS==2) {

			uint64_t l2n;
			uint64_t l2;

			if (ALPHA>1) {
				l2n = tree_result;
				l2 = l2n * WORDS_IN_L2_SYMMLEVEL_HEADER + WORDS_IN_L2_SYMMLEVEL_HEADER;
				while(k > L2L1[l2]>>L2_SHIFT) {
					l2 += WORDS_IN_L2_SYMMLEVEL_HEADER;
				}
				l2 -= WORDS_IN_L2_SYMMLEVEL_HEADER;
				l2n = l2 / WORDS_IN_L2_SYMMLEVEL_HEADER;
			} else {
				l2n = tree_result;
				l2 = l2n * WORDS_IN_L2_SYMMLEVEL_HEADER;
			}

				  __m256i dL1s_A = _mm256_loadu_si256((__m256i*)&L2L1[l2]);
			const __m256i dL1s_B = _mm256_loadu_si256((__m256i*)&L2L1[l2+4]);


			#if DEBUGINFO_L2_SUMMARY
			// print256_epi32(dL1s_A, "dL1s_A raw");
			// print256_epi32(dL1s_B, "dL1s_B raw");
			#endif
			// #if PREFETCH_T0
			// _mm_prefetch((char*)(&L2L1[l2]), _MM_HINT_T0);
			// #endif
						  dL1s_A = _mm256_and_si256(dL1s_A, summary_mask_hide_L2);
			uint64_t L2		= L2L1[l2 + L2_POS] >> L2_SHIFT;
			// uint64_t L2next = L2L1[l2 + L2_POS + WORDS_IN_L2_SYMMLEVEL_HEADER] >> L2_SHIFT;
			// uint64_t L1_16  = L2next - L2;
			const uint64_t k_cmp = k - L2;
			#if DEBUGINFO_L2_SUMMARY
			if (k_cmp>L2_SIZE) {
				cout << "k0 " << k0 << endl;
				cerr << "ERROR, k too big. " << k << endl;
			}
			#endif
			const __m256i		  K2   = _mm256_set1_epi32(k_cmp);
			const __m256i K2_grt_dL1sA = _mm256_cmpgt_epi32(K2, dL1s_A);
			const __m256i K2_grt_dL1sB = _mm256_cmpgt_epi32(K2, dL1s_B);
			const uint32_t A_mvmsk  = _mm256_movemask_epi8(K2_grt_dL1sA);
			const uint32_t B_mvmsk  = _mm256_movemask_epi8(K2_grt_dL1sB);
			const int A_pop = std::popcount(A_mvmsk)/4;
			const int B_pop  = std::popcount(B_mvmsk)/4;
			//L1s_A has one entry with only zeros, where L2 was masked out, where the compare then always gives a positive result k>0. To compensate that, we count -1
			//Also, we dont need to check against the end of the 16-th superblock, because that check would be equal to k>L2next,
			//which we already know is not the case, bc. then we would be in the wrong L2 Block altogether
			const uint64_t superblocks_skipped = A_pop + B_pop - 1;

			const uint64_t l1n_direct = l2n * L1_IN_L2 + superblocks_skipped;
			uint64_t l1_direct = l1n_direct * WORDS_IN_L1HEADER;
			l1 = l1_direct + WORDS_IN_L1HEADER + L1_POS; //let it point directly on L1 of the next superblock

			#if DEBUGINFO_L2_SUMMARY
			print256_epi32(dL1s_A, "dL1s_A");
			print256_epi32(dL1s_B, "dL1s_B");

			cout<<"l2n "	<< l2n << endl;
			cout<<"L2 "	<< L2 << endl;
			cout<<"k' "	<< k << endl;
			cout<<"skipped "	<< superblocks_skipped << endl;
			cout<<"l1n_direct "	<< l1n_direct << endl;
			cout<<"l1_direct "	<< l1_direct << endl;
			cout<<"l1 "	<< l1 << endl;
			#endif
		}
		else { //no summary level
			l1 = tree_result * WORDS_IN_L1HEADER  + WORDS_IN_L1HEADER + L1_POS; //Index of the L1 value (+L1_POS) of the next superblock (+WORDS_IN_HEADER)

			//Todo: replace with explicit ALPHA loop
			while(k > L1L0[l1]>>L1_SHIFT) {
				l1 += WORDS_IN_L1HEADER;
			}
		}

		//We now found the correct superblock
		uint64_t const l1_next = l1;
		//Currently l1 points on the L1 value of the next L1-block, bring it back to the start of the current L1-block
		l1 -= (WORDS_IN_L1HEADER + L1_POS);
		uint64_t const l1n = l1/WORDS_IN_L1HEADER;
		//Read that L1-block
		uint64_t const L1  = L1L0[l1+L1_POS] >> L1_SHIFT;
		k -= L1;
		bool const minflag = (L1L0[l1+MINFLAG_POS]>>MINFLAG_SHIFT)&1ULL;


		uint64_t l1n_scansteps = l1n - tree_result;
		#if DEBUGINFO_T3_WALK
			printf("\n");
			printf("l1n_start %lu\n", l1n_start);
			printf("l1n final %lu\n", l1n);
		#endif
		#if DEBUGINFO_T3_SCANLENGTH
			printf("          >>> STEPS %lu <<< \n", l1n_scansteps);
		#endif
		#if ALERT_T3_ALPHA_BREAKING
			if (l1n_scansteps >= ALPHA) {
				printf(" Error: Scanned through too many superblocks, broke the T3-alpha bound\n");
				printf("      k = %lu \n", k);
				printf("      alpha    = %lu \n", ALPHA);
				printf("      l1n-DIFF = %lu \n", l1n_scansteps);
			}
		#endif
		#if DEBUGINFO_SB_SELECT
			cout << "ß post scan l1  " << l1 << endl;
			cout << "ß post scan l1n " << l1n << endl;
			cout << "      L1=" << L1 << endl;
			cout << "in SB: k=" << k << endl;
			cout << "minflag =" << minflag << endl;
		#endif
		/**
		 * The Superblock Structure
		 *		[L1    sparse_w0    dense_blocks   d0...d7    c0...c7  b0...b7   A1...A7  minflag  A0    sparse_w0]
		 *       43    21           32             8*12       8*12     8*12      6*16     1        15    16
		 */

		//Unpack the A,b,c,d values
		__m128i A_raw = _mm_loadu_si128((__m128i*)&L1L0[l1 + L0A_POS]);
		__m128i b_raw = _mm_loadu_si128((__m128i*)&L1L0[l1 + L0b_POS]);
		__m128i c_raw = _mm_loadu_si128((__m128i*)((uint8_t*)&L1L0[l1 + L0c_POS] + 4)); //c bits start in the middle of a 64-bit word -> shift start by 4 byte
		__m128i d_raw = _mm_loadu_si128((__m128i*)&L1L0[l1 + L0d_POS]);
		//could also use loadu, compiler then knows that doesnt need to care about alignment, doesnt need to prepare error message
		//Alternatively, _mm_stream_load_si128 could also be interesting, reading with a non-temporal prefetch hint. However, there is no non-temporal unaligned load, would need some work-around

		__m128i	    A = _mm_and_si128(A_raw, remove_clutter_around_A);
		//could load c_raw directly 4 bits later, to save this shift

		__m128i	b_low = _mm_shuffle_epi8(b_raw,unpack12to16_shuff_low); //unpack packed 12 bit into high and low bytes
		__m128i	c_low = _mm_shuffle_epi8(c_raw,unpack12to16_shuff_low);
		__m128i	d_low = _mm_shuffle_epi8(d_raw,unpack12to16_shuff_low);

		__m128i	b_hi  = _mm_shuffle_epi8(b_raw,unpack12to16_shuff_hi);
		__m128i	c_hi  = _mm_shuffle_epi8(c_raw,unpack12to16_shuff_hi);
		__m128i	d_hi  = _mm_shuffle_epi8(d_raw,unpack12to16_shuff_hi);

				b_low = _mm_and_si128(b_low, keepLower12_epi16);
				c_low = _mm_and_si128(c_low, keepLower12_epi16);
				d_low = _mm_and_si128(d_low, keepLower12_epi16);

				b_hi  = _mm_srli_epi16(b_hi,4);
				c_hi  = _mm_srli_epi16(c_hi,4);
				d_hi  = _mm_srli_epi16(d_hi,4);

		__m128i    b  = _mm_or_si128(b_hi, b_low);  //combine high and low bytes of each 12bit value into 16bit slots
		__m128i    c  = _mm_or_si128(c_hi, c_low);
		__m128i    d  = _mm_or_si128(d_hi, d_low);

		__m128i	   B  = _mm_add_epi16(A,b); //convert relative counts b into absolute L0-values B
		__m128i	   C  = _mm_add_epi16(B,c); //convert relative counts c into absolute L0-values C
		__m128i	   D  = _mm_add_epi16(C,d); //convert relative counts d into absolute L0-values D

		//if the L0's counted 0-bits (minflag==0), need to convert into the complement count of 1-bits
		if(minflag==0) [[unlikely]] {
			A = _mm_sub_epi16(MAX_A, A);
			B = _mm_sub_epi16(MAX_B, B);
			C = _mm_sub_epi16(MAX_C, C);
			D = _mm_sub_epi16(MAX_D, D);
		}

		//The last A value is not stored as an L0-value, but can be inferred by the neighbouring L1 values
		//Since the L1 values always track 1-bits, A31 also always tracks 1-bits, and must be inserted only after the complement count
		 /*
		  * Technically we dont need A31 in dense case!! Its always k<=A31, thus the compare is always positive.
		  * Could just do 31 - pop - pop
		  * We would save these instructions here, AND we would not need to load L1next!!
		  * Only downside: A31 would not be available as potential L2prev -- but that only affects 1/32 of cases, and we just dont reverse in that rare case
		  * Instead, in sparse case, we need L2next to determine whether we are sparse
		  * To save neighbouring cacheline-access, maybe store A31 explicitly in header, in case of SUMMARY_LEVEL 2
		  */
		// uint64_t const A31 = L1next - L1;
		//We had until now A shifted by 16 bit to the left, such that we could add/subtract it directly from b,c,d. Now we need it in line with the rest
		//[A7 A6 A5 A4 A3 A2 A1 _ ] --> [ _ A7 A6 A5 A4 A3 A2 A1]
		A = _mm_srli_si128(A,2);
		// A = _mm_insert_epi16(A,A31,7);

		//We now prepare the L0-values A,B,C,D for two incoming uses.
		//First, we merge them into pairwise 256-bit vectors, to be able to treat them with 2x mm256 instructions instead of 4x mm128
		//Second, we store these values in a scalar uint16_t array, such that we have direct random access to any single one
		__m128i L0_evens_lo = _mm_unpacklo_epi16(B, D);
		__m128i L0_evens_hi = _mm_unpackhi_epi16(B, D);
		__m256i L0_evens    = _mm256_set_m128i(L0_evens_hi, L0_evens_lo); //directly store 128vectors, adapt read index to non-interleaved store
		__m128i L0_odds_lo = _mm_unpacklo_epi16(C, A);
		__m128i L0_odds_hi = _mm_unpackhi_epi16(C, A);
		__m256i L0_odds    = _mm256_set_m128i(L0_odds_hi, L0_odds_lo);
		alignas(32) uint16_t tmp_L0_evens[16];
		alignas(32) uint16_t tmp_L0_odds[16];
		_mm256_store_si256((__m256i*)tmp_L0_evens, L0_evens);
		_mm256_store_si256((__m256i*)tmp_L0_odds,  L0_odds);


		const __m256i		  K_rs = _mm256_set1_epi16((k-1)-32768); //shift down into range [32767,-32768] range where cmpgt works
		const __m256i L0evns_rs    = _mm256_sub_epi16(L0_evens, range_shift_32768_epi16);
		const __m256i L0odds_rs    = _mm256_sub_epi16(L0_odds,  range_shift_32768_epi16);
		const __m256i L0evns_grt_K = _mm256_cmpgt_epi16(L0evns_rs, K_rs);
		const __m256i L0odds_grt_K = _mm256_cmpgt_epi16(L0odds_rs, K_rs);
		uint32_t even_mvmsk = _mm256_movemask_epi8(L0evns_grt_K);
		uint32_t odd_mvmsk  = _mm256_movemask_epi8(L0odds_grt_K);
		const int even_pop = std::popcount(even_mvmsk)/2;
		const int odd_pop  = std::popcount(odd_mvmsk)/2;
		//found the block index l0
		//the last 32-th block is always filled with zeros, so it never triggers L0>k, so we only have 31 participating L0 values in the compare
		uint64_t const l0 = 31 - even_pop - odd_pop;

		#if DEBUGINFO_SB_SELECT
			print256_epi16(L0_evens, "L0 evens");
			print256_epi16(L0_odds, "L0 odds ");
			printWord(even_mvmsk,"even mvmsk");
			printWord(odd_mvmsk,"odd mvmsk");
		#endif


		uint64_t const i_upto_block = l1n * L1_SIZEz + l0 * L0_SIZEz;

		#if DEBUGINFO_SB_SELECT
			cout << "ß l0 " << l0 << endl;
		#endif

		//Now that we now l0, extract the bounding L0 values L0prev and L0next of our block
		//every block is bounded by one even L0 and one odd L0, only their order might be swapped.
		uint16_t const nearby_L0_even = tmp_L0_evens[l0/2];
		uint16_t const nearby_L0_odd  = tmp_L0_odds[((l0-1)/2)%16]; //wrap l0=0 to a dummy entry

		uint32_t L0prev = (l0%2==0) ? nearby_L0_odd : nearby_L0_even;
		uint32_t L0next = (l0%2==0) ? nearby_L0_even : nearby_L0_odd;
				 L0prev = (l0>0) ? L0prev : 0;


		//Edge case: If we are in the last block, then L0next is not stored in our current superblock, but we must reconstruct it via the next L1' value
		if (l0==31) {
			uint64_t const L1next = L1L0[l1_next] >> L1_SHIFT;
			L0next = L1next - L1;
			// L0next = A31;
			// A = _mm_insert_epi16(A,A31,7);
		}


		#if DEBUGINFO_SB_SELECT
			cout << "nearby_L0_even " << nearby_L0_even << endl;
			cout << "nearby_L0_odd  " << nearby_L0_odd << endl;
			cout << "L0prev " << L0prev << endl;
			cout << "L0next " << L0next << endl;
		#endif


		bool block_is_dense;
		bool alldense;
		bool allsparse;
		uint64_t mbits_in_block;
		__m256i blocks_even;
		__m256i blocks_odd;

		if(BV_COMPRESSION==1) {
			uint32_t const onebits_in_block = L0next - L0prev;
			mbits_in_block   = minflag==1 ? onebits_in_block : L0_SIZEz - onebits_in_block;
			block_is_dense   = mbits_in_block > MAX_SPARSE;

			//block values a,b,c,d count how many m-bits there are in an individual L0-block
			__m128i blocks_even_lo = _mm_unpacklo_epi16(b, d);
			__m128i blocks_even_hi = _mm_unpackhi_epi16(b, d);
					blocks_even    = _mm256_set_m128i(blocks_even_hi, blocks_even_lo);
			__m128i             a  = _mm_sub_epi16(A,D);  //number of 1-bits in each A-block individually (not cumulative)
			if (minflag==0) [[unlikely]] {
				a = _mm_sub_epi16(mm_L0size_epi16, a); //convert to 0-bit count
			}
			__m128i blocks_odd_lo = _mm_unpacklo_epi16(c, a);
			__m128i blocks_odd_hi = _mm_unpackhi_epi16(c, a);
					blocks_odd    = _mm256_set_m128i(blocks_odd_hi, blocks_odd_lo);


			#if DEBUGINFO_SB_SELECT
			print256_epi16(blocks_even, "m-bit count blocks even", 2, 0);
			print256_epi16(blocks_odd, "m-bit count blocks odd", 2, 1);
			cout << "i_upto_block " << i_upto_block << endl;
			#endif


			//Can try to take a shortcut if every L0 block is sparse or dense. Saves instructions but costs a new branch
			#if CHECK_ALL_DENSE_SPARSE
			alldense  = L1L0[l1+ALLDENSE_W]>>ALLDENSE_SHIFT & 1ULL;
			allsparse = L1L0[l1+ALLSPARSE_W]>>ALLSPARSE_SHIFT & 1ULL;
			#else
			alldense = false;
			allsparse = false;
			#endif
		} else {
			//The vanilla option with Bitvector-compression fully disabled.
			//Thus the full bitvector is stored dense
			block_is_dense = true;
			alldense = true;
			allsparse = false;
			// blocks_even = _mm256_setzero_si256();
			// blocks_odd = _mm256_setzero_si256();
			mbits_in_block = 0;
		}

		if (block_is_dense) {
			//Our block is dense, need to know how many previous blocks in this SB are also dense
			uint64_t dense_blocks_in_SB = 0;
			#if COUNT_HIT_TYPES
			hit.dense_block++;
			#endif
			if (alldense) {
				dense_blocks_in_SB = l0;
				#if COUNT_HIT_TYPES
				hit.alldense_block++;
				#endif
			} else {
				//Only consider dense blocks
				__m256i const blocks_marked_dense_even = _mm256_cmpgt_epi16(blocks_even, sparse_limit_si256);
				__m256i const blocks_marked_dense_odd  = _mm256_cmpgt_epi16(blocks_odd,  sparse_limit_si256);
				uint32_t intmask_even = _mm256_movemask_epi8(blocks_marked_dense_even);
				uint32_t intmask_odd  = _mm256_movemask_epi8(blocks_marked_dense_odd);

				#if DEBUGINFO_SB_SELECT
				uint32_t intmask_even_debug = intmask_even;
				uint32_t intmask_odd_debug = intmask_odd;
				#endif

				uint32_t const even_shift_amount = 32 - ((l0+1)/2)*2;
				uint32_t const  odd_shift_amount = 32 - (l0/2)*2;

				//Only consider dense blocks coming *before* our block
				intmask_even = ((uint64_t)intmask_even) <<  even_shift_amount; //temporarily work on 64-bits to allow 32-bit-shifts, which would be UB for 32-bit integers
				intmask_odd  = ((uint64_t)intmask_odd)  <<  odd_shift_amount;

				#if DEBUGINFO_SB_SELECT
				printWord(intmask_even_debug,"intmask even");
				printWord(intmask_even);
				cout << even_shift_amount << " even_shift_amount " << endl;
				printWord(intmask_odd_debug, "intmask odd");
				printWord(intmask_odd);
				cout <<  odd_shift_amount << "  odd_shift_amount " << endl;
				#endif

				//Calculate an all-reduction to find the number of dense blocks before our block
				uint32_t const count_even = std::popcount(intmask_even)/2;
				uint32_t const count_odd  = std::popcount(intmask_odd)/2;
				dense_blocks_in_SB		  = l0>0 ? count_even + count_odd : 0;
			}
			uint64_t const dense_blocks_before_SB	 = L1L0[l1+DENSE_BLOCK_COUNT_POS]>>DENSE_BLOCK_COUNT_SHIFT;
			uint64_t const dense_blocks_before_block = dense_blocks_before_SB + dense_blocks_in_SB;
			uint64_t const dense_w0 = dense_blocks_before_block * W_IN_L0z;
			uint64_t w = dense_w0;

			#if REVERSE==2 && POP==2 && CUSTOM_L0!=512
			const bool reverse = (L0next - k) < (k - L0prev);
			#else
			constexpr bool reverse = false;
			#endif

			//treat k as signed now, because it can temporarily become negative in the coming loops
			int kk = k;

			#if DEBUGINFO_SB_SELECT
				cout << "dense blocks before SB       " << dense_blocks_before_SB<< endl;
				cout << "dense blocks before locally: " << dense_blocks_in_SB<< endl;
				cout << "dense blocks before total  : " << dense_blocks_before_block<< endl;
			#endif


			#if PREFETCH_NT
				//we do not want the raw bitvector bits in the cache, because they have very low information-density
				//thus we use the Non-Temporal Access Hint to tell cache system that it can quickly kick them out again
				//This way, more of the cache is reserved for high information-dense structures like the headers and the sample-tree
				#if CUSTOM_L0==512
				_mm_prefetch((char*)(&dense[w]), _MM_HINT_NTA);
				#elif CUSTOM_L0==1024
				_mm_prefetch((char*)(&dense[w]), _MM_HINT_NTA);
				_mm_prefetch((char*)(&dense[w+8]), _MM_HINT_NTA);
				#elif CUSTOM_L0==0 || CUSTOM_L0==2048
				_mm_prefetch((char*)(&dense[w]), _MM_HINT_NTA);
				_mm_prefetch((char*)(&dense[w+8]), _MM_HINT_NTA);
				_mm_prefetch((char*)(&dense[w+16]), _MM_HINT_NTA);
				_mm_prefetch((char*)(&dense[w+24]), _MM_HINT_NTA);
				#endif
			#endif

			if(!reverse) {
				#if POP==2 && CUSTOM_L0!=512
				kk -= L0prev;
				int pop1 = 0;
				int pop2 = 0;
				while(kk>0) [[likely]]{
					pop1 = std::popcount(dense[w]);
					pop2 = std::popcount(dense[w+1]);
					kk -= pop1 + pop2;
					w+=2;
				}
				w--;
				kk += pop2;
				if(kk<=0) {
					w--;
					kk += pop1;
				}
				#else
				kk -= L0prev;
				int pop0 = 0;
				while(kk>0) {
					pop0 = std::popcount(dense[w++]);
					kk-=pop0;
				}
				w--;
				kk += pop0;
				#endif
			}
			else {

				w += W_IN_L0z-1;
				kk -= L0next;
				int pop1 = 0;
				int pop2 = 0;
				while(kk<=0) [[likely]] {
					pop1 = std::popcount(dense[w]);
					pop2 = std::popcount(dense[w-1]);
					kk += pop1 + pop2;
					w-=2;
				}
				w++;
				int64_t const k2 = kk - pop2;
				if(k2>0) {
					w++;
					kk = k2;
				}
			}
			uint64_t const w_steps_forward = w - dense_w0;
			uint64_t const i_block_to_word = w_steps_forward * 64;
			//select into single word
			uint64_t const p = _pdep_u64(1ULL << (kk-1), dense[w]);
			uint64_t const i_in_word = _tzcnt_u64(p);
			return i_upto_block + i_block_to_word + i_in_word;
		}
		else {
			//Our block is sparse, need to know how many sparse bits there are in previous blocks in this SB
			#if COUNT_HIT_TYPES
			hit.sparse_block++;
			#endif
			uint64_t mbits_before_our_block = 0;
			if (allsparse) {
				#if COUNT_HIT_TYPES
				hit.allsparse_block++;
				#endif
				//lucky case: every single block in the SB is encoded sparse, so we can read off the number of previous mbits just from L0prev
				//up to computing the complement in case the mbits are zeros
				mbits_before_our_block = (minflag==1) ? L0prev : l0*L0_SIZEz - L0prev;
			}
			else {
				//Only consider sparse blocks
				__m256i const sparse_even_mask = _mm256_cmpgt_epi16(sparse_limit_plus1_si256,blocks_even);
				__m256i const sparse_odd_mask  = _mm256_cmpgt_epi16(sparse_limit_plus1_si256,blocks_odd);
				__m256i const sparse_even = _mm256_and_si256(blocks_even, sparse_even_mask);
				__m256i const sparse_odd  = _mm256_and_si256(blocks_odd, sparse_odd_mask);
				__m256i       sparse_packed = _mm256_packus_epi16(sparse_even, sparse_odd); //32 entries á 8 bit
				__m256i const l0_epi8 = _mm256_set1_epi8(l0);
				//Only consider sparse blocks *before* our block
				__m256i const packed_blocks_before_l0_mask = _mm256_cmpgt_epi8(l0_epi8, block_order_in_packus);
				sparse_packed = _mm256_and_si256(sparse_packed, packed_blocks_before_l0_mask); //kept only the blocks before our block
				//Calculate an all-reduction: The sum of m-bits in all sparse blocks before our block
				__m256i const sum_s3_s2_s1_0 = _mm256_sad_epu8(sparse_packed, _mm256_setzero_si256());
				__m256i const shifted_s3_s1  = _mm256_bsrli_epi128(sum_s3_s2_s1_0, 8);
				__m256i const sum_s32_s10    = _mm256_add_epi64(sum_s3_s2_s1_0, shifted_s3_s1);
				uint64_t const s32 = _mm256_extract_epi64(sum_s32_s10, 2);
				uint64_t const s10 = _mm256_extract_epi64(sum_s32_s10, 0);
				mbits_before_our_block = s32 + s10; //before this block

				#if DEBUGINFO_SB_SELECT
					print256_epi16(sparse_even_mask, "sparse even mask", 2, 0);
					print256_epi16(sparse_odd_mask, "sparse odd mask", 2, 1);
					print256_epi16(sparse_even, "sparse even", 2, 0);
					print256_epi16(sparse_odd, "sparse odd", 2,1);
				#endif
			}
			//The superblock knows where its sparse representation starts
			uint64_t const sparse_w0 = ((L1L0[l1+S_W0_HIGH_POS] << S_W0_LOWSIZE) | (L1L0[l1+S_W0_LOW_POS] & S_W0_LOWMASK)) & S_W0_TOTALMASK;
			k -= L0prev;

			#if DEBUGINFO_SB_SELECT
				printf("k'' %lu \n", k);
				printf("mbits before '' %lu \n", mbits_before_our_block);
				printf("sparse_w0    '' %lu \n", sparse_w0);
			#endif

			if (minflag==1) {
				#if COUNT_HIT_TYPES
				hit.sparse_match++;
				#endif
				//The block is stored in compressed form, and the deltas are stored for each 1-bit. We can directly access the k-th 1-bit.
				uint64_t const z = extract_z_value(&sparse[sparse_w0], mbits_before_our_block + k - 1);
				uint64_t const i_in_block = z + k-1 ;
				#if DEBUGINFO_SB_SELECT
					printf("z %lu \n", z);
					printf("i_in_block '' %lu \n", i_in_block);
				#endif
				return i_upto_block + i_in_block;
			}
			else {
				#if COUNT_HIT_TYPES
				hit.sparse_nomatch++;
				#endif
				//The block is stored in compressed form, but the deltas represent 0-bits. (because those are in the minority)
				//We need to scan through the 0-bit deltas to find our spot.
				#if DEBUGINFO_SB_SELECT
				cout << "mbits in block = " << mbits_in_block << endl;
				#endif
				// if (Z_SEARCH_STRAT==Z_LINEAR_SEARCH) {

					constexpr uint64_t STRIDE = 10;
					uint64_t z_index=0;
					uint64_t z = 0;
					//Linear Scan in strides of 10
					while (z_index + STRIDE < mbits_in_block) {
						z = extract_z_value(&sparse[sparse_w0], mbits_before_our_block + z_index + STRIDE);
						if (z >= k) {
							break;
						}
						z_index += STRIDE;
					}
					//Linear Scan in steps of 1
					while (z_index < mbits_in_block) {
						z = extract_z_value(&sparse[sparse_w0], mbits_before_our_block + z_index);
						if (z >= k) {
							break;
						}
						z_index++;
					}
					uint64_t const i_in_block = k-1 + z_index;
					return i_upto_block + i_in_block;
				// }
				// else if (Z_SEARCH_STRAT==Z_BINARY_SEARCH) {
					// uint64_t low = 0, high = mbits_in_block;
					// while (low < high) {
						// uint64_t mid = (low + high) / 2;
						// uint64_t z = extract_z_value(&sparse[sparse_w0], mbits_before_our_block + mid);
						// if (z < k) {
							// low = mid + 1;
						// } else {
							// high = mid;
						// }
					// }
					// uint64_t z_index = low;
					// uint64_t const i_in_block = k-1 + z_index;
					// return i_upto_block + i_in_block;

			}
		}
	}


	[[nodiscard]] uint64_t select_0(int64_t k) const override { return 0;}


	static inline uint64_t extract_12bit_value(const aligned_vector &vec, const uint64_t start_offset_bits, uint64_t k) {
		//interpreting the 64-bit vector as array of bytes allows more precise access to the 12 bits of interest
		const auto *vec_bytes = reinterpret_cast<const uint8_t*>(vec.data());
		uint64_t const bit_pos = k * 12 + start_offset_bits;
		uint64_t const byte_pos = bit_pos / 8;
		uint64_t const bit_offset = bit_pos % 8;
		uint64_t const value = (uint64_t)vec_bytes[byte_pos] | ((uint64_t)vec_bytes[byte_pos + 1] << 8);
		return (value >> bit_offset) & 0xFFF;
	}



	uint64_t rank_1(const uint64_t i) const override {
		uint64_t const l1n = i/L1_SIZEz;
		uint64_t const l1  = l1n * WORDS_IN_L1HEADER;

		uint64_t const pos_in_L0 =  i % L0_SIZEz;

		#if REVERSE==2 && POP==2 && CUSTOM_LO != 512
		bool const reverse = pos_in_L0 > L0_SIZEz/2;
		#else
		constexpr bool reverse = false;
		#endif

		#if PREFETCH_T0
		_mm_prefetch((char*)(&L1L0[l1]), _MM_HINT_T0);
		#endif
		uint64_t const L1  = L1L0[l1+L1_POS] >> L1_SHIFT;
		bool const minflag = (L1L0[l1+MINFLAG_POS]>>MINFLAG_SHIFT)&1ULL;


		#if DEBUGINFO_RANK
		cout << "l1n = " << l1n << endl;
		cout << "l1	= " << l1 << endl;
		cout << "minflag = " << minflag << endl;
		cout << "reverse = " << reverse << endl;
		#endif


		// const bool last_L0_is_sparse = L1L0[l1+LAST_L0_IS_SPARSE_W]>>LAST_L0_IS_SPARSE_SHIFT & 1ULL;

		//For comments on the SIMD instructions, see select_1
		__m128i A_raw = _mm_loadu_si128((__m128i*)&L1L0[l1 + L0A_POS]);
		__m128i b_raw = _mm_loadu_si128((__m128i*)&L1L0[l1 + L0b_POS]);
		__m128i c_raw = _mm_loadu_si128((__m128i*)((uint8_t*)&L1L0[l1 + L0c_POS] + 4));
		__m128i d_raw = _mm_loadu_si128((__m128i*)&L1L0[l1 + L0d_POS]);

		__m128i	    A = _mm_and_si128(A_raw, remove_clutter_around_A);

		__m128i	b_low = _mm_shuffle_epi8(b_raw,unpack12to16_shuff_low);
		__m128i	c_low = _mm_shuffle_epi8(c_raw,unpack12to16_shuff_low);
		__m128i	d_low = _mm_shuffle_epi8(d_raw,unpack12to16_shuff_low);

		__m128i	b_hi  = _mm_shuffle_epi8(b_raw,unpack12to16_shuff_hi);
		__m128i	c_hi  = _mm_shuffle_epi8(c_raw,unpack12to16_shuff_hi);
		__m128i	d_hi  = _mm_shuffle_epi8(d_raw,unpack12to16_shuff_hi);

				b_low = _mm_and_si128(b_low, keepLower12_epi16);
				c_low = _mm_and_si128(c_low, keepLower12_epi16);
				d_low = _mm_and_si128(d_low, keepLower12_epi16);

				b_hi  = _mm_srli_epi16(b_hi,4);
				c_hi  = _mm_srli_epi16(c_hi,4);
				d_hi  = _mm_srli_epi16(d_hi,4);

		__m128i    b  = _mm_or_si128(b_hi, b_low);
		__m128i    c  = _mm_or_si128(c_hi, c_low);
		__m128i    d  = _mm_or_si128(d_hi, d_low);

		__m128i	   B  = _mm_add_epi16(A,b);
		__m128i	   C  = _mm_add_epi16(B,c);
		__m128i	   D  = _mm_add_epi16(C,d);


		if(minflag==0) [[unlikely]] {
			A = _mm_sub_epi16(MAX_A, A);
			B = _mm_sub_epi16(MAX_B, B);
			C = _mm_sub_epi16(MAX_C, C);
			D = _mm_sub_epi16(MAX_D, D);
		}

		//maybe need to shift A such that its zero-entry comes into the first array slot
		//A = _mm_srli_si128(A,2);

		alignas(32) uint16_t tmp_L0s[32];
		_mm_store_si128((__m128i*)&tmp_L0s[0], A);
		_mm_store_si128((__m128i*)&tmp_L0s[8], B);
		_mm_store_si128((__m128i*)&tmp_L0s[16], C);
		_mm_store_si128((__m128i*)&tmp_L0s[24], D);

		const uint64_t l0 = (i%L1_SIZEz)/L0_SIZEz; //(0..31)
		const uint64_t l0_next = (l0 + 1); //(1..32)
		const uint64_t l0_section        = (l0%4)*8;
		const uint64_t l0_next_section   = (l0_next%4)*8;
		const uint64_t l0_nr_in_section       = l0/4;
		const uint64_t l0_next_nr_in_section  = l0_next/4;

		uint64_t L0_prev = tmp_L0s[l0_section + l0_nr_in_section];
		uint64_t L0_next = tmp_L0s[l0_next_section + l0_next_nr_in_section];


		// if (l0==0 && L0_prev!=0) {
			// cout << "error, L0 expected to be 0, but is " << L0_prev << endl;
		// }

		if (l0_next==32) {
			uint64_t l1_next = l1 + WORDS_IN_L1HEADER;
			uint64_t const L1next = L1L0[l1_next + L1_POS] >> L1_SHIFT;
			L0_next = L1next - L1;
		}



		#if DEBUGINFO_RANK
		cout << "l0 " << l0 << endl;
		cout <<	"L0prev " << L0_prev << endl;
		cout << "l0_next " << l0_next << endl;
		cout << "L0next " << L0_next << endl;
		#endif

		bool block_is_dense;
		bool alldense;
		bool allsparse;
		uint64_t mbits_in_block;
		__m256i blocks_even;
		__m256i blocks_odd;

		if (BV_COMPRESSION==1) {

			uint64_t const onebits_in_block	= L0_next - L0_prev;
			mbits_in_block   = minflag==1 ? onebits_in_block : L0_SIZEz - onebits_in_block;
			block_is_dense   = mbits_in_block > MAX_SPARSE;

			#if CHECK_ALL_DENSE_SPARSE
			alldense  = L1L0[l1+ALLDENSE_W]>>ALLDENSE_SHIFT & 1ULL;
			allsparse = L1L0[l1+ALLSPARSE_W]>>ALLSPARSE_SHIFT & 1ULL;
			#else
			alldense = false;
			allsparse = false;
			#endif

			//block values a,b,c,d count how many m-bits there are in an individual L0-block
			__m128i blocks_even_lo = _mm_unpacklo_epi16(b, d);
			__m128i blocks_even_hi = _mm_unpackhi_epi16(b, d);
					blocks_even    = _mm256_set_m128i(blocks_even_hi, blocks_even_lo);
			A = _mm_srli_si128(A,2); //need to align A and D to now subtract backwards
			__m128i             a  = _mm_sub_epi16(A,D);  //number of 1-bits in each A-block individually (not cumulative)
			if (minflag==0) [[unlikely]] {
				a = _mm_sub_epi16(mm_L0size_epi16, a); //convert to 0-bit count
			}
			__m128i blocks_odd_lo = _mm_unpacklo_epi16(c, a);
			__m128i blocks_odd_hi = _mm_unpackhi_epi16(c, a);
					blocks_odd    = _mm256_set_m128i(blocks_odd_hi, blocks_odd_lo);
		}
		else {
			block_is_dense = true;
			alldense = true; //In case of no compression we know that all blocks are stored dense (uncompressed)
			allsparse = false;
			mbits_in_block = 0;
			// blocks_even = _mm256_setzero_si256();
			// blocks_odd = _mm256_setzero_si256();
		}



		uint64_t  L0_tracked = reverse ? L0_next : L0_prev;

		#if DEBUGINFO_RANK
		cout << "L0_tracked " << L0_tracked << endl;
		#endif

		if (block_is_dense) {
			//Our block is dense, need to know how many previous blocks in this SB are also dense
			uint64_t prev_dense_blocks_in_SB = 0;
			#if COUNT_HIT_TYPES
			hit.dense_block++;
			#endif
			if (alldense) {
				prev_dense_blocks_in_SB = l0;
				#if COUNT_HIT_TYPES
				hit.alldense_block++;
				#endif
			} else {
				//Only consider dense blocks
				__m256i const blocks_marked_dense_even = _mm256_cmpgt_epi16(blocks_even, sparse_limit_si256);
				__m256i const blocks_marked_dense_odd  = _mm256_cmpgt_epi16(blocks_odd,  sparse_limit_si256);
				uint32_t intmask_even = _mm256_movemask_epi8(blocks_marked_dense_even);
				uint32_t intmask_odd  = _mm256_movemask_epi8(blocks_marked_dense_odd);

				uint32_t const even_shift_amount = 32 - ((l0+1)/2)*2;
				uint32_t const  odd_shift_amount = 32 - (l0/2)*2;

				//Only consider dense blocks coming *before* our block
				intmask_even = ((uint64_t)intmask_even) <<  even_shift_amount; //temporarily work on 64-bits to allow 32-bit-shifts, which would be UB for 32-bit integers
				intmask_odd  = ((uint64_t)intmask_odd)  <<  odd_shift_amount;

				//Calculate an all-reduction to find the number of dense blocks before our block
				uint32_t const count_even = std::popcount(intmask_even)/2;
				uint32_t const count_odd  = std::popcount(intmask_odd)/2;
				prev_dense_blocks_in_SB		  = l0>0 ? count_even + count_odd : 0;

				#if DEBUGINFO_RANK
				printWord(intmask_even);
				cout << even_shift_amount << " even_shift_amount " << endl;
				printWord(intmask_odd);
				cout <<  odd_shift_amount << "  odd_shift_amount " << endl;
				cout << "prev_dense_blocks_in_SB: " << prev_dense_blocks_in_SB << endl;
				#endif

			}
			uint64_t const dense_blocks_before_SB	 = L1L0[l1+DENSE_BLOCK_COUNT_POS]>>DENSE_BLOCK_COUNT_SHIFT;
			uint64_t const dense_blocks_before_block = dense_blocks_before_SB + prev_dense_blocks_in_SB;
			uint64_t const dense_w0 = dense_blocks_before_block * W_IN_L0z;
			uint64_t w_border = reverse ? dense_w0 + (W_IN_L0z -1) : dense_w0;

			#if DEBUGINFO_RANK
			cout << "dense blocks before SB       " << dense_blocks_before_SB<< endl;
			cout << "dense blocks before locally: " << dense_blocks_before_block<< endl;
			cout << "dense blocks before total  : " << dense_blocks_before_block<< endl;
			cout << "dense_w0: " << dense_w0 << endl;
			cout << "w_border: " << w_border << endl;
			#endif


			#if PREFETCH_NT
			#if CUSTOM_L0==512
						_mm_prefetch((char*)(&dense[dense_w0]), _MM_HINT_NTA);
			#elif CUSTOM_L0==1024
						_mm_prefetch((char*)(&dense[dense_w0]), _MM_HINT_NTA);
						_mm_prefetch((char*)(&dense[dense_w0+8]), _MM_HINT_NTA);
			#elif CUSTOM_L0==0 || CUSTOM_L0==2048
						_mm_prefetch((char*)(&dense[dense_w0]), _MM_HINT_NTA);
						_mm_prefetch((char*)(&dense[dense_w0+8]), _MM_HINT_NTA);
						_mm_prefetch((char*)(&dense[dense_w0+16]), _MM_HINT_NTA);
						_mm_prefetch((char*)(&dense[dense_w0+24]), _MM_HINT_NTA);
			#endif
			#endif

			//Linear popcounting from L0 border to target word Wt

			uint64_t const w_target = dense_w0 + (i % L0_SIZEz)/64;
			uint64_t const pos_target = i % 64;

			#if DEBUGINFO_RANK
			cout << "w_target: " << w_target << endl;
			cout << "pos_target: " << pos_target << endl;
			cout << "L0_tracked=" << L0_tracked << endl;
			#endif

			uint64_t const target_W = dense[w_target];
			const bool w_target_is_odd = w_target%2;
			int pop1 = 0;
			int pop2 = 0;
			if(reverse) {
				for(uint64_t v=w_border; v>w_target; v-=2) {
					pop1 = std::popcount(dense[v]);
					pop2 = std::popcount(dense[v-1]);
					L0_tracked -= pop1 + pop2;
				}
				if(!w_target_is_odd) {
					L0_tracked += pop2;
				}
				L0_tracked -= std::popcount(target_W >> pos_target);
			}
			else { //forwards scanning
				for(uint64_t v=w_border; v<w_target; v+=2) {
					pop1 = std::popcount(dense[v]);
					pop2 = std::popcount(dense[v+1]);
					L0_tracked += pop1 + pop2;
					#if DEBUGINFO_RANK
					cout << "v=" << v<<", L0_tracked=" << L0_tracked << endl;
					#endif
				}
				if(w_target_is_odd) {
					L0_tracked -= pop2;
					#if DEBUGINFO_RANK
					cout << "target odd, L0_tracked=" << L0_tracked << endl;
					#endif
				}
				if(pos_target>0) {
					L0_tracked += std::popcount(target_W << (64-pos_target));
					#if DEBUGINFO_RANK
					printWord(target_W, "target_W");
					printWord(target_W << (64-pos_target), "target_W << (64-pos_target)");
					cout << "partial popcnt, L0_tracked=" << L0_tracked << endl;
					#endif
				}
			}
			return L1 + L0_tracked;
		} else {
			//Our block is sparse. need to know how many sparse bits there are in previous blocks in this SB
			#if COUNT_HIT_TYPES
			hit.sparse_block++;
			#endif
			uint64_t mbits_before_our_block = 0;
			if (allsparse) {
				#if COUNT_HIT_TYPES
				hit.allsparse_block++;
				#endif
				//lucky case: every single block in the SB is encoded sparse, so we can read off the number of previous mbits just from L0prev
				//up to computing the complement in case the mbits are zeros
				mbits_before_our_block = (minflag==1) ? L0_prev : l0*L0_SIZEz - L0_prev;
			}
			else {
				//Only consider sparse blocks
				__m256i const sparse_even_mask = _mm256_cmpgt_epi16(sparse_limit_plus1_si256,blocks_even);
				__m256i const sparse_odd_mask  = _mm256_cmpgt_epi16(sparse_limit_plus1_si256,blocks_odd);
				__m256i const sparse_even = _mm256_and_si256(blocks_even, sparse_even_mask);
				__m256i const sparse_odd  = _mm256_and_si256(blocks_odd, sparse_odd_mask);
				__m256i       sparse_packed = _mm256_packus_epi16(sparse_even, sparse_odd); //32 entries á 8 bit
				__m256i const l0_epi8 = _mm256_set1_epi8(l0);
				//Only consider sparse blocks *before* our block
				__m256i const packed_blocks_before_l0_mask = _mm256_cmpgt_epi8(l0_epi8, block_order_in_packus);
				sparse_packed = _mm256_and_si256(sparse_packed, packed_blocks_before_l0_mask); //kept only the blocks before our block
				//Calculate an all-reduction: The sum of m-bits in all sparse blocks before our block
				__m256i const sum_s3_s2_s1_0 = _mm256_sad_epu8(sparse_packed, _mm256_setzero_si256());
				__m256i const shifted_s3_s1  = _mm256_bsrli_epi128(sum_s3_s2_s1_0, 8);
				__m256i const sum_s32_s10    = _mm256_add_epi64(sum_s3_s2_s1_0, shifted_s3_s1);
				uint64_t const s32 = _mm256_extract_epi64(sum_s32_s10, 2);
				uint64_t const s10 = _mm256_extract_epi64(sum_s32_s10, 0);
				mbits_before_our_block = s32 + s10; //before this block

				#if DEBUGINFO_RANK
					// print256_epi16(sparse_even_mask, "sparse even mask", 2, 0);
					// print256_epi16(sparse_odd_mask, "sparse odd mask", 2, 1);
					// print256_epi16(sparse_even, "sparse even", 2, 0);
					// print256_epi16(sparse_odd, "sparse odd", 2,1);
				#endif
			}
			//The superblock knows where its sparse representation starts
			uint64_t const sparse_w0 = ((L1L0[l1+S_W0_HIGH_POS] << S_W0_LOWSIZE) | (L1L0[l1+S_W0_LOW_POS] & S_W0_LOWMASK)) & S_W0_TOTALMASK;

			#if DEBUGINFO_RANK
				printf("mbits before '' %lu \n", mbits_before_our_block);
				printf("sparse_w0    '' %lu \n", sparse_w0);
			#endif
			#if DEBUGINFO_RANK
			cout << "mbits in block = " << mbits_in_block << endl;
			#endif
			constexpr uint64_t STRIDE = 10;
			uint64_t z_count=0;
			uint64_t z = 0;
			uint64_t pos_probed = 0;
			//Scan in strides of 10
			while (z_count + STRIDE < mbits_in_block & pos_probed < pos_in_L0) {
				z = extract_z_value(&sparse[sparse_w0], mbits_before_our_block + z_count + STRIDE);
				pos_probed = z + z_count + STRIDE;
				#if DEBUGINFO_RANK
				cout << (z_count+STRIDE) << ": z=" << z << ", i_probed="<<pos_probed <<endl;
				#endif
				z_count +=STRIDE;
			}
			if (pos_probed>=pos_in_L0 & pos_in_L0>0) {
				z_count-=STRIDE;
			}
			pos_probed = 0;
			//Scan in steps of 1
			while (z_count < mbits_in_block & pos_probed < pos_in_L0) {
				z = extract_z_value(&sparse[sparse_w0], mbits_before_our_block + z_count);
				pos_probed = z + z_count;
				#if DEBUGINFO_RANK
				cout << z_count << ": z=" << z << ", i_probed="<<pos_probed <<endl;
				#endif
				z_count++;
			}
			if (pos_probed>=pos_in_L0 & pos_in_L0>0) {
				z_count--;
			}
			uint64_t minority_bits = z_count;
			uint64_t majority_bits = pos_in_L0 - minority_bits;
			uint64_t rank_in_L0 = (minflag==0) ? majority_bits : minority_bits;
			#if DEBUGINFO_RANK
			cout << "after loop:" << endl;
			cout << "minority_bits=" << minority_bits << endl;
			cout << "majority_bits=" << majority_bits << endl;
			cout << "pos_in_L0=" << pos_in_L0 << endl;
			cout << "rank_in_L0=" << rank_in_L0 << endl;
			#endif

			return L1 + L0_prev + rank_in_L0;
		}

	}

};
















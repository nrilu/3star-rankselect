//
// Created by nicco on 7/30/24.
//

#pragma once

#include "alignedallocator.hpp"
#include <cstdint>
#include <vector>
#include <chrono>

#include <x86intrin.h>  //_pdep, _tzcnt, __mm256i
#include <bit>          //std::popcount

//#include <immintrin.h>

//Structure constants
//static constexpr size_t HUGEPAGE_SIZE = 1<<21; //2MiB
static constexpr uint64_t WORDSIZE = 64;  //unsigned version
static constexpr int64_t  WORDSIZE_s = 64; //signed version, might be faster when adding to other signed values
static constexpr uint64_t L2SIZE = 1024;
static constexpr uint64_t L2_IN_L1 = 8;
static constexpr uint64_t WORDS_IN_L2 = L2SIZE / WORDSIZE; //16
static constexpr uint64_t L1SIZE = L2SIZE * L2_IN_L1; //8192
static constexpr uint64_t WORDS_IN_L1 = WORDS_IN_L2 * L2_IN_L1; //128
static constexpr uint64_t WORDS_IN_MERGED = 130;
extern uint64_t SAMPLEDIST_EXP;   //used for SAMPLEDIST = 2^ SAMPLEDIST_EXP
extern uint64_t SAMPLEDIST;
//constexpr uint64_t SAMPLEDIST_EXP = 15;


static constexpr uint64_t MIN_SPACING_FOR_PLA = 1000;

class Bitvector {
public:
    virtual ~Bitvector() = default;
    virtual void build_auxiliaries() = 0;
    //Total space usage of the data-structure
    virtual uint64_t space_in_bits() const = 0;
    virtual bool access(uint64_t i) const = 0;
    virtual uint64_t rank_1(uint64_t i) const = 0;
    virtual uint64_t rank_0(uint64_t i) const = 0;
    virtual uint64_t select_1(int64_t k) const = 0;
    virtual uint64_t select_0(int64_t k) const = 0;
};

std::string get_unixtimestamp();

//printing toggles
extern bool show_overhead;

//global print functions
static void print_overhead_line(const std::string &label, int64_t bits_used);
static void print_separator();

//choice of Job
extern std::string BITGEN;
extern std::string BVNAME;
extern bool doBench;
extern bool doValid;
extern bool doValidw;
extern bool doMultibench;
constexpr uint64_t validw_updaterate = 8192;

//choice of bits
extern uint64_t synth_bits;
extern bool useSynth;
extern bool do_manual_bits;
extern bool do_fractals;
extern bool do_chunky;
extern bool do_smooth;
extern bool do_quicksmooth;
extern bool do_gidneysmooth;
extern bool do_alternate;
extern bool do_every_kth;
extern bool do_bimodal;
extern int  bimodal_factor;
extern bool do_patches;
extern bool do_sinus;
extern bool do_gap;
//extern double bimodal_short;
//extern double bimodal_long;

extern int every_kth_spacing;
extern int synth_seqsize;
extern int synth_01ratio;
extern bool invert_all_bits;

//Choice of queries
extern bool do_random_queries;
extern uint64_t synth_accesses;
extern uint64_t synth_ranks;
extern uint64_t synth_selects;
extern uint64_t synth_select1s;
extern uint64_t synth_select0s;
extern std::vector<uint64_t> query_type_counter;

//Input counters
extern uint64_t N_queries;
extern uint64_t N_bits;
extern uint64_t N_ones;
extern uint64_t N_zeros;
extern uint64_t N_words;
extern uint64_t seed;
extern double eff_01ratio;

//Input information for m3 benchmarking
extern int parameter_cycles;
extern int curr_parameter_cycle;
extern int curr_instance;
extern double sin_threshhold;

//For the gap gadget instance
inline extern uint64_t curr_gap_size_bits=500;
constexpr uint64_t generated_dummy_queries = 1 << 26;

//input storage
extern std::vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>> bits;
extern std::vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>> queries;

//Bits stored in reversed 64-bit order to comply with Sux and Poppy standard
extern std::vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>> bits_64bit_reversed;
extern uint64_t total_used_space_in_bits;

//timings
typedef std::chrono::steady_clock::time_point tp;
typedef std::chrono::high_resolution_clock::time_point highres_tp;
extern tp tCopyFinished;
extern double dRead;
extern double dCopy;
extern double dBuild;
extern double dQueries;
extern double dPerQuery;


//Tree3 construction strategies
#define GET_CUBIC_THEORY_PARAMS 1
#define GET_LIFTED_CUBIC_THEORY_PARAMS 2
#define GET_TFOCUS_THEORY_PARAMS 3
#define GET_TOPHEAVY_PARAMS 4
#define GET_SMALLEST_TREE_PARAMS 5
#define GET_FASTEST_TREE_PARAMS 6

#define Z_LINEAR_SEARCH 1
#define Z_BINARY_SEARCH 2

inline extern int TREE3_STRATEGY = GET_TOPHEAVY_PARAMS;
inline extern int ALPHA = 16;
inline extern int SUMMARY_LEVELS = 1;
inline extern int Z_SEARCH_STRAT=Z_LINEAR_SEARCH;
inline extern int BV_COMPRESSION=1;

//pass information
struct tree_stats {
	double theo_a_exp;
	double theo_b_exp;
	uint64_t a_exp;
	uint64_t b_exp;
	uint64_t topgroups, midgroups, botgroups;
	uint64_t theo_worstcase_size_bits;
	uint64_t theo_actual_size_bits;
	uint64_t final_tree_size_bits;
	double theo_worstcase_size_factor;
	double theo_actual_size_factor;
	double final_tree_size_factor;
};

//runtime counters
struct hit_counter {
	uint64_t dense_block = 0;
	uint64_t alldense_block = 0;
	uint64_t sparse_block = 0;
	uint64_t allsparse_block = 0;
	uint64_t sparse_match = 0;
	uint64_t sparse_nomatch = 0;
};
struct cpu_cycles {
	uint64_t in_tree = 0;
	uint64_t scan_to_L1 = 0;
	uint64_t in_SIMD = 0;
	uint64_t popcounting = 0;
	uint64_t zreads = 0;
};
inline hit_counter hit{};
inline cpu_cycles cycles{};
inline tree_stats treestats{};


#include "bv.hpp"
#include "file_io.hpp"
#include "synth.hpp"
#include "naive.hpp"
#include "pastaflat.hpp"
#include "m3.hpp"
#include "Tree3.hpp"

#include "external_sdsl_sd.hpp"
#include "external_sdsl_rrr.hpp"
#include "external_sdsl_mcl.hpp"
#include "external_sdsl_v1.hpp"
#include "external_sdsl_v5.hpp"
#include "external_rank9sel.hpp"
#include "external_simple_s.hpp"
#include "external_poppy.hpp"
#include "external_poppy2.hpp"


#include <iostream>
#include <string>
#include <cpuid.h>
#include <fstream>

using namespace std;
//------------------------------MANUAL PARAMETERS---------------------------------------//
bool run_testsection=0;
bool show_overhead  =0;
bool even_more_silent = EVEN_MORE_SILENT;
std::vector<std::pair<std::string, std::string>> CURR_RESULT_LINE;

const std::vector<std::string> names =
    {"m2","m3","m5","m7","naive",
    "minflag","tight","merged","tight_mm","tight_sjump", "minflag_ssimd","mincount_simd", "mincount2_tests",
    "pastaflat",
    "sdsl_sd", "sdsl_rrr", "sdsl_mcl", "sdsl_v1", "sdsl_v5",
    "rank9sel", "simple_s0","simple_s1","simple_s2","simple_s3", "simple_s4", "poppy", "poppy2", "s18vec"
    };


vector<string> BVNAMES_LIST = {"simple_s0","simple_s1","simple_s2","simple_s3","simple_s4"};//a name from the list
string JOB          = "multibench";   //multibench, valid, validw                 :: Bench or validate the structure
string BITGEN       = "indep"; //synth, file, man, alternate, fractal, chunky, smooth, indep, every_kth, bimodal, patches:: How input should be generated
string QTYPE        = "random";     //man, random                                 :: How queries should be generated

std::vector<double> N_BITS = {1e9};

vector<int> synth_01ratios = {128};
vector<int> TREE3_STRATEGY_LIST = {GET_TOPHEAVY_PARAMS};
vector<int> ALPHA_LIST = {ALPHA};
vector<int> SUMMARY_LEVEL_LIST = {SUMMARY_LEVELS};
vector<int> Z_SEARCH_STRAT_LIST = {Z_SEARCH_STRAT};
vector<int> BV_COMPRESSION_LIST = {BV_COMPRESSION};
bool invert_all_bits=0;

int iterations = 1; //no. of iterations per single instance. Higher = more correct timings
int instances = 1; //no. of instances to generate per parameter setting. Higher = more diverse inputs

//Random Queries
uint64_t synth_accesses = 0;
uint64_t synth_ranks    = 0;
uint64_t synth_selects  = 0;
uint64_t synth_select1s = 2e7;
uint64_t synth_select0s = 0;


//number of sequential identical 64 bits (only workes with 50:50 ratio)
vector<int> synth_seqsizes = {0};
int every_kth_spacing = 3;
int bimodal_factor = 500;
//double bimodal_short = 0.25;
//double bimodal_long = 20;

//For the gap gadget instance
vector<uint64_t> GAP_SIZES = {curr_gap_size_bits};

constexpr int GAP_CRITICAL_QUERIES = 20;
constexpr int GAP_DUMMY_QUERIES = 2e7;

constexpr uint64_t POPPY_MAX_N = 2e9;
//Old sample distances when not using tree
std::vector<uint64_t> SAMPLEDIST_EXPS = {19}; //only for m2/m5/m7 when without tree. In tree gets replaced by TREE_EXP_LVL0

//-----------INTERNAL VARIABLES----------------------------------------------------------//


//forward declarations
int shown_bit_preview_rows = 30;

//Bitvector config
string BVNAME;
uint64_t SAMPLEDIST_EXP;
uint64_t SAMPLEDIST;

//job selection
bool doBench;
bool doValid;
bool doValidw;
bool doMultibench;

int curr_parameter_cycle;
int curr_instance;
int parameter_cycles = 16; //no of repetitions of the same instance in case we benchmark different Tree3 configurations

// int extra_stretch_bits_iterations = 1;
// int curr_extra_stretch_bits = 0;
// constexpr uint64 EXTRA_STRECH_BITS_STEP = 1000;

double gap_avg_dummy=0;
double gap_critical_1st=0;
double gap_critical_2nd=0;
double gap_critical_3rd=0;
double gap_critical_avg_cached=0;
double gap_single_dense_1st=0;
double gap_single_dense_2nd=0;
double gap_single_dense_3rd=0;
double gap_single_dense_avg_cached=0;

//bit generation
uint64_t synth_bits;
bool useSynth;
bool do_manual_bits;
bool do_fractals;
bool do_chunky;
bool do_smooth;
bool do_quicksmooth;
bool do_gidneysmooth;
bool do_alternate;
bool do_every_kth;
bool do_bimodal;
bool do_patches;
bool do_sinus;
bool do_gap;
int synth_seqsize;
int synth_01ratio;
uint64_t seed;
double eff_01ratio;
uint64_t N_bits;
uint64_t N_ones;
uint64_t N_zeros;
uint64_t N_words;
bool shadow_current_result_line;
uint64_t total_used_space_in_bits;

//query generation
uint64_t N_queries;
bool do_random_queries;
vector<uint64_t> query_type_counter;

//storing input and queries
vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>> bits;
vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>> queries;

//Reversed bitvector for Poppy
vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>> bits_64bit_reversed;

//timing
tp tCopyFinished;
double dRead;
double dCopy;
double dBuild;
double dQueries;
double dPerQuery;
double averaged_avg_query;
double averaged_build;
double averaged_iteration;
uint64_t used_space;

//data structure space usage
double current_overhead_percent;

void init_counters() {
    dRead      = 0;
    dCopy      = 0;
    dBuild     = 0;
    dQueries   = 0;
    dPerQuery  = 0;
    averaged_avg_query  =0;
    averaged_build      =0;
    averaged_iteration  =0;
    N_bits   =0;
    N_queries=0;
    N_ones   =0;
    N_zeros  =0;
    N_words  =0;
    used_space=0;
    seed = 0;
    curr_parameter_cycle=0;
    curr_instance=0;
}

void init_benchvalues() {
    BVNAME          = BVNAMES_LIST[0];
    synth_bits      = N_BITS[0];
    synth_seqsize   = synth_seqsizes[0];
    synth_01ratio   = synth_01ratios[0];
    #if TREE_TYPE
    //if using BTREE, the array SAMPLEDIST_EXPS is not used
    //keep one dummy entry such that loop of multibench just passes by
    SAMPLEDIST_EXPS.clear();
    SAMPLEDIST_EXPS.push_back(TREE_EXP_LVL0);
    #endif
}

void init_selections() {
    doBench     = JOB == "bench";
    doValid     = JOB == "valid" || JOB == "validw";
    doValidw    = JOB == "validw";
    doMultibench = JOB == "multibench";
    do_manual_bits = BITGEN == "man";
    do_chunky      = BITGEN == "chunky";
    do_smooth      = BITGEN == "smooth";
    do_quicksmooth = BITGEN == "quicksmooth";
    do_gidneysmooth= BITGEN == "indep";
    do_fractals    = BITGEN == "fractal";
    do_alternate   = BITGEN == "alternate";
    do_every_kth   = BITGEN == "every_kth";
    do_bimodal     = BITGEN == "bimodal";
    do_patches     = BITGEN == "patches";
    do_sinus       = BITGEN == "sin";
    do_gap         = BITGEN == "gap";
    useSynth       = BITGEN == "synth" || do_manual_bits || do_fractals || do_chunky || do_smooth \
                                       || do_alternate || do_quicksmooth || do_gidneysmooth || do_every_kth || do_bimodal \
                                        || do_patches || do_sinus || do_gap;
    do_random_queries = QTYPE == "random";
}

void set_sample_exp(uint64_t const my_sample_exp) {
    #define POW2(x) (1 << (x))
    #if TREE_TYPE
        SAMPLEDIST_EXP  = TREE_EXP_LVL0;
        SAMPLEDIST      = POW2(SAMPLEDIST_EXP);
    #elif CMAKE_SAMPL_EXP>0
        SAMPLEDIST_EXP  = CMAKE_SAMPL_EXP;
        SAMPLEDIST      = POW2(CMAKE_SAMPL_EXP);
    #else
        SAMPLEDIST_EXP  = my_sample_exp;
        SAMPLEDIST      = std::pow(2,SAMPLEDIST_EXP);
    #endif
}

string in_path;
string result_path;
bool write_individual_query_results;

string mb_write_folder = "multibench_out/";
string mb_name;
string mb_path;

void init_paths() {
    string read_folder = "../generated/";
    string read_filename = "extremes";
    write_individual_query_results = false;
    in_path = read_folder + read_filename + ".txt";
    result_path = read_folder + read_filename + "_out.txt";
}

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

using std::cout;
using std::endl;

//--------------FUNCTIONS------------------------------------------------------------//

std::unique_ptr<Bitvector> createBitvector(const std::string &type) {
    if (type=="m3")            return std::make_unique<m3>();
    if (type=="naive")         return std::make_unique<naive>();
    #if LOAD_PASTA
    if (type=="pastaflat")     return std::make_unique<pastaflat>();
    #endif
    #if LOAD_SDSL
    if (type=="sdsl_sd")       return std::make_unique<external_sdsl_sd>();
    if (type=="sdsl_rrr")      return std::make_unique<external_sdsl_rrr>();
    if (type=="sdsl_mcl")      return std::make_unique<external_sdsl_mcl>();
    if (type=="sdsl_v5")       return std::make_unique<external_sdsl_v5>();
    if (type=="sdsl_v1")       return std::make_unique<external_sdsl_v1>();
    #endif
    #if LOAD_SUX
    if (type=="rank9sel")      return std::make_unique<external_rank9sel>();
    if (type=="simple_s0")     return std::make_unique<external_simple_select<0>>();
    if (type=="simple_s1")     return std::make_unique<external_simple_select<1>>();
    if (type=="simple_s2")     return std::make_unique<external_simple_select<2>>();
    if (type=="simple_s3")     return std::make_unique<external_simple_select<3>>();
    if (type=="simple_s4")     return std::make_unique<external_simple_select<4>>();
    #endif
    #if LOAD_POPPY
    if (type=="poppy")         return std::make_unique<external_poppy>();
    if (type=="poppy2")        return std::make_unique<external_poppy2>();
    #endif

    throw std::runtime_error("BVNAME is set to "+type+" but its external library is not loaded");
}

double get_total_overhead_in_percent(int64_t total_space_in_bits) {
    int64_t bits_overhead = total_space_in_bits - N_bits; //N_bits == #bits in the input
    double ovhd_percent = 100.0 * static_cast<double>(bits_overhead) / static_cast<double>(N_bits);
    return ovhd_percent;
}

static double get_specific_ovhd_in_percent(int64_t specific_bits_used) {
    return 100.0 * specific_bits_used / (double)N_bits;
}

static void print_overhead_line(const string &label, int64_t bits_used) {
    cout << setw(11) << label << setw(13) << right << bits_used << setw(3) << " " << get_specific_ovhd_in_percent(bits_used) << "%" << endl;
}

static void print_separator() {
    cout << "                      ---" << endl;
}


std::string spacebar(int n) {
    std::string s="";
    for(int i=0; i<n; i++) {s+=" ";}
    return s;
}

void print_synth_command() {
    cout << "generating synthetic bits and queries" << endl;
}

void print_filename() {
    std::cout << "reading file: " << in_path << endl << endl;
}


void print_validation_tablenames() {
    cout << endl << "Queries: " << queries.size() << "/" <<  N_queries  << endl;
    cout << "Query" << spacebar(17) << "naive" << spacebar(10) << BVNAME << spacebar(8) << "mismatch" << endl;
}

void print_validation_row(uint64_t code, uint64_t i, uint64_t naive_out, uint64_t complex_out) {
    std::vector<std::string> qNames = {"access ", "rank_0 ", "rank_1 ", "select_0 ", "select_1 "};
    std::string s1 = qNames[code] + std::to_string(i);
    std::string s2 = spacebar(22 - s1.length()) + std::to_string(naive_out);
    std::string s3 = spacebar(14 - std::to_string(naive_out).length());
    if(complex_out!=0xFFFFFFFFFFFFFFFF) {s3 += std::to_string(complex_out);}
    std::string s4 = spacebar(12 - std::to_string(complex_out).length());
    (naive_out == complex_out) || code==0 ? s4 +="." : s4 +="###";
    if (naive_out != complex_out) s4+= "   "+to_string((int64_t)complex_out-(int64_t)naive_out);
    cout << s1 << s2 << s3 << s4 << endl;
}

std::string pad(uint64_t n, int maxlength) {
    if(n==0) return spacebar(maxlength)+"0";
    int len = log10(n);
    return spacebar(maxlength - len) + std::to_string(n);
}


double round_x(double x) {
    if(x < 10) return ((int)(x * 10))/10;
    else return (int)(x);
}

vector<uint64_t> get_query_counts() {
    query_type_counter = std::vector<uint64_t>(5,0);
    for(uint64_t query : queries) {
        uint64_t code = query >> 61;
        query_type_counter[code]++;
    }
    return query_type_counter;
}

void print_query_types(int padlength) {
    auto query_type_counter = get_query_counts();
    //int qlength = log10(N_queries);
    int qlength = padlength;
    //cout << N_queries << " Queries" << endl;
    cout << pad(query_type_counter[0], qlength) << " Accesses" << endl;
    cout << pad(query_type_counter[1], qlength) << " Rank 0s" << endl;
    cout << pad(query_type_counter[2], qlength) << " Rank 1s" << endl;
    cout << pad(query_type_counter[3], qlength) << " Select 0s" << endl;
    cout << pad(query_type_counter[4], qlength) << " Select 1s" << endl;
}

void print_received_input() {
    if(!doMultibench) {
        file_io::printInitialBits(shown_bit_preview_rows);
        cout << "---" << endl;
    }
    int bitlength = log10(N_bits);
    cout << N_bits << " (" << round_x(N_bits/1e9) << "*10^9) Bits" << endl;
    cout << pad(N_zeros,bitlength) << " (" << (((double)N_zeros / N_bits)*100) << "%) 0s " << endl;
    cout << pad(N_ones,bitlength)  << " (" << (((double)N_ones / N_bits)*100) << "%) 1s " << endl;
    print_query_types(bitlength);
}

string get_cpu_info() {
    vector<uint32_t> regs(4);
    string cpu_name;
    for (int i = 0; i <= 2; ++i) {
        __get_cpuid(0x80000002 + i, &regs[0], &regs[1], &regs[2], &regs[3]);
        for (int j = 0; j < 4; ++j) {
            cpu_name += std::string(reinterpret_cast<char*>(&regs[j]), 4);
        }
    }
    // Trim leading and trailing spaces
    cpu_name.erase(std::find(cpu_name.begin(), cpu_name.end(), '\0'), cpu_name.end());
    cpu_name.erase(0, cpu_name.find_first_not_of(" \t\n\r\f\v"));
    cpu_name.erase(cpu_name.find_last_not_of(" \t\n\r\f\v") + 1);
    return cpu_name;
}

void print_cpu_info() {
    cout << get_cpu_info() << endl;
}

string get_compiler_info() {
    #if defined(__clang__)
        return "Clang compiler: " +(string)__clang_version__;
    #elif defined(__GNUC__) || defined(__GNUG__)
        return "GCC compiler: " + (string) __VERSION__;
    #elif defined(_MSC_VER)
        return "MSVC compiler: " +(string) _MSC_FULL_VER;
    #else
        return "Unknown compiler";
    #endif
}

void print_compiler_info() {
    cout << get_compiler_info() << endl;
}



void print_flag_info() {
    //cout << "---------------------------------" << endl;
    print_cpu_info();
    print_compiler_info();
    string const benching_axes[3] = {"sampledist","ratio","nbits"};
        cout << "            -O " << OPTIMIZE << endl;
        cout << "      no debug " << NDEBUG << endl;
        cout << "        x-axis " << benching_axes[BENCHING_AXIS] << endl;
        cout << "   # iter/inst " << iterations << endl;
        cout << "  # inst/param " << instances << endl;
        cout << "             - " << endl;
        cout << "        bitgen " << BITGEN << endl;
        cout << "          bits " << (N_bits/1000000000) << "e9" << endl;
        if(invert_all_bits) {
        cout << "       10ratio " << synth_01ratio << endl;
        } else {
        cout << "       01ratio " << synth_01ratio << endl;
        }
        cout << "             - " << endl;
        if (do_bimodal) {

        cout << "bimodal_factor " << bimodal_factor << endl;
        }

    if(BVNAME=="minflag_ssimd" or BVNAME=="mincount_simd" or BVNAME=="m2" or BVNAME=="m5" or BVNAME=="m7" or BVNAME=="m3") {
        cout << "           pop " << POP << endl;
        cout << "       reverse " << REVERSE << endl;
        cout << "   l1 miniskip " << L1MINISKIP << endl;
        #if PREFETCH_NT
        cout << "   prefetch nt " << PREFETCH_NT << endl;
        #endif
        #if PREFETCH_T0
        cout << "   prefetch t0 " << PREFETCH_T0 << endl;
        #endif
        cout << "     hugepages " << HUGEPAGES << endl;
        cout << "         align " << ALIGN << endl;
        cout << "     alignment " << ALIGNMENT << endl;
        #if CUSTOM_L0>0
        cout << "     custom L0 " << CUSTOM_L0 << endl;
        #endif
        if(synth_seqsize>0) {
            cout << "       seqsize " << synth_seqsize << endl;
        }
        if (BVNAME!="m3") {
        cout << "             - " << endl;
        cout << "  sampling_exp " << SAMPLEDIST_EXP << "  (" << SAMPLEDIST << ")" <<endl;
        cout << "          scan " << SCAN << endl;
        cout << "       extract " << EXTRACT << endl;
        cout << "       reverse " << REVERSE << endl;
        cout << "           pop " << POP << endl;
        }
        if(BVNAME=="m2") {
        cout << "  incl_last_l2 " << INCL_LAST_L2 << endl;
        cout << "     s0_invpop " << INV_POP_SELECT0 << endl;
        }
        if(BVNAME=="m7" || BVNAME=="m3") {
        cout << "    max_sparse " << m3::MAX_SPARSE_PUBLIC << endl;
        }
        cout << "     tree_type " << TREE_TYPE<< endl;
        cout << "tree3_p_choice " << TREE3_STRATEGY<< endl;
        #if TREE_TYPE==5
        cout << "             - " << endl;
        cout << "      tree_exp " << TREE_EXP_LVL0 << endl;
        cout << "      tree_div " << TREE_DIVISION_EXP << endl;
        cout << "max_superblcks " << TREE_MAX_SUPERBLOCKS << endl;
        cout << " inv_expl_ovhd " << TREE_INV_EXPL_OVHD << endl;
        cout << "     pack_expl " << PACK_EXPL << endl;
        #endif
        cout << "             - " << endl;
    }
    if(BVNAME=="sdsl_rrr") {
        cout << " rrr_blocksize " << RRR_BLOCKSIZE << endl;
    }
    // if(BVNAME=="simple_s") {
        // cout << "simple_s_param " << SIMPLE_SELECT_PARAM << endl;
    // }

    else std::cout << BVNAME << endl;
    cout << "---------------------------------" << endl;
}


/*  Validate results against naive implementation */
void validation() {
    useSynth ? synth::generate_synthetic_input() : file_io::readFile(in_path);
    useSynth ? print_synth_command() : print_filename();

    print_received_input();
    print_flag_info();

    auto bv = createBitvector(BVNAME);
    auto naiveBv = createBitvector("naive");
    bv->build_auxiliaries();
    used_space = bv->space_in_bits();
    double ovhd_percent = get_total_overhead_in_percent(used_space);
    cout << "Overhead: " << ovhd_percent << "%" << endl;

    uint64_t constexpr indexMask = (1ULL << 61)-1;
    print_validation_tablenames();
    uint64_t qcount = 0;
    if(doValidw) cout << "validw: only printing wrong queries" << endl << "progress:"<<endl;;

    for(uint64_t query : queries) {
        uint64_t const code = query >> 61;
        uint64_t const i = query & indexMask;
        uint64_t naive_out = 0;
        uint64_t complex_out = 0;
        if(code==0)     {naive_out = naiveBv->access(i);   complex_out = bv->access(i); }
        else if(code==1){naive_out = naiveBv->rank_0(i);   complex_out = bv->rank_0(i);}
        else if(code==2){naive_out = naiveBv->rank_1(i);   complex_out = bv->rank_1(i);}
        else if(code==3){naive_out = naiveBv->select_0(i); complex_out = bv->select_0(i);}
        else if(code==4){naive_out = naiveBv->select_1(i); complex_out = bv->select_1(i);}
        if (!doValidw || naive_out != complex_out) {
            print_validation_row(code, i, naive_out, complex_out);
        }
        qcount++;
        if(qcount % validw_updaterate == 0) {
            cout << qcount << "  " << 100*qcount/(double)queries.size() << "%" << endl;
        }
    }
}



void print_timing_headers() {
    cout << endl;
    cout << "bits";
    cout << (useSynth ? "   Generate" : "   Read    ");
    cout << "      Copy       Build     Solve      query     overhead      01ratio    density  seed" << endl;
}

void print_timings_row() {
    cout << (N_bits/1000000000) << "e9";
    cout << setw(9) << round_x(dRead) << " s";
    cout << setw(9) << round_x(dCopy) << " ms";
    cout << setw(8) << round_x(dBuild) << " ms";
    cout << setw(9) << round_x(dQueries) << " ms";
    cout << setw(6) << round_x(dPerQuery) << " ns";
    cout << setw(11) << round_x(current_overhead_percent*10000)/10000 << "%";
    cout << setw(13) << eff_01ratio << "";
    cout << setw(13) << synth::bit_density << " ";
    cout << setw(10) << seed;
    cout << setw(10) << BVNAME;
    if(BENCHING_AXIS==0) {
        cout << "   SD=" << SAMPLEDIST_EXP;
    }
    cout << endl;
}

void print_summary() {
    cout << "---------------------------------" << endl;
    cout << spacebar(12) << averaged_build << " ms";
    cout << spacebar(5) << averaged_iteration << " ms";
    cout << spacebar(5) << averaged_avg_query << " ns";
    cout << "  " << BVNAME << " algo" << endl;

    cout << "finished benchmark" << endl;
}



std::string get_date(bool prefix=true) {
    auto t = std::chrono::system_clock::now();
    std::time_t date = std::chrono::system_clock::to_time_t(t);
    if (prefix) {
        return "DATE="+(std::string)std::ctime(&date);
    }
    else {
        return std::ctime(&date);
    }
}

void fill_spaces(string &s) {
    std::replace(s.begin(), s.end(), ' ', '_');
}

void strip_linebreaks(string &s) {
    s.erase(std::remove_if(s.begin(), s.end(),
    [](char c) { return c == '\n' || c == '\r'; }),
    s.end());
}

template<typename T>
void res(const std::string& key, T value) {
    std::ostringstream oss;
    oss << value;
    std::string new_val = oss.str();
    fill_spaces(new_val);
    strip_linebreaks(new_val);
    // Try to find and update the key
    for (auto& [k, v] : CURR_RESULT_LINE) {
        if (k == key) {
            v = new_val;
            return;
        }
    }
    // If not found, insert new
    CURR_RESULT_LINE.emplace_back(key, new_val);
}


#define mres(flag) res(#flag, TOSTRING(flag));


void update_RESULT_single_benchmark(bool shadow=false) {
    res("DATE",get_date(false));
    res("BVNAME",BVNAME);
    res("BITGEN", BITGEN);
    res("synth_bits", synth_bits);
    res("curr_instance", curr_instance);
    mres(CUSTOM_L0)
    res("ALPHA",ALPHA);
    res("tree3-strat", TREE3_STRATEGY);
    res("generate(s)", dRead);
    res("build(ms)", averaged_build);
    res("t(ns)", averaged_avg_query);
    res("overhead_percent", current_overhead_percent);

    res("synth_01ratio", synth_01ratio);
    res("density", synth::bit_density);
    res("sinperiods", synth::get_sinus_periods_as_string());
    res("sinthresh", synth::sin_threshhold);
    res("total_used_space_in_bits", total_used_space_in_bits);
    res("eff_01ratio", eff_01ratio);
    res("summary-levels",SUMMARY_LEVELS);
    res("instances", instances);
    res("BV_COMPRESSION",BV_COMPRESSION);

    res("z-search-strat", Z_SEARCH_STRAT);
    // Tree/structure stats
    auto ts = treestats;
    // mres(TREE3_PARAM_CHOICE)
    res("theo_a_exp", ts.theo_a_exp);
    res("theo_b_exp", ts.theo_b_exp);
    res("a_exp", ts.a_exp);
    res("b_exp", ts.b_exp);
    res("topgroups", ts.topgroups);
    res("midgroups", ts.midgroups);
    res("botgroups", ts.botgroups);
    res("parameter_cycles", parameter_cycles);
    res("curr_parameter_cycle", curr_parameter_cycle);
    res("theo_worstcase_tree_size_bits", ts.theo_worstcase_size_bits);
    res("theo_actual_tree_size_bits", ts.theo_actual_size_bits);
    res("final_tree_size_bits", ts.final_tree_size_bits);
    res("theo_worstcase_tree_size_factor", ts.theo_worstcase_size_factor);
    res("theo_actual_tree_size_factor", ts.theo_actual_size_factor);
    res("final_tree_size_factor", ts.final_tree_size_factor);

    res("N_zeros", N_zeros);
    res("N_ones", N_ones);
    res("N_bits", N_bits);
    res("BIMODAL_FACTOR",bimodal_factor);
    res("MAX_SPARSE", m3::MAX_SPARSE_PUBLIC);
    res("seed", seed);
    res("invert_all_bits", invert_all_bits);
    res("bench_iterations", iterations);
    res("shadow",(int)shadow);

    if (do_gap) {
        res("GAP_CRITICAL_QUERIES", GAP_CRITICAL_QUERIES);
        res("GAP_DUMMY_QUERIES", GAP_DUMMY_QUERIES);
        res("gap_size_bits",curr_gap_size_bits);
        res("gap_avg_dummy", gap_avg_dummy);
        res("gap_critical_1st", gap_critical_1st);
        res("gap_critical_2nd", gap_critical_2nd);
        res("gap_critical_3rd", gap_critical_3rd);
        res("gap_critical_avg_cached", gap_critical_avg_cached);
        res("gap_single_dense_1st", gap_single_dense_1st);
        res("gap_single_dense_2nd", gap_single_dense_2nd);
        res("gap_single_dense_3rd", gap_single_dense_3rd);
        res("gap_single_dense_avg_cached", gap_single_dense_avg_cached);
    }
    // Hit stats
    // res("hit_dense", hit.dense_block);
    // res("hit_alldense", hit.alldense_block);
    // res("hit_sparse", hit.sparse_block);
    // res("hit_allsparse", hit.allsparse_block);
    // res("hit_sparse_match", hit.sparse_match);
    // res("hit_sparse_nomatch", hit.sparse_nomatch);


    // Query counts
    auto qc = get_query_counts();
    res("access", qc[0]);
    res("rank0", qc[1]);
    res("rank1", qc[2]);
    res("select0", qc[3]);
    res("select1", qc[4]);

    // Cpu cycle stats
    // res("cycles_in_tree", cycles.in_tree);
    // res("cycles_scan_to_L1", cycles.scan_to_L1);
    // res("cycles_in_SIMD", cycles.in_SIMD);
    // res("cycles_popcounting", cycles.popcounting);
    // res("cycles_zreads", cycles.zreads);
}

void update_RESULT_back_cmake_flags() {

    mres(TREE_TYPE)

    mres(COUNT_HIT_TYPES)
    // mres(REFINED_SIMD)
    // mres(CUSTOM_MAX_SPARSE)
    mres(CHECK_ALL_DENSE_SPARSE)

    mres(BENCHING_AXIS)
    mres(FORCE_CRIT_PATH)
    mres(WARMUP_ROUNDS)

    mres(SCAN)
    mres(EXTRACT)
    mres(REVERSE)
    mres(POP)
    // mres(L1MINISKIP)
    mres(PREFETCH_NT)
    mres(PREFETCH_T0)
    mres(HUGEPAGES)
    mres(ALIGN)
    // mres(AVOID_PDEP)

    // mres(PACK_TREE)
    // mres(TREE_EXP_LVL0)
    // mres(TREE_DIVISION_EXP)
    // mres(TREE_MAX_SUPERBLOCKS)
    // mres(TREE_INV_EXPL_OVHD)

    mres(RRR_BLOCKSIZE)
    // mres(SIMPLE_SELECT_PARAM)

    res("-march","native");
    #if NDEBUG
    res("-DNDEBUG","1");
    #elif
    res("-DNDEBUG","0");
    #endif
    mres(OPTIMIZE)

    res("CPU",get_cpu_info());
    string compiler;
    #if defined(__clang__)
    compiler="Clang "+(string)__clang_version__;
    #elif defined(__GNUC__) || defined(__GNUG__)
    compiler="GCC " + (string) __VERSION__;
    #else
    cout << "Compiler is neither gcc nor clang" << endl;
    compiler="not_recognized"
    #endif
    res("Compiler",compiler);

}


bool is_curr_result_shadowed() {
    for (auto& [k, v] : CURR_RESULT_LINE) {
        if (k=="shadow") {
            return v=="1";
        }
    }
    return false;
}


void write_RESULT_line() {
    // update_RESULT_header_flags();
    update_RESULT_single_benchmark();
    update_RESULT_back_cmake_flags();
    std::ofstream file(mb_path,std::ios::app);
    if (!is_curr_result_shadowed()) {
        file << "RESULT";
    }
    for (auto& [k, v] : CURR_RESULT_LINE) {
        file << " " << k <<"=" << v;
    }
    file << std::endl;
    file.close();

}

void write_custom_line(const string s) {
    std::ofstream file(mb_path,std::ios::app);
    file << s << std::endl;
    file.close();
}


std::string get_unixtimestamp() {
    const auto p1 = std::chrono::system_clock::now();
    return to_string(std::chrono::duration_cast<std::chrono::seconds>(p1.time_since_epoch()).count());
}


bool mb_file_exists() {
    ifstream f(mb_path.c_str());
    return f.good();
}

void mb_set_path() {
    string mb_file_name = (mb_name.empty() ? get_unixtimestamp() : mb_name);
    mb_path = mb_write_folder + mb_file_name + ".txt";
}



string mb_get_cpu_flags() {
    string details;
    #if defined(__clang__)
    details+=">COMP_VERB=Clang compiler: " +((string)__clang_version__)+"\n";
    details+=">COMP_SHORT=clang\n";
    #elif defined(__GNUC__) || defined(__GNUG__)
    details+=">COMP_VERB=GCC compiler: " + (string) __VERSION__+"\n";
    details+=">COMP_SHORT=gcc\n";
    #else
    cout << "Compiler is neither gcc nor clang" << endl;
    details+="NO RECOGNIZED COMPILER\n"
    #endif
    details+=">Opt=-O"+string(TOSTRING(OPTIMIZE));
    return details;
}



double get_duration(tp tStart, tp tEnd) {
    double dur = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tStart).count()/1000.0; //ms
    dur = ((int)(dur * 10.0))/10.0; //round to one digit after decimal point
    return dur;
}


//After a benchmark, run some sanity check queries, maybe uncovers some rare edgecases
void mini_postrun_validation(std::unique_ptr<Bitvector> &bv) {
    // printf("mini validation\n");
    auto naiveBv = createBitvector("naive");
    constexpr uint64_t random_queries = 50;
    // constexpr uint64_t sequential_queries = 20000;
    uint64_t constexpr indexMask = (1ULL << 61)-1;
    for (uint64_t q=0; q < min(random_queries, queries.size()); q++) {
        uint64_t query = queries[q];
        uint64_t const code = query >> 61;
        uint64_t i = query & indexMask;
        // if (q>random_queries) {
        //     uint64_t inject = (q+1)*63;
        //     switch (code) {
        //         case 0: i = min(N_bits-1, inject); break;
        //         case 1: i = min(N_bits-1, inject); break;
        //         case 2: i = min(N_bits-1, inject); break;
        //         case 3: i = min(N_zeros, inject); break;
        //         case 4: i = min(N_ones, inject); break;
        //         default: ;
        //     }
        // }
        // else {
            // cout << q << endl;
        // }
        uint64_t naive_out = 0;
        uint64_t complex_out = 0;
        if(code==0)     {naive_out = naiveBv->access(i);   complex_out = bv->access(i); }
        else if(code==1){naive_out = naiveBv->rank_0(i);   complex_out = bv->rank_0(i);}
        else if(code==2){naive_out = naiveBv->rank_1(i);   complex_out = bv->rank_1(i);}
        else if(code==3){naive_out = naiveBv->select_0(i); complex_out = bv->select_0(i);}
        else if(code==4){naive_out = naiveBv->select_1(i); complex_out = bv->select_1(i);}
        if (naive_out != complex_out) {
            print_validation_row(code, i, naive_out, complex_out);
            string s="Error! query k="+to_string(i)+" has mismatch. naive="+to_string(naive_out)+",  bv="+to_string(complex_out);
            write_custom_line(s);
        }
    }
}




uint64_t global_xor_dummy = 0x12345678;
/* Solve all queries */
std::vector<uint64_t> solve_queries() {
    auto bv = createBitvector(BVNAME);
    hit = {};//reset
    cycles = {};

    //Within build_auxiliaries() an additional time point "tCopyFinished" is set
    //to distinguish between stupid copying and interesting building
    tp tAux_Begin = std::chrono::steady_clock::now();
    bv->build_auxiliaries();
    tp tAux_End = std::chrono::steady_clock::now();
    dCopy  = get_duration(tAux_Begin, tCopyFinished);
    dBuild = get_duration(tCopyFinished, tAux_End);

    used_space = bv->space_in_bits();
    current_overhead_percent = get_total_overhead_in_percent(used_space);
    total_used_space_in_bits = used_space;

    std::vector<uint64_t> results(0);
    //skipping the result storing, if we are interested in validation we run that separately, see function
    //Storing also not needed to enforce calculation, this is already enforced via XORing

    uint64_t q_counter=0;
    uint64_t constexpr indexMask = (1ULL << 61)-1;
    uint64_t result = 0;

    tp t0 = std::chrono::steady_clock::now();
    for(uint64_t const query : queries) {
        uint64_t const code = query >> 61;
        uint64_t const i = query & indexMask;
        if     (code==0) {result = bv->access(i);}
        else if(code==1) {result = bv->rank_0(i);}
        else if(code==2) {result = bv->rank_1(i);}
        else if(code==3) {result = bv->select_0(i);}
        else if(code==4) {result = bv->select_1(i);}
        //results[q_counter] = result;
        ++q_counter;
        #if FORCE_CRIT_PATH
        global_xor_dummy ^= result;
        #endif
    }
    tp t1 = std::chrono::steady_clock::now();
    dQueries = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()/1000.0; //ms
    dPerQuery = 1e6 * dQueries/N_queries;//ns
    mini_postrun_validation(bv);
    return results;
}


void gap_many_dummy_queries(uint64_t q_counter, Bitvector *bv) {
}

void gap_single_critical_query(Bitvector *bv) {

}

std::vector<uint64_t> solve_gap_gadget() {
    auto bv = createBitvector(BVNAME);

    tp tAux_Begin = std::chrono::steady_clock::now();
    bv->build_auxiliaries();
    tp tAux_End = std::chrono::steady_clock::now();
    dCopy  = get_duration(tAux_Begin, tCopyFinished);
    dBuild = get_duration(tCopyFinished, tAux_End);

    used_space = bv->space_in_bits();
    current_overhead_percent = get_total_overhead_in_percent(used_space);
    total_used_space_in_bits = used_space;

    uint64_t q_counter=0;
    uint64_t constexpr indexMask = (1ULL << 61)-1;
    uint64_t result = 0;

    //To know the number of the first 1-bit after the gap, we need to know those before the gap
    //Not all bitvectors support rank, so for uniformity we just do it manually here
    uint64_t ones_before_gap = 0;
    for (uint64_t w=0; w<N_words/2; w++) {
       ones_before_gap += popcount(bits[w]);
    }
    uint64_t first_m_after_gap = ones_before_gap + 1;


    cout << "ones_before_gap = " << ones_before_gap << endl;

    //quickly validate that we have hit the gap
    uint64_t select_before = bv->select_1(first_m_after_gap-1);
    uint64_t select_after = bv->select_1(first_m_after_gap);
    uint64_t diff = select_after - select_before;
    cout << "select on last before gap = " << select_before << endl;
    cout << "expected gap center       = " << N_bits/2 << endl;
    cout << "select on first after gap = " << select_after << endl;
    cout << "expected gap size         = " << curr_gap_size_bits << endl;
    cout << "measured gap size         = " << diff << endl;


    // std::vector<double> critical_query_time_list = std::vector<double>(GAP_QUERIES);
    // std::vector<double> avg_dummy_query_time_list = std::vector<double>(GAP_QUERIES);
    // std::vector<double> single_dummy_time_list = std::vector<double>(GAP_QUERIES);
    highres_tp t0;
    highres_tp t1;

    for (int gap_iter = 0; gap_iter < GAP_CRITICAL_QUERIES; gap_iter++) {
        //Many random dummy queries into the left dense part, to empty the cache (and get some baseline times)
        t0 = std::chrono::high_resolution_clock::now();
        for (int dummy=0; dummy<GAP_DUMMY_QUERIES; dummy++) {
            q_counter++;
            q_counter = q_counter % generated_dummy_queries; //we will loop multiple times through the pre-generated query list
            auto Q = queries[q_counter];
            uint64_t const i = Q & indexMask;
            uint64_t result = bv->select_1(i);
            global_xor_dummy ^= result;
        }
        t1 = std::chrono::high_resolution_clock::now();
        gap_avg_dummy = std::chrono::duration_cast<std::chrono::nanoseconds>(t1- t0).count()/GAP_DUMMY_QUERIES;


        //Now one critical query right after the gap
        t0 = std::chrono::high_resolution_clock::now();
        uint64_t index_after_gap = first_m_after_gap + (global_xor_dummy & 0xF); //use the xor_dummy to quickly give us a pseudo-random offset between 0 and 15
        global_xor_dummy ^= bv->select_1(index_after_gap);
        t1 = std::chrono::high_resolution_clock::now();
        gap_critical_1st = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();


        //Repeat this critical query once. Now it is cached, so it should be much faster
        //We only repeat it once for now, to also measure how much influence, if any, the fact has that we only measure a single query
        t0 = std::chrono::high_resolution_clock::now();
        global_xor_dummy ^= bv->select_1(index_after_gap);
        t1 = std::chrono::high_resolution_clock::now();
        gap_critical_2nd = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

        t0 = std::chrono::high_resolution_clock::now();
        global_xor_dummy ^= bv->select_1(index_after_gap);
        t1 = std::chrono::high_resolution_clock::now();
        gap_critical_3rd = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

        //Repeat the critical query many times. This should ideally give on average the same time as the query above.
        //If not, it gives us an indication how inprecise our timings are for single queries
        constexpr int CACHED_REPS = 10000;
        t0 = std::chrono::high_resolution_clock::now();
        for (int cached_rep=0; cached_rep<CACHED_REPS; cached_rep++) {
            global_xor_dummy ^= bv->select_1(index_after_gap);
        }
        t1 = std::chrono::high_resolution_clock::now();
        gap_critical_avg_cached = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()/CACHED_REPS;


        //As a next baseline, also query into the right dense part. This should not be cached, so rather slow, but it is still an easy query structure-wise, so this gives us more information about cache vs. structure influences
        //The right dense part starts at m/2, we can safely go somewhere in the middle of that, i.e. 3/4m, plus some small offset, to randomize where exactly we end up
        t0 = std::chrono::high_resolution_clock::now();
        uint64_t right_dense_index = (3*N_ones)/4 + (global_xor_dummy & 0x7FFFF); //again mis-use xor_dummy as a random-number-generator, now for 0...500k
        global_xor_dummy ^= bv->select_1(right_dense_index);
        t1 = std::chrono::high_resolution_clock::now();
        gap_single_dense_1st = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();


        //Repeat once, now cached
        t0 = std::chrono::high_resolution_clock::now();
        global_xor_dummy ^= bv->select_1(right_dense_index);
        t1 = std::chrono::high_resolution_clock::now();
        gap_single_dense_2nd = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

        t0 = std::chrono::high_resolution_clock::now();
        global_xor_dummy ^= bv->select_1(right_dense_index);
        t1 = std::chrono::high_resolution_clock::now();
        gap_single_dense_3rd = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

        //Repeat many times, still cached, hopefully with similar measurements as the single-cached result, otherwise reveals limits in measurement-technique
        t0 = std::chrono::high_resolution_clock::now();
        for (int cached_rep=0; cached_rep<CACHED_REPS; cached_rep++) {
            global_xor_dummy ^= bv->select_1(right_dense_index);
        }
        t1 = std::chrono::high_resolution_clock::now();
        gap_single_dense_avg_cached = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()/CACHED_REPS;

        write_RESULT_line();

    }
    return {};
}


//Answer a given instance of bitvector+queries
void bench(bool silent=false) {
    if(!silent) (useSynth ? print_synth_command() : print_filename());
    double total_query_time = 0;
    double total_build_time = 0;
    std::vector<uint64_t> results{};
    int iteration = 0;
    if(!silent) { print_received_input();}
    if(!even_more_silent) {print_flag_info();}
    while(iteration < iterations) {
        if (do_gap) {
            results = std::move(solve_gap_gadget()); //specialized benchmarking on the gap gadget
        } else {
            results = std::move(solve_queries()); //normal benchmarking
        }
        print_timings_row();
        total_query_time += dQueries;
        total_build_time += dBuild;
        iteration++;
        //if(doMultibench) mb_write_single_tQueries();
    }
    if (write_individual_query_results && !useSynth) {
        file_io::write_to_file(results, result_path);
    }
    averaged_avg_query = round_x((total_query_time /(iterations * N_queries))*1e9)/1000.0; //result is in ns
    averaged_build = round_x(total_build_time / iterations);
    averaged_iteration = round_x(total_query_time/iterations);
    if(!silent) print_summary();
}

//Main benchmarking loop
void multibench() {
    mb_set_path();
    cout << "writing to: " << mb_path << endl;
    print_timing_headers();
    for(auto n: N_BITS) {
        synth_bits = n;
        seed = 90;
        for(int i=0; i<instances; i++) {
            curr_instance = i;
            seed += 10; //New seed for random instance generation
            for(auto ratio01 : synth_01ratios) {
                synth_01ratio = ratio01;
                //Now we know everything there is to know about the instance and can generate it
                //This input generation is the most time-consuming aspect of the whole benchmark, so we lifted it as far up as possible
                bool silent = (i>0 && instances>1);
                #if SILENT
                silent=true;
                #endif
                if (BITGEN!="sin") { //Mix in new seeds also at every ratio step, except for sinusoidal where we want to keep the same periods over one full ratio-scan
                    seed += 1;
                }
                for (uint64_t gap : GAP_SIZES) {
                    curr_gap_size_bits = gap; //needs to be set before the instance is generated, in case of gap_gadget instance
                    useSynth ? synth::generate_synthetic_input() : file_io::readFile(in_path);
                    //Code here is really not elegant. Should refactor to a cleaner combinatorial expression
                    //Cycle through competitors
                    for(string const bv : BVNAMES_LIST) {
                        BVNAME = bv;
                        if ((BVNAME=="poppy" || BVNAME=="poppy2") && n >= POPPY_MAX_N) {
                            continue; //poppy only supports instances up to n=2^32 bits
                        }
                        //if m3: try different strategies how to select (a,b) parameters for the sample tree
                        for (int bv_cp : BV_COMPRESSION_LIST) {
                            BV_COMPRESSION = bv_cp;
                            for (int strat : TREE3_STRATEGY_LIST) {
                                TREE3_STRATEGY = strat;
                                //if m3: try different alpha values for the sample tree
                                for (int alpha : ALPHA_LIST) {
                                    ALPHA = alpha;
                                    if (TREE3_STRATEGY==GET_SMALLEST_TREE_PARAMS && ALPHA<8) {
                                        continue;
                                    }
                                    //if m3: try different numbers of summary levels
                                    for (int summlevels : SUMMARY_LEVEL_LIST) {
                                        SUMMARY_LEVELS=summlevels;
                                        //if m3: try either linear or binary search in z-values
                                        for (int z_search_strat : Z_SEARCH_STRAT_LIST) {
                                            Z_SEARCH_STRAT=z_search_strat;
                                            //if m3: cycle through many different sample tree (a,b) choices, benchmark each
                                            for (int c=0; c<parameter_cycles; c++) {
                                                curr_parameter_cycle = c;
                                                //if doing gap gadget: cycle through different gap sizes
                                                curr_gap_size_bits = gap;
                                                bench(silent);
                                                write_RESULT_line();
                                                if (TREE3_STRATEGY!=GET_FASTEST_TREE_PARAMS || BVNAME!="m3") { //parameter cycles only relevant for m3 when searching for the fastest tree-params
                                                    break;
                                                }
                                            }
                                            if (BVNAME!="m3") {break;} //parameter scans only relevant for m3, skip for competition
                                        }
                                        if (BVNAME!="m3") {break;}
                                    }
                                    if (BVNAME!="m3") {break;}
                                }
                                if (BVNAME!="m3") {break;}
                            }
                            if (BVNAME!="m3") {break;}
                        }
                    }
                    if (!do_gap) {
                        break;
                    }
                }

            }
        }
    }
    cout << "written to: " << mb_path << endl;
    #if FORCE_CRIT_PATH
    //Printing dummy to enforce that it not optimized away
    cout << "global XOR-dummy " << global_xor_dummy << endl;
    #endif
}



bool isCorrectBvName() {
    if (any_of(names.begin(), names.end(), [&](const std::string& elem) { return elem == BVNAME; })) {
        return true;
    } else {
        cout << endl << "ERROR" << endl << "Incorrect bvname spelling '" << BVNAME << "'" << endl;
        cout << "Must be one of these:" << endl;
        for(auto name : names) {
            cout << "'" << name <<"'"<< endl;
        }
        return false;
    }
}

//Read the n_bits values from command-line and store times *10^9 in N_BITS
void store_nbits_list(const std::string& n_bits_list) {
    N_BITS.clear();
    std::string current;
    for (char ch : n_bits_list) {
        if (ch != ',') {
            current += ch;
        }
        else {
            //
            uint64_t n_bits = stoi(current) * 1e9;
            cout << "read nbits: " << n_bits << endl;
            N_BITS.push_back(n_bits);
            current.clear();
        }
    }
    // Add the last token if any
    if (!current.empty()) {
        uint64_t n_bits = stoi(current) * 1e9;
        cout << "read nbits: " << n_bits << endl;
        N_BITS.push_back(n_bits);
    }
}


template <typename T>
void store_integer_list(std::vector<T> &V, const std::string &printname, const std::string &param_list) {
    std::stringstream ss(param_list);
    std::string token;
    V.clear();
    while (std::getline(ss, token, ',')) {
        cout << printname << " " << token << endl;
        V.push_back(std::stoi(token));
    }
}


void store_bvnames(const std::string& bvnames_list) {
    std::stringstream ss(bvnames_list);
    std::string token;
    BVNAMES_LIST.clear();
    while (std::getline(ss, token, ',')) {
        cout << "BV NAME " << token << endl;
        BVNAMES_LIST.push_back(token);
    }
}



bool is_official_run(int argc, char * argv[]) {
    return (argc==3 && argv[1][0]!='-');
}

int temp_ratio_factor=0;
int temp_ratio_start=0;

void read_flags(int argc, char * argv[]) {
    cout << argc << " argc" << endl;
    if(argc>1) {
        for(int i=1; i<argc; i+=2) {
            const string flag = string(argv[i]);
            const char* param = argv[i+1];
            const string string_param = string(param);
            cout << flag << setw(11) << string_param << endl;
            if(flag=="-bvname") {
                store_bvnames(string_param);
            }
            if(flag=="-job") {
                JOB = string_param;
            }
            if(flag=="-bitgen") {
                BITGEN = string_param;
            }
            if(flag=="-bmfactor") {
                bimodal_factor = stoi(param);
            }
            if(flag=="-qtype") {
                QTYPE = string_param;
            }
            if(flag=="-iterations") {
                iterations = stoi(param);
            }
            if(flag=="-instances") {
                instances = stoi(param);
            }
            if(flag=="-parameter_cycles") {
                parameter_cycles = stoi(param);
            }
            if (flag=="-tree3-strat") {
                store_integer_list(TREE3_STRATEGY_LIST, "tree strategy", string_param);
            }
            if (flag=="-alpha") {
                store_integer_list(ALPHA_LIST, "alpha", string_param);
            }
            if (flag=="-summary-levels") {
                store_integer_list(SUMMARY_LEVEL_LIST, "summary levels", string_param);
            }
            if (flag=="-z-search-strat") {
                store_integer_list(Z_SEARCH_STRAT_LIST, "z search strategy", string_param);
            }
            if (flag=="-bv-compression") {
                store_integer_list(BV_COMPRESSION_LIST, "bv compression", string_param);
            }
            if(flag=="-nbits") {
                store_nbits_list(param);
            }
            if (flag=="-gap-size") {
                store_integer_list(GAP_SIZES, "gap size", string_param);
            }

            if(flag=="-invert_all") {
                invert_all_bits=stoi(param);
            }
            if(flag=="-ranks") {
                synth_ranks = stod(param);
            }
            if(flag=="-selects") {
                synth_selects = stod(param);
            }
            if(flag=="-select0s") {
                synth_select0s = stod(param);
            }
            if(flag=="-select1s") {
                synth_select1s = stod(param);
            }
            if (flag=="-01ratio") {
                store_integer_list(synth_01ratios, "01ratio", string_param);
            }
            if(flag=="-filename_out") {
                mb_name=string_param;
            }
            if(flag=="-folder_out") {
                mb_write_folder=string_param;
            }
            if(flag=="-progress") {
                cout << "Progress: " << string_param << endl;
            }
        }
    }
}

int main(int argc, char * argv[]) {
    // if(is_official_run(argc,argv)) {
        // official_run(argc, argv);
        // return 0;
    // }
    if(run_testsection) {
        //mincount2_tests();
        //ORourke::demo();
        //mincount_simd::tests();
        //playground::submin();
        //playground::identicalL2();
        //playground::kmax();
        //synth::generate_synthetic_input();
        //Tree::tests();
        //playground::delta13();
        //playground::modify_assembly();
        //playground::pack();
        // playground::test_pdep_alternative();

    }
    else {
        cout << "---------------" << endl;
        read_flags(argc, argv);
        init_benchvalues();
        init_selections();
        init_counters();
        init_paths();
        set_sample_exp(SAMPLEDIST_EXPS[0]);
        print_cpu_info();
        print_compiler_info();
        cout << JOB << endl;
        if (not isCorrectBvName()) return 1;
        // if(doBench) bench();
        if(doValid) validation();
        if(doMultibench) multibench();
    }
}


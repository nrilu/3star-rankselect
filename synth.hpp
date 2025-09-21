//
// Created by nicco on 7/31/24.
//

#pragma once

#include "bv.hpp"
#include "alignedallocator.hpp"
// #include "plotter.hpp"

#include <vector>
#include <random>
//#include <bits/ranges_algo.h>
#include <unordered_map>
#include <cmath>
#include <unordered_set>

#include <sys/mman.h> //for huge tables via madvise(...)
//#include <unistd.h> //sysconf()
//#include <bits/confname.h> //_SC_PAGESIZE


//-------------------- GENERATE INPUT ------------------------------------------//
//generate random bitvector during runtime

constexpr bool print_chunk_info = false;

namespace synth {
    using namespace std; //not nice programming style, but we now already commited for it - refactor sometime
    typedef uint64_t u64;

    const std::unordered_map<std::string, uint64_t> map = {
        {"access", 0},
        {"rank 0", 1},
        {"rank 1", 2},
        {"select 0", 3},
        {"select 1", 4},
    };

    //Generate whole random 64-Bit words for fast 1:1 ratio
    inline void fast_1to1_ratio() {
        // std::random_device rd;
        std::mt19937_64 gen(seed);
        std::uniform_int_distribution<uint64_t> rand64Bit(0, UINT64_MAX);

        for (uint64_t w = 0; w < N_words; ++w) {
            uint64_t W = rand64Bit(gen);
            bits[w] = W;
            uint64_t pop = std::popcount(W);
            N_ones += pop;
            N_zeros += 64 - pop;
        }
    }

    //Generate sequences of only 1s or only 0s, each with 64*synth_sequencesize bits
    inline void generate_sequences() {
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<int> randBinary(0, 1);
        u64 w=0;
        while(w < N_words) {
            //Random choice of either 64 Zeros or 64 Ones
            u64 W = 0; //W = 64 Zeros
            int pop = 0;
            if(randBinary(gen) == 1) {
                W = 0xffffffffffffffffULL; //W = 64 Ones
                pop = 64;
            }
            //Fill bitvector with choosen bit-type W
            for(int i=0; i<synth_seqsize && w<N_words; i++) {
                bits[w++] = W;
                N_ones += pop;
                N_zeros += 64-pop;
            }
        }
    }

    // //Generate random bits with bias towards 0s as given by synth_ratio0sTo1s
    // inline void generate_biased_distr() {
    //     std::random_device rd;
    //     std::mt19937_64 gen(rd());
    //     std::uniform_int_distribution<> distr(0, synth_01ratio);
    //
    //     uint64_t W=0; //current word to fill
    //     for (u64 w = 0; w < N_words; w++) {
    //         for(int i=0; i<64; i++) {
    //             if(distr(gen) == 0) {
    //                 W |= 1ULL << i;
    //             }
    //         }
    //         bits[w] = W;
    //         int pop = std::popcount(W);
    //         N_ones += pop;
    //         N_zeros += 64 - pop;
    //         W=0;
    //         }
    // }


    //Combine query parts and save it
    inline void writeQuery(const std::string &name, const uint64_t index) {
        uint64_t query_coded = 0;
        uint64_t code = map.at(name);
        query_coded |= code << 61;
        query_coded |= index;
        queries.push_back(query_coded);
    }



    inline void fillchunk(int dist, uint64_t count, uint64_t& w) {
        if(dist<64) {
            int b = 64/dist;
            uint64_t W = (1ULL << b)-1;
            if(b==64) {
                W = 0xFFFFFFFFFFFFFFFF;
            }
            while(count>0 && w < N_words) {
                bits[w++]=W;
                count -= b;
                N_ones += b;
                N_zeros += 64-b;
            }
        } else {
            constexpr uint64_t W0 = 0;
            constexpr uint64_t W1 = 1;
            while(count>0 && w < N_words) {
                int skips = dist/64-1;
                bits[w++]=W1;
                N_ones+=1;
                N_zeros+=63;
                while(skips>0 && w < N_words) {
                    bits[w++]=W0;
                    N_zeros+=64;
                    skips--;
                }
                count--;
            }
        }
    }

    inline void chunky() {
        std::random_device rd;
        //uint64_t seed = 12345;
        std::mt19937_64 rn(seed);

        int avgdist = synth_01ratio + 1; //converting from ratio to average distance

        int min_count_exp = 6;
        int max_count_exp = 17;
        std::uniform_int_distribution<int> count_exp(min_count_exp, max_count_exp);
        std::uniform_real_distribution<> randreal(0.0,1.0);

        //define deviation from avgdist
        vector<double> mult = {1.0/32, 1.0/8, 1, 8, 64};
        vector<int> dist{};
        if(print_chunk_info) cout << "d=" << endl;
        for(auto m : mult) {
            int d = std::max(1,(int)(avgdist*m));
            dist.push_back(d);
            if(print_chunk_info) cout << "  " <<d << endl;
        }

        //probabilities s.t. average distance == avgdist
        double pair_prob = 0.4;
        double p0 = (dist[4]-(double)avgdist) / (dist[4]-dist[0]);
        double p1 = (dist[3]-(double)avgdist) / (dist[3]-dist[1]);
        double p2 = (1 - 2*pair_prob)/pair_prob;
        double p3 = 1-p1;
        double p4 = 1-p0;

        vector<double> p = {p0,p1,p2,p3,p4};
        if(print_chunk_info) cout << "p= " << endl;
        for(auto &i : p) {
            i*=pair_prob;
            if(print_chunk_info) cout << "  " << i << endl;
        }

        vector<int> chunktype(5,0);
        uint64_t w = 0;
        while(w < N_words) {
            int count = 1 << count_exp(rn);
            double const r = randreal(rn);
            int i=0;
            double psum = 0;
            while(i<4 && r > (psum += p[i])) {
                i++;
            }
            chunktype[i]++;
            fillchunk(dist[i],count,w);
        }
        if(print_chunk_info) {
            cout << "chunktype" << endl;
            for(auto c : chunktype) {
                cout << "   "<< c << endl;
            }
        }
        double avg_1_spacing = N_words*64 / (double)N_ones;
        //cout << "                                               Avg Spacing = " << avg_1_spacing << endl;

    }

    //each bit is flipped with probability 1/(1+synth_01ratio) to a 1, independent from all other bits
    inline void smooth() {
        std::random_device rd;
        std::mt19937_64 rn(seed);
        std::uniform_int_distribution<int> ratio_range(0, synth_01ratio);
        uint64_t w = 0;
        while(w < N_words) {
            uint64_t W = 0;
            for(int i=0; i<64; i++) {
                if(ratio_range(rn) == 0) {
                    W |= (1ULL << i);
                }
            }
            bits[w]=W;
            uint64_t pop = popcount(W);
            N_ones += pop;
            N_zeros += 64-pop;
            w++;
        }

    }

    //A version of smooth that has factor 64 less calls to the random number generator, increasing bit generation speed
    //It samples entire 64-bit words at once, giving only an approximation to a true independent smooth distribution
    inline void quick_smooth() {
        std::random_device rd;
        std::mt19937_64 gen(rd());
        int r = synth_01ratio;
        if(r<=100) {
            //Given some desired avg number n of 1s, sample a 64-bit domain of length 2n bits
            //On avg, half of the 2n bits will be 1s, matching the desired target.
            int n_1s = 64/(1+r);
            int bits_to_sample = 2*n_1s;
            std::uniform_int_distribution<uint64_t> distr(0, (1ULL<<bits_to_sample)-1);
            for(uint64_t w= 0; w<N_words; w++) {
                uint64_t W = distr(gen);
                bits[w] = W;
                int pop = std::popcount(W);
                N_ones += pop;
                N_zeros += 64 -pop;
            }
        }
        else {
            int n_zero_words = lround(r/64);
            std::uniform_int_distribution<uint64_t> distr(0, 3); //on avg. one 1-bit
            for(uint64_t w=0; w<N_words; w++) {
                if(w%(n_zero_words+1)==0) {
                    bits[w] = distr(gen);
                    int pop = std::popcount(bits[w]);
                    N_ones += pop;
                    N_zeros += 64 -pop;
                }
                else {
                    bits[w] = 0;
                    N_zeros +=64;
                }
            }
        }
    }


    inline void fractal() {

        std::random_device rd;
        uint64_t seed_fractal = 12348; //12345
        std::mt19937_64 rn(seed_fractal);
        int max_len = 13;
        int min_spacing = 0;
        int max_spacing = 7;
        double p64 = 0.1;
        uint64_t w = 0;
        std::uniform_int_distribution<int> rand_exp(0, max_len);
        std::uniform_int_distribution<int> rand_spacing(min_spacing, max_spacing);
        std::uniform_int_distribution<int> rand_p64(0, 1/p64);

        uint64_t W0 = 0;
        uint64_t W1 = 1;
        uint64_t W64 = 0xFFFFFFFFFFFFFFFF;
        while(w < N_words) {
            int len = 1 << rand_exp(rn);
            int spacing = 1 << rand_spacing(rn);
            if(rand_p64(rn)==0) {
                for(int i=0; i<20*spacing && w < N_words; i++) {
                    bits[w++] = W64;
                    N_ones += 64;
                }
            }
            for(int v=0; v<len && w < N_words; v++) {
                bits[w++] = W1;
                N_ones += 1;
                N_zeros += 63;
                for(int i=0; i<spacing && w < N_words; i++) {
                    bits[w++] = W0;
                    N_zeros+=64;
                }
            }
        }

        double avg_1_spacing = N_words*64 / N_ones;

        cout << "avg spacing " << avg_1_spacing << endl;
        //plotter::bitscatter();
    }


    //Helperfunction for Craig Gidneys sampling
    inline void or_low_probability_bits(double probability, std::mt19937_64 &rng) {
        if (probability == 0) {
            return;
        }
        double d = log1p(-probability);
        size_t i = 0;
        size_t n = N_bits;
        while (true) {
            double u = ((double)rng() + 0.5) / pow(2, 64);
            i += (uint64_t)floor(log(u) / d);
            if (i >= n) {
                return;
            }
            bits[i / 64] |= 1ULL << (i & 63);
            i++;
        }
    }

    //Using Craig Gidneys fast biased bit generation code
    //https://algassert.com/post/2200
    inline void gidney_generate_biased_random_bits(double probability, std::mt19937_64 &rng) {
        constexpr size_t COIN_FLIPS = 8;
        constexpr double BUCKETS = (double)(1 << COIN_FLIPS);
        double raised = probability * BUCKETS;
        double raised_floor = floor(raised);
        double raised_leftover = raised - raised_floor;
        double p_truncated = raised_floor / BUCKETS;
        double p_leftover = raised_leftover / BUCKETS;
        uint64_t p_top_bits = (uint64_t)raised_floor;

        // Flip 8 coins for each output bit, using the position of the first HEADS
        // result to select a bit from the probability's binary representation.
        for(uint64_t w = 0; w< N_words; w++) {
            uint64_t alive = rng();
            uint64_t result = 0;
            for (size_t k_bit = COIN_FLIPS - 1; k_bit--;) {
                uint64_t shoot = rng();
                result ^= shoot & alive & -((p_top_bits >> k_bit) & 1);
                alive &= ~shoot;
            }
            bits[w] = result;
        }

        // De-truncate the probability by refining the output.
        or_low_probability_bits(p_leftover / (1 - p_truncated), rng);
    }

    //Using Craig Gidneys fast biased bit generation code
    //https://algassert.com/post/2200
    inline void start_gidney_fast_independent_bits() {
        //For 1:1 ratio we dont need fancy techniques can just sample whole random 64-bit words
        if(synth_01ratio==1) {
            fast_1to1_ratio();
        }
        else {
            //For probabilities <=0.5 use Gidney
            double probability = 1.0/(1+synth_01ratio);
            std::mt19937_64 rng(seed);
            gidney_generate_biased_random_bits(probability, rng);
            for(auto W : bits) {
                int pop = std::popcount(W);
                N_ones += pop;
                N_zeros += 64-pop;
            }
        }

    }




    inline void manual_queries() {

        // writeQuery("select 1",281929);
        // writeQuery("select 1",281930);
        // writeQuery("select 1",281931);
        // writeQuery("select 1",2043);
        // writeQuery("select 1",47893);
        // writeQuery("select 1",47894);
        // writeQuery("select 1",65535);
        // writeQuery("select 1",65536);
        // writeQuery("select 1",65537);
        // writeQuery("select 1",262143);
        // writeQuery("select 1",262144);
        // writeQuery("select 1",262145);
        // writeQuery("select 1",4095);
        // writeQuery("select 1",4096);
        // writeQuery("select 1",4097);
        // writeQuery("select 1",7);
        // writeQuery("select 1",2113);
        // writeQuery("select 1",2114);

        for(int i=0; i<2000000; i+=1) {
            writeQuery("rank 1",i);
        }


        // writeQuery("rank 1",149624);
        // writeQuery("rank 1",149625);

        // writeQuery("rank 1",16384);
        // writeQuery("rank 1",16385);
        // writeQuery("rank 1",16386);

        // writeQuery("rank 1",0);
        // writeQuery("rank 1",488);
        // writeQuery("rank 1",489);

        // writeQuery("rank 1",1);
        // writeQuery("rank 1",2);

        // writeQuery("rank 1",1024);
        // writeQuery("rank 1",1025);

        // for(int i=1; i<2000000; i+=1) {
        // writeQuery("select 1",i);
        // }

        // int start = 4750000;
        // for(int i=start; i<start+100000; i+=1) {
        // writeQuery("select 1",i);
        // }


        // writeQuery("select 1",999122574);

        // writeQuery("select 1",544216);

        // writeQuery("select 1",97395);

        // writeQuery("select 1",1976730);
        // writeQuery("select 1",1976731);

        // writeQuery("select 1",150277600);
        // writeQuery("select 1",150277601);

        // writeQuery("select 1",127008);
        // writeQuery("select 1",127009);


        // writeQuery("select 1",4757056);
        // writeQuery("select 1",4757057);

        // writeQuery("select 1",9);
        // writeQuery("select 1",69347);

        // writeQuery("select 1",75777);

        // writeQuery("select 1",73088);
        // writeQuery("select 1",73089);

        // writeQuery("select 1",111843);
        // writeQuery("select 1",111844);

        // writeQuery("select 1", 800000);
        // writeQuery("select 1", 900000);
        // writeQuery("select 1",1000000);
        // writeQuery("select 1",198020);
        // writeQuery("select 1",684);
        // writeQuery("select 1",1366);
        // writeQuery("select 1",1367);

        //writeQuery("select 0", 2945 );
        //writeQuery("select 1",304475629);


    }

    inline void alternating_bits() {
        uint64_t const n = N_bits;
        uint64_t const W = 0xAAAAAAAAAAAAAAAA; //= 0b10101010...
        for (int w = 0; w < N_words; w++) {
            bits[w] = W;
        }
        N_ones = n/2;
        N_zeros = n/2;
        cout << "creating alternating bits " << endl;

    }


    inline void every_kth() {
        uint64_t const n = N_bits;
        uint64_t W = 0;
        for (int i=0; i<n; i++) {
            if (i%every_kth_spacing == 0) {
                W |= 1ULL << (i%64);
                N_ones+=1;
            }
            else {
                N_zeros+=1;
            }
            if ((i+1)%64==0) {
                uint64_t w=i/64;
                bits[w]=W;
                W=0;
            }
        }
        cout << "creating bits every " << every_kth_spacing << " k-th place" << endl;
    }

    inline double get_duration(tp tStart, tp tEnd) {
        double dur = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tStart).count()/1000.0; //ms
        dur = ((int)(dur * 10.0))/10.0; //round to one digit after decimal point
        return dur;
    }

    inline void bimodal_geometric() {
        double d_avg = synth_01ratio + 1; //ratio 9:1 leads to the average distance of 10
        double d_short = std::max(1.0,(double)d_avg/(double)bimodal_factor);
        //double d_short = std::max(1.0,d_avg/2.0);
        double d_long  = bimodal_factor *d_avg;
        double p_short = (d_long - d_avg)/(d_long - d_short);

        //std::random_device rd;
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> uniform01;
        std::geometric_distribution<int> geom_short(1/d_short);
        std::geometric_distribution<int> geom_long(1/d_long);

        uint64_t const n = N_bits;
        uint64_t i = 0;
        cout << "Starting bimodal generation" << endl;

        //few 1-bits, first remember and then insert individually
        if (synth_01ratio >= 16) {
            // Store positions where 1s should be placed
            std::vector<uint64_t> bit_positions(0);
            //cout << "reserve" << endl;
            bit_positions.reserve(1.2* ((double)n / d_avg));
            while (i < n) {
                double p_choice = uniform01(gen);
                uint64_t di;
                if (p_choice < p_short) {
                    di = geom_short(gen) + 1;
                } else {
                    di = geom_long(gen) + 1;
                }
                i += di;
                if (i < n) {
                    bit_positions.push_back(i);
                }
            }

            // Insert the 1s
            for (uint64_t pos : bit_positions) {
                bits[pos / 64] |= 1ULL << (pos % 64);
            }
            N_ones = bit_positions.size();
            N_zeros = n - N_ones;
        }
        //rather dense 1-bits, create and insert on the fly, without remembering
        else {
            uint64_t W = 0;
            uint64_t di = 0;
            while (i<n) {
                double p_choice = uniform01(gen);
                if (p_choice < p_short) {
                    di = geom_short(gen)+1;
                } else {
                    di = geom_long(gen)+1;
                }
                uint64_t j=i+di;
                while (i<j) {
                    if ((i+1)%64==0) {
                        bits[i/64]=W;
                        int pop = std::popcount(W);
                        N_ones += pop;
                        N_zeros += 64-pop;
                        W=0;
                    }
                    i++;
                }
                W |= 1ULL << (i%64);
            }
        }

        //cout << "finished bimodal distribution" << endl;
    }


    inline void patches() {
        std::random_device rd;
        //uint64_t seed = 12345;
        std::mt19937_64 rn(seed);
        //Within one patch, allow up to 2^max_wordcount_exp words to be set at once
        constexpr int max_wordcount_exp = 13;
        constexpr uint64_t wordcount_MASK = (1ULL<<max_wordcount_exp)-1;
        std::uniform_real_distribution<float> rand_p(0.0, 0.5);
        //Each word can be filled with 0 bits, 1 bit, 10 bits or 60 bits
        constexpr double d1 = 1.0/64;
        constexpr double d10 = 10.0/64;
        constexpr double d60 = 60.0/64;
        constexpr uint64_t W1 = 1ULL;
        constexpr uint64_t W10 = 0x3FF;
        constexpr uint64_t W60 = 0xFFFFFFFFFFFFFFF;

        //Each fill type receives probabilities such that the overall desired density is approximated
        const double dgoal = 1.0/(synth_01ratio+1);
        double p1 = dgoal/d1;
        double p2 = dgoal/d10;
        double p3 = dgoal/d60;
        double psum = p1 + p2 + p3;
        //normalize
        p1 = p1/psum;
        p2 = p2/psum;
        p3 = p3/psum;

        p1 = 0.01;
        p2 = 0.1;
        p3 = 0.39;
        // cout << "p1 = " << p1 << endl;
        // cout << "p21 = " << p2 << endl;
        // cout << "p3 = " << p3 << endl;
        // cout << "p1+p2+p3= " << p1+p2+p3 << endl;
        vector<double> P = {p1, p1+p2, p1+p2+p3};
        uint64_t w  = 0;
        uint64_t COIN_BITS = rn();
        uint64_t coinflips_used = 0;
        while (w<N_words) {
            bool fill = (COIN_BITS>>coinflips_used)&1ULL;
            coinflips_used++;
            if (coinflips_used==64) {
                COIN_BITS = rn();
                coinflips_used = 0;
            }
            //Sample how many words should be in this patch
            //We Want to bias this towards having often lower counts and rarely high counts
            //A quick heuristic to achieve this is to take the minimum of two random numbers
            uint64_t wc_sample1 = rn() & wordcount_MASK;
            uint64_t wc_sample2 = rn() & wordcount_MASK;
            uint64_t wordcount = std::min(wc_sample1,wc_sample2);
            if (!fill) {
                //skip filling
                w+=wordcount;
                continue;
            }
            double p = rand_p(rn);
            uint64_t W = 0;
            if (p<P[0]) {
                W = W1 << coinflips_used;
            } else if (p<P[1]) {
                W = W10;
            } else  {
                W = W60;
            }
            while (wordcount>0 && w<N_words) {
                wordcount--;
                bits[w++] = W;
            }
        }

        for(auto W : bits) {
            int pop = std::popcount(W);
            N_ones += pop;
            N_zeros += 64-pop;
        }
    }


    static constexpr uint64_t SIN_LOOKUP_SIZE = 2048;
    static constexpr double TWO_PI = 6.28318530718f;

    // Generate 1024-point sine lookup table
    static const std::array<double, SIN_LOOKUP_SIZE> sin_lut = []() {
        std::array<double, SIN_LOOKUP_SIZE> lut{};
        for (int i = 0; i < SIN_LOOKUP_SIZE; ++i) {
            lut[i] = std::sin((TWO_PI * i) / SIN_LOOKUP_SIZE);
            // cout << "lut[" << i << "] = " << lut[i] << endl;
        }
        return lut;
    }();

    // Fast sine lookup with linear interpolation
     __attribute__ ((always_inline)) inline double fast_sin(uint64_t i, double freq) {
         double phase = i * freq;  // Unwrapped phase (cycles)
         double index_f = phase * SIN_LOOKUP_SIZE;  // Map to LUT index space
         uint64_t index = static_cast<uint64_t>(index_f) % SIN_LOOKUP_SIZE;
         // if (sin_lut[index]==0) {
             // cout << i <<  " " << phase << " " <<index_f << " " << static_cast<uint64_t>(index_f) << endl;
         // }
         return sin_lut[index];
    }

    static constexpr int num_sinuses = 4;
    inline std::array<double, num_sinuses> make_freqs(const std::array<uint64_t, num_sinuses>& periods) {
        std::array<double, num_sinuses> result{};
        for (std::size_t i = 0; i < result.size(); ++i)
            result[i] = 1 / (double)periods[i];
        return result;
    }

    // inline std::array<uint64_t,num_sinuses> periods = {100000,  97777, 7678, 33}; //exemplary first settings
    inline std::array<uint64_t,num_sinuses> periods = {71,90,81,18}; //exemplary first settings
    static constexpr std::array<double,num_sinuses> weights = {   0.4,   0.3,  0.2, 0.1}; //should sum to 1
    inline std::array<double,num_sinuses> freqs   = make_freqs(periods);

    __attribute__ ((always_inline)) inline double get_sin_sum(uint64_t i) {
        double sum = 0;
        for (size_t k=0; k<num_sinuses; k++) {
            sum += weights[k] * fast_sin(i,freqs[k]);
            // sum += weights[k] * std::sin(i*freqs[k]);
        }
        return sum;
    }

    inline std::string get_sinus_periods_as_string() {
        std::string s;
        for (uint64_t p : periods) {
            s += std::to_string(p) + ",";
        }
        return s;
    }


    inline void print_periods_overview() {
        cout << "  Periods     Weights" << endl;
        for (int k=0; k<num_sinuses; k++) {
            cout << setw(10) << periods[k] << setw(10) << weights[k] << endl;
        }
    }

    inline void generate_periods() {
        //Generate four random periods spanning multiple orders of magnitude
        std::mt19937 rn(seed);
        // int min_order = 1; //Min period: 10^1
        // int max_order = 6; //Max period: 10^6
        std::vector<int> magnitudes = {0,0,0,0,0,1,1,1,1,2,2,2,3,3,4,4,5,6}; //Choose periods from 10^0 to 10^6, with bias towards lower orders, more similar to independent bits
        std::uniform_int_distribution<int> rand_order_index(0, magnitudes.size()-1);

        if (curr_instance>0) { //keep exemplary setting for first instance for easier comparison during development
            for (int i = 0; i<num_sinuses; i++) {
                //Choose an order of magnitude
                int rand_magnitude = magnitudes[rand_order_index(rn)];
                //Choose a random number within this magnitude
                uint64_t lower = static_cast<int>(std::pow(10, rand_magnitude));
                uint64_t upper = static_cast<int>(std::pow(10, rand_magnitude + 1)) - 1;
                if (rand_magnitude==0) { //so few options, extend them a bit
                    lower = 4; //minimum sensible period
                    upper = 20; //allow some more variety at very low periods
                }
                std::uniform_int_distribution<uint64_t> period_distr(lower, upper);
                periods[i] = period_distr(rn);
            }
            freqs = make_freqs(periods);
        }

        if (synth_01ratio==1) {
            print_periods_overview();
        }
        // double weight_sum = 0;
        // for (int k=0; k<num_sinuses; k++) {
            // weight_sum += weights[k];
        // }
        // if (weight_sum < 0.99 || weight_sum > 1.01) {
            // cout << "\n \n Warning: Sinusoidal weights don't sum to 1" << endl;
        // }
    }


    inline double find_sinus_threshhold() {
        constexpr uint64_t BUCKETS_PER_HALF = 500;
        constexpr uint64_t SAFETY_PADDING = 6;
        constexpr uint64_t TOTAL_BUCKETS = 2*BUCKETS_PER_HALF + 2*SAFETY_PADDING;
        constexpr uint64_t MID = TOTAL_BUCKETS/2;
        constexpr double float_BPH = BUCKETS_PER_HALF;
        std::array<uint64_t, TOTAL_BUCKETS> histogram = {};
        for (uint64_t i=0; i<N_bits; i++) {
            double sum = get_sin_sum(i); // in [-1,1]
            // cout << sum << endl;
            int index = sum*float_BPH; // in [-NUM_BUCKETS, NUM_BUCKETS]
            int centered_index = index + MID; //center sum=0 in the middle of the histogram
            histogram[centered_index]++;
            // if (centered_index==MID) {
                // cout << sum << endl;
            // }
            // if (i%(1ULL<<29)==0) {
                // cout << i << " " << i/(double)N_bits << endl;
            // }

        }
        const uint64_t goal_1bits = N_bits / (synth_01ratio+1);
        uint64_t accumulated_1bits = 0;
        double thresh = 0;
        // cout << "goal_1bits = " << goal_1bits << endl;
        // for (int i=TOTAL_BUCKETS-1; i>=0; i--) {
            // cout << i << ": " <<histogram[i] << endl;
        // }
        for (int i=TOTAL_BUCKETS-1; i>=0; i--) {
            accumulated_1bits+=histogram[i];
            if (accumulated_1bits >= goal_1bits) {
                thresh = (i-MID) /float_BPH;
                // cout << "hit tresh at i=" << i << endl;
                // cout <<"threshhold = " << thresh << " to include density "  << accumulated_1bits/(double)N_bits << endl;
                break;
            }
        }
        return thresh;
    }

    inline double sin_threshhold = 0;

    inline void sinusoidal() {
        generate_periods();
        // cout << "Finding threshhold" << endl;
        double thresh = find_sinus_threshhold();
        sin_threshhold = thresh;
        // cout << "Generating sinusoidal distribution" << endl;
        for (uint64_t w=0; w<N_words; w++) {
            uint64_t W = 0;
            const uint64_t i0 = w*64;
            for (uint64_t i=0; i<64; i++) {
                const double sinsum = get_sin_sum(i+i0);
                if (sinsum > thresh) {
                    W |= (1ULL<<i);
                }
            }
            bits[w] = W;
            // if (w%(1ULL<<23)==0) {
                // cout << w*64 << endl;
            // }
        }
        for(auto W : bits) {
            int pop = std::popcount(W);
            N_ones += pop;
            N_zeros += 64-pop;
        }
    }

    inline void gap_gadget_instance() {
        //We first generate a uniformly indepent bit instance and then zero out a gap in the middle.
        start_gidney_fast_independent_bits();
        uint64_t gap_words = curr_gap_size_bits / 64;
        cout << "gap_words = " << gap_words << endl;
        uint64_t dense_words_before_gap = (N_words - gap_words)/2;
        for (uint64_t w=0; w<gap_words; w++) {
            bits[w + dense_words_before_gap] = 0; //clear away any set bits
        }
        //Update the bit counts. Could be done more elegantly centralized, but leave for now this.
        N_ones = 0;
        N_zeros = 0;
        for(auto W : bits) {
            int pop = std::popcount(W);
            N_ones += pop;
            N_zeros += 64-pop;
        }


    }

    // inline void every64() {
    //     uint64_t const n = N_bits;
    //     uint64_t const W = 1; //= 0b10101010...
    //     for (int w = 0; w < N_words; w++) {
    //         bits[w] = W;
    //         N_ones += 1;
    //         N_zeros += 63;
    //     }
    //     cout << "creating every64 bits " << endl;
    //
    // }

    inline void manual_bits() {
        uint64_t n = N_bits;

        N_ones = 0;
        N_zeros = 0;

        //repeating:
        // 2 blocks with alternating bits
        // 2 blocks with few 1-bits

        uint64_t const ALT = 0xAAAAAAAAAAAAAAAA; //= 0b10101010...

        for(int w=0; w< N_words; w++) {
            if((w/64)%2==0) {
                bits[w] = ALT;
                // cout << ALT << endl;
            } else {
                bits[w] = 0xF;
                // cout << 1ULL << endl;
                N_ones += 1;
            }
            int pop = std::popcount(bits[w]);
            N_ones += pop;
            N_zeros += 64 - pop;
        }


        // uint64_t w = 0;
        // uint64_t pad = 10;
        //
        // uint64_t GROUP = 200000;
        //
        // while(N_ones<GROUP) {
        //     bits[w] = 0b1000;
        //     N_ones += 1;
        //     N_zeros += 63;
        //
        //     w += pad;
        //     N_zeros += (pad-1)*64;
        // }
        // while(N_ones < 2*GROUP) {
        //     bits[w++] = 0xFFFFFFFFFFFFFFFF;
        //     N_ones += 64;
        // }
        // while(N_ones<3*GROUP) {
        //     bits[w] = 0b1000;
        //     N_ones += 1;
        //     N_zeros += 63;
        //
        //     w += pad;
        //     N_zeros += (pad-1)*64;
        // }
    }


    //Reverse the order of bits in a 64-bit word
    inline uint64_t reverse_single_word(uint64_t v) {
        //Taken from:
        //http://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel
        //and expanded to 64-bit

        // swap odd and even bits
        v = ((v >> 1) & 0x5555555555555555ULL) | ((v & 0x5555555555555555ULL) << 1);
        // swap consecutive pairs
        v = ((v >> 2) & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
        // swap nibbles
        v = ((v >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
        // swap bytes
        v = ((v >> 8) & 0x00FF00FF00FF00FFULL) | ((v & 0x00FF00FF00FF00FFULL) << 8);
        // swap 2-byte long pairs
        v = ((v >> 16) & 0x0000FFFF0000FFFFULL) | ((v & 0x0000FFFF0000FFFFULL) << 16);
        // swap 4-byte long pairs
        v = ((v >> 32) & 0x00000000FFFFFFFFULL) | ((v & 0x00000000FFFFFFFFULL) << 32);
        return v;
    }


    //Poppy wants the input 64-bit words in reversed order (relative to our standard)
    inline void create_64bit_reversed_bitvector() {
        bits_64bit_reversed = std::vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>>(N_words);
        for(uint64_t w=0; w<N_words; w++) {
            uint64_t const W = bits[w];
            bits_64bit_reversed[w] = reverse_single_word(W);
        }
    }


    inline void generate_bits() {
        if     (do_manual_bits){manual_bits();}
        else if(do_fractals)    {fractal();}
        else if(do_chunky)      {chunky();}
        else if(do_smooth)      {smooth();}
        else if(do_quicksmooth) {quick_smooth();}
        else if(do_gidneysmooth){start_gidney_fast_independent_bits();}
        else if(do_alternate)   {alternating_bits();}
        else if(do_every_kth)   {every_kth();}
        else if(do_bimodal)     {bimodal_geometric();}
        else if (do_patches)    {patches();}
        else if (do_sinus)      {sinusoidal();}
        else if (do_gap)        {gap_gadget_instance();}
        else {smooth();}

        if(invert_all_bits) {
            for(auto &W : bits) {
                W = ~W;
            }
            std::swap(N_ones, N_zeros);
        }
        eff_01ratio = N_zeros /(double)N_ones;

        //unordered_set<string> need_reversed_bitvector = {"poppy", "poppy2"};
        // if(BVTYPE=="poppy" || BVTYPE=="poppy2") {
            // cout << " ###   Reversing bit vector for Poppy standard ###" << endl;
            // create_64bit_reversed_bitvector();
            // cout << " ###   Reversing finished ###" << endl;
        // }

        // plotter::bitscatter();
    }


    //Generated random queries, with amounts specified by main file
    inline void random_queries() {
        //uint64_t seed = 12345;
        std::random_device rd;
        std::mt19937 gen(rd());

        std::uniform_int_distribution<uint64_t> randRank(0, N_bits-1); //Rank distribution
        std::uniform_int_distribution<uint64_t> randSel0(1, N_zeros-1);//Select_0 distribution
        std::uniform_int_distribution<uint64_t> randSel1(1, N_ones-1); //Select_1 distribution

        for(auto i=0; i<synth_accesses; i++) {
            writeQuery("access",randRank(gen));
        }
        for(auto i=0; i<synth_ranks/2; i++) {
            writeQuery("rank 0",randRank(gen));
            writeQuery("rank 1",randRank(gen));
        }
        for(auto i=0; i<synth_selects/2; i++) {
            writeQuery("select 0",randSel0(gen));
            writeQuery("select 1",randSel1(gen));
        }
        for(auto i=0; i<synth_select1s; i++) {
            writeQuery("select 1",randSel1(gen));
        }
        for(auto i=0; i<synth_select0s; i++) {
            writeQuery("select 0",randSel0(gen));
        }

    }

    inline void dummy_gap_queries() {
        //Generate 2^26 queries within [0, m/4], i.e. well within the left dense half, to flood the cache with them when queried
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<uint64_t> first_quarter(1, N_ones/4);
        // std::uniform_int_distribution<uint64_t> last_quarter((3*N_ones)/4, N_ones);
        for (int i=0; i<generated_dummy_queries; i++) {
            writeQuery("select 1",first_quarter(gen));
        }
    }

    inline void generate_queries() {
        queries = std::vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>>(0);
        if (do_gap) {
            dummy_gap_queries();
        }
        else if(do_random_queries) {
            random_queries();
        }
        else {
            manual_queries();
        }
    }

    inline double bit_density = 0;

    inline void generate_synthetic_input() {
        tp t0 = std::chrono::steady_clock::now();

        N_words = synth_bits/WORDSIZE;
        if(synth_bits%WORDSIZE!=0) N_words++;
        bits = std::vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>>(N_words,0);
        N_bits = N_words * 64;
        N_ones=0; N_zeros=0; dRead = 0;

        generate_bits();
        generate_queries();
        N_queries = queries.size();
        bit_density = N_ones / (double)N_bits;

        tp t1 = std::chrono::steady_clock::now();
        dRead = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count(); //s
    }
}

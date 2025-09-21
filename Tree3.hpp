
#pragma once
#include <cmath>
#include <filesystem>

#include "bv.hpp"


using namespace std; //clean up eventually

class Tree3 {


public:
	static constexpr bool SHOW_T3_BUILDING		=0;
	static constexpr bool SHOW_T3_VERBOSE_STATISTICS=0;
	static constexpr bool WRITE_RESULT_FILE		=0;
	static constexpr double build_space_cutoff_factor = 1.2;


	 /*
	  *   Layout:
	  *   T			[l1 j kappa]
	  *   M_s		[dl c | dl c | ... | rho_max h]
	  *   K_rho_max	[C_rho_min  ... C_rho_max]
	  *   B_rho	    [ddl ddl ddl ...]
	  */
	uint64_t L_exp = 16; // to be set from outside

	uint64_t a_exp;
	uint64_t b_exp;

	uint64_t a, b;
	uint64_t mod_a, mod_b;
	uint64_t MIDENTRIES;
	uint64_t MIDENTRIES_w_SENTINEL;
	uint64_t gamma_min;
	uint64_t s_min;
	uint64_t rho_min;
	uint64_t i_END_OF_BITVECTOR;

	uint64_t j_max, kappa_max, gamma_max, s_max;
	uint64_t rhomax_max, h_max, C_max;

	uint64_t count_topgroups;
	uint64_t count_midgroups;
	uint64_t count_botgroups;
	uint64_t estimated_packed_size;

	struct topgroup {
		uint64_t l1{0};
		uint64_t j{0};
		uint64_t kappa{0};
		uint64_t s{0};
		uint64_t gamma{0};
	};

	struct midgroup {
		std::vector<uint32_t> DL, c_local;
		uint64_t i_TG_next{0};
		uint64_t rhomax{0};
		uint64_t c_max{0};
		uint64_t kappa{0};
		uint64_t h{0};
		uint64_t s{0};
		uint64_t gamma{0};
	};


	uint64_t glob_l1_bits      = 0;
	uint64_t glob_j_bits       = 0;
	uint64_t glob_kappa_bits   = 0;
	uint64_t glob_topgroup_bits= 0;
	uint64_t glob_rhomax_bits = 0;
	uint64_t glob_h_bits      = 0;
	uint64_t glob_C_bits        = 0;
	// uint64_t glob_C_single_size= 0;
	// uint64_t glob_C_tuple_size = 0;

	// struct packed_sizes {
		// uint64_t l1_bits      = 0;
		// uint64_t j_bits       = 0;
		// uint64_t kappa_bits   = 0;
		// uint64_t topgroup_bits= 0;
		// uint64_t rho_max_bits = 0;
		// uint64_t h_bits       = 0;
		// uint64_t C_single_size   = 0;
		// uint64_t C_tuple_size = 0;
		// uint64_t C_slots	  = 0;
	// };


	static constexpr size_t LOGNMAX= 64;
	static constexpr uint64_t SCANNED_BEYOND_ALPHA = 10000000000000000000;
	typedef std::vector<uint64_t, AlignedAllocator<uint64_t,ALIGNMENT>> aligned_vector;

	std::array<uint64_t, LOGNMAX> num_MID_s;
	std::array<uint64_t, LOGNMAX> num_K_rhomax;
	std::array<uint64_t, LOGNMAX> num_BOT_rho;

	std::array<uint64_t, LOGNMAX> num_MID_kappa;
	std::array<uint64_t, LOGNMAX> num_MID_gamma;

	std::vector<topgroup>                      tmp_TOP;
	std::array<std::vector<midgroup>, LOGNMAX> tmp_MID;
	std::array<std::vector<uint32_t>, LOGNMAX> tmp_KK;
	std::array<std::vector<uint32_t>, LOGNMAX> tmp_BOT;

	aligned_vector                      pack_TOP;
	std::array<aligned_vector, LOGNMAX> pack_MID;
	std::array<aligned_vector, LOGNMAX> pack_KK;
	std::array<aligned_vector, LOGNMAX> pack_BOT;


	void reset_tree() {
		//We only set a_exp and b_exp, the rest follows
		a = 1ULL << a_exp;
		b = 1ULL << b_exp;
		mod_a = a-1;
		mod_b = b-1;
		MIDENTRIES = a/b;
		//In theory we have a/b entries per midgroup, in practive we need one more sentinel at the end, for the l1'-l1 calculation
		MIDENTRIES_w_SENTINEL = MIDENTRIES + 1;
		gamma_min = rlog2(ALPHA);
		s_min     = rlog2(ALPHA);
		rho_min   = rlog2(ALPHA);
		//Already set some sensible value, in case rhomax_max is then never needed
		rhomax_max = rho_min;

		i_END_OF_BITVECTOR = N_bits;

		//Running counters to track the largest occurring value for each parameter, to determine their required size in bits
		j_max     = 0;
		kappa_max = 0;
		gamma_max = 0;
		s_max     = 0;
		h_max     = 0;
		C_max = 0;

		//Logging Statistics
		count_topgroups = 0;
		count_midgroups = 0;
		count_botgroups = 0;
		estimated_packed_size = 0;

		//track the current number of s-sized midgroups, tuples in K_rhomax, and rho-sized botgroups
		num_MID_s = {};
		num_K_rhomax = {};
		num_BOT_rho = {};

		//For statistics: Track the detailed s = gamma + kappa composition
		num_MID_kappa = {};
		num_MID_gamma = {};

		glob_l1_bits      = rlog2(N_bits >> L_exp);


		//Temporary vectors to extract the tree-data in an organized manner
		std::vector<topgroup>().swap(tmp_TOP);
		for (auto& M : tmp_MID) {
			for (auto& mb : M) {
				std::vector<uint32_t>().swap(mb.DL);
				std::vector<uint32_t>().swap(mb.c_local);
				// std::vector<uint32_t>().swap(mb.C_global);
			}
			std::vector<midgroup>().swap(M);
		}
		for (auto& K : tmp_KK)  { std::vector<uint32_t>().swap(K); }
		for (auto& B : tmp_BOT) { std::vector<uint32_t>().swap(B); }

		//Final vectors to put in the packed tree data
		aligned_vector().swap(pack_TOP);
		for (auto& pM : pack_MID ) {aligned_vector().swap(pM);}
		for (auto& pK : pack_KK )  {aligned_vector().swap(pK);}
		for (auto& pB : pack_BOT)  {aligned_vector().swap(pB);}

	}


	//Get the number of bits needed to store values up to x
	//Detail: We always allocate at least 1 bit, even if "no" data needs to be represented. Makes construction easier, increases space only marginally
	static uint64_t rlog2(uint64_t x) {
		if(x==0) return 1;
		return 64 - _lzcnt_u64(x);
	}

	//Faster version for query, assumes x>0 to get rid of branch/check. During query we only take log of values >0 anyways
	static uint64_t quick_rlog2(const uint64_t x) {
		return 64 - _lzcnt_u64(x);
	}


	static uint64_t count_1s_from_position_and_higher(const uint64_t p, const uint64_t W) {
		return std::popcount(W >> p);
    }

	static uint64_t count_1s_before_position(const uint64_t p, const uint64_t W) {
    	if(p==0) {return 0; }
		return std::popcount(W << (64-p));
    }

	static uint64_t get_kth_1_position(const int k, const uint64_t W) {
		return _tzcnt_u64(_pdep_u64(1ULL << (k-1), W));
    }


	[[nodiscard]] uint64_t scan_forward_to_nth_1(const uint64_t i, uint64_t n, const uint64_t alpha) const {
		if (i == i_END_OF_BITVECTOR) {
			return i_END_OF_BITVECTOR;
		}
		const uint64_t l1_start = i >> L_exp;
		uint64_t w = i/64;
		uint64_t count = count_1s_from_position_and_higher(i%64, bits[w]);
		w++;
		while (count < n) {
			const uint64_t l1_curr = w >> (L_exp - 6);
			if (l1_curr - l1_start >= alpha) [[unlikely]] {
				return SCANNED_BEYOND_ALPHA;
			}
			if (w >= N_words) [[unlikely]] {
				return i_END_OF_BITVECTOR;
			}
			count += std::popcount(bits[w]);
			w++;
		}
		w--;
		count -= std::popcount(bits[w]);
		const uint64_t k = n - count;
		const uint64_t p = get_kth_1_position(k, bits[w]);
		return w*64 + p;
	}

		void get_and_store_BOT_indices(const uint64_t i_MIDENTRY, const uint64_t l1_ME, const uint64_t rho) {
	    uint64_t w			= i_MIDENTRY/64;
    	uint64_t p_start	= i_MIDENTRY%64;
    	int k_start = 1 + count_1s_before_position(p_start,bits[w]);
    	int count = 0;
		// std::vector<uint64_t> I(b);
		const uint64_t w_max = N_words;
		auto &B_rho = tmp_BOT[rho];
		// aligned_vector &pack_B_rho = pack_BOT[rho];
		// uint64_t i_Brho = curr_BOT_rho_index[rho];
    	while(count < b) {
    		if (w > w_max) {
    			//reached end of bitvector, exit
    			// uint64_t i = -1;
    			// store_single_ddl_BOT(i, l1_MG, B);
    			return;
    		}
			const uint64_t W = bits[w];
			int pop = std::popcount(W);
    		const uint64_t l1 = (w*64) >> L_exp;
    		const uint64_t ddl = l1- l1_ME;

			// _mm_prefetch(reinterpret_cast<const char*>(&B.back()), _MM_HINT_NTA); // Prefetch highest tree level into cache

			for(int k=k_start; k<=pop && count < b; k++) {
				// uint64_t i = w*64 + get_kth_1_position(k,W);
				// store_single_ddl_BOT(i, l1_MG, B);
				// i_Brho = insert_precise(i_Brho, ddl, rho, pack_B_rho);
				B_rho.emplace_back(ddl);
				count++;
			}
			k_start=1;
			w++;
    	}
		// curr_BOT_rho_index[rho] = i_Brho;
    }


	void build_tree() {
		uint64_t i_TOP = scan_forward_to_nth_1(0,1,N_bits);
		//Build the tree top-down
		while (i_TOP < i_END_OF_BITVECTOR) {
			i_TOP = new_topgroup(i_TOP);
		}
		new_topgroup(i_END_OF_BITVECTOR); //sentinel
	}


	uint64_t new_topgroup(const uint64_t i_TOP) {
		topgroup TOP;
		TOP.l1 = i_TOP >> L_exp;
		//Scan forwards, but abort if we pass beyond ALPHA superblocks, in that case we need to zoom in more
		uint64_t i_TOP_next = scan_forward_to_nth_1(i_TOP, a+1, ALPHA);
		if (i_TOP_next == SCANNED_BEYOND_ALPHA) {
			//Bits spanned more than ALPHA superblocks, split top-group into a midgroup
			midgroup MID = new_midgroup(i_TOP);
			if (is_tree_too_large()) {
				return i_END_OF_BITVECTOR;
			}
			i_TOP_next = MID.i_TG_next;
			TOP.j      = num_MID_s[MID.s]++;
			TOP.kappa  = MID.kappa;
			TOP.s	   = MID.s;
			TOP.gamma  = MID.gamma;
			j_max      = std::max(TOP.j, j_max);
		}
		tmp_TOP.emplace_back(TOP);
		if (SHOW_T3_BUILDING )	print_TOP(TOP);
		if (is_tree_too_large()) {
			return i_END_OF_BITVECTOR;
		}
		count_topgroups++;
		//If we scan (a,b)-Parameter choices we also run into some bad ones, abort them quickly by detecting a too large growing tree
		estimated_packed_size += glob_l1_bits + rlog2(j_max) + 1;
		return i_TOP_next;
	}


	midgroup new_midgroup(const uint64_t i_MIDGROUP) {
		midgroup MGRP;
		//Make space for all (dl, c) tuples
		MGRP.DL.reserve(MIDENTRIES+10);
		MGRP.c_local.reserve(MIDENTRIES+10);
		//We need to know the rho-count within this midgroup for the local c-values
		std::array<uint64_t, LOGNMAX> local_BOT_rho_count = {};
		const uint64_t l1_MGRP = i_MIDGROUP >> L_exp;
		uint64_t i_MIDENTRY = i_MIDGROUP;
		uint64_t i_MIDENTRY_next = 0;
		//Create all the a/b midentries
		for (int ME_nr=0; ME_nr < MIDENTRIES; ME_nr++) {
			const uint64_t l1_MIDENTRY = i_MIDENTRY >> L_exp;
			const uint64_t dl    = l1_MIDENTRY - l1_MGRP;
			//Scan completely until the next midgroup (to the b+1-th bit), here there is no point in stopping earlier
			i_MIDENTRY_next = scan_forward_to_nth_1(i_MIDENTRY, b+1, N_bits);
			const uint64_t l1_MIDENTRY_next = i_MIDENTRY_next >> L_exp;
			const uint64_t r_prime    = l1_MIDENTRY_next - l1_MIDENTRY;
			uint64_t c = 0;
			if (r_prime >= ALPHA) {
				//Midentry splits into the bot-level
				const uint64_t rho = rlog2(r_prime);
				//Offset values to this midgroup need rho bits space each
				c = local_BOT_rho_count[rho]++;
				//We need to know the largest rho in this midgroup, that determines how many C-counters need to be stored in K
				MGRP.rhomax = std::max(MGRP.rhomax, rho);
				//c_max determines how many bits kappa we need to allocate for each tuple
				MGRP.c_max = std::max(MGRP.c_max, c);
				//Store all the single sample offsets ddl
				get_and_store_BOT_indices(i_MIDENTRY, l1_MIDENTRY, rho);
				estimated_packed_size += b * rho;
				//Remember globally the largest occuring rhomax over all midgroups
				rhomax_max = std::max(MGRP.rhomax, rhomax_max);
				count_botgroups++;
			}
			MGRP.DL.emplace_back(dl);
			MGRP.c_local.emplace_back(c);
			i_MIDENTRY = i_MIDENTRY_next;
			if (is_tree_too_large()) {
				return midgroup{};
			}
		}
		//After b midgroups we need one additional sentinel midgroup
		//This is to such that there exists a "next" l1' value to read when accessing the very last midentry per midgroup
		const uint64_t l1_MIDG_sentinel = i_MIDENTRY >> L_exp;
		const uint64_t dl_sentinel      = l1_MIDG_sentinel - l1_MGRP;
		const uint64_t c_sentinel = 0;
		MGRP.DL.emplace_back(dl_sentinel);
		MGRP.c_local.emplace_back(c_sentinel);

		//We need pointers into the K-level if there is at least one split in this midgroup
		if (MGRP.rhomax > 0) {
			MGRP.h = num_K_rhomax[MGRP.rhomax];
		} else {
			//can leave BLCK.h at zero, because we will never need to access K
		}
		//Update that we now will create one new tuple in K_rhomax
		num_K_rhomax[MGRP.rhomax]++;
		h_max = std::max(MGRP.h, h_max);
		store_Cs_in_K(MGRP.rhomax);

		//Now that the midgroup is finished, update the count of botgroups with size rho
		for (int rho=0; rho<LOGNMAX; rho++) {
			num_BOT_rho[rho] += local_BOT_rho_count[rho];
		}

		//calc how many bits we need to reference this midgroup
		MGRP.kappa = rlog2(MGRP.c_max);
		MGRP.i_TG_next = i_MIDENTRY_next;
		const uint64_t l1_TG_next = MGRP.i_TG_next >> L_exp;
		const uint64_t r = l1_TG_next - l1_MGRP;
		MGRP.gamma = rlog2(r);
		MGRP.s = MGRP.gamma + MGRP.kappa;
		tmp_MID[MGRP.s].push_back(MGRP);

		// more counters
		gamma_max = std::max(MGRP.gamma, gamma_max);
		kappa_max = std::max(MGRP.kappa, kappa_max);
		s_max     = std::max(MGRP.s, s_max);
		num_MID_gamma[MGRP.gamma]++;
		num_MID_kappa[MGRP.kappa]++;
		count_midgroups++;
		if constexpr (SHOW_T3_BUILDING) {
			print_rho_count(local_BOT_rho_count, MGRP.rhomax);
		}
		estimated_packed_size += MIDENTRIES_w_SENTINEL * MGRP.s + rlog2(rhomax_max) + rlog2(h_max);
		return MGRP;
	}

	void store_Cs_in_K(const uint64_t rho_max) {
		std::vector<uint32_t> &K = tmp_KK[rho_max];
		for (uint64_t rho=rho_min; rho <= rho_max; rho++) {
			//Store the counts from *before* we entered this midgroup
			uint64_t C_rho = num_BOT_rho[rho];
			K.emplace_back(C_rho);
			//remember the globally largest C, this determines the space we need to allocate for all C's
			C_max = std::max(C_rho, C_max);
		}
		estimated_packed_size += (rho_max - rho_min + 1) * rlog2(C_max);
	}


	// void store_C_glob(const uint64_t rho_max) {
	// 	std::vector<uint32_t> &K = tmp_KK[rho_max];
	// 	for (uint64_t rho=rho_min; rho <= rho_max; rho++) {
	// 		uint64_t C_rho = counter_BOT_rho[rho];
	// 		K.emplace_back(C_rho);
	// 		C_rho_max = std::max(C_rho, C_rho_max);
	// 	}
	// }


	// __attribute__((always_inline)) static void store_single_ddl_BOT(const uint64_t i, const uint64_t l1_MG_start, std::vector<uint32_t> &B) {
	// 	uint64_t l1 = i >> L1_exp;
	// 	uint64_t ddl = l1- l1_MG_start;
	// 	B.emplace_back(ddl);
	// }

	// void store_BOT_ddl(const uint64_t rho, const std::vector<uint64_t> &I, const uint64_t l1_MG) {
	// 	std::vector<uint32_t> &B = tmp_BOT[rho];
	// 	for (uint64_t i : I) {
	// 		if (i==-1) {
	// 			break; //reached end of bitvector
	// 		}
	// 		uint64_t l1 = i >> L1_exp;
	// 		uint64_t ddl = l1 - l1_MG;
	// 		B.emplace_back(ddl);
	// 	}
	// }

	// void update_BOT_rho_count(std::array<uint64_t, LOGNMAX> &local_BOT_rho_count) {
	// }

	void determine_global_sizes() {
		glob_l1_bits      = rlog2(N_bits >> L_exp);
		glob_j_bits       = rlog2(j_max);
		glob_kappa_bits   = rlog2(kappa_max);
		glob_topgroup_bits = glob_l1_bits + glob_j_bits + glob_kappa_bits;
		glob_rhomax_bits = rlog2(rhomax_max);
		glob_h_bits       = rlog2(h_max);
		glob_C_bits       = rlog2(C_max);
		// uint64_t C_slots   = glob_rho_max_max - rho_min + 1;
		// glob_C_tuple_size  = glob_C_single_size * C_slots;
	}




	void pack_tree() {
		determine_global_sizes();
		uint64_t i=0;
		for (auto tg : tmp_TOP) {
			//have local order [l1 j kappa], i.e. l1 in the higher bits per tuple
			i = insert_precise(i, tg.kappa, glob_kappa_bits,pack_TOP);
			i = insert_precise(i, tg.j,		glob_j_bits,	pack_TOP);
			i = insert_precise(i, tg.l1,	glob_l1_bits,	pack_TOP);
		}


		for (int s=0; s<LOGNMAX; s++) {
			i = 0;
			uint64_t midgroup_count = 0;
			aligned_vector &pack_Ms = pack_MID[s];
			for (midgroup mgrp : tmp_MID[s]) {
				//store a/b midentries, each containing the info [dl, c]
				for (int g_nr = 0; g_nr < MIDENTRIES_w_SENTINEL; g_nr++) {
					uint64_t dl = mgrp.DL[g_nr];
					uint64_t c  = mgrp.c_local[g_nr];
					//to have local order [dl c] with dl in the higher bits per tuple, must insert c first
					i = insert_precise(i, c , mgrp.kappa, pack_Ms);
					i = insert_precise(i, dl, mgrp.gamma, pack_Ms);
				}

				// store h first, such that we can imagine the local order as [rhomax h], with rhomax in the higher bits per tuple
				i = insert_precise(i, mgrp.h     , glob_h_bits     , pack_Ms);
				i = insert_precise(i, mgrp.rhomax, glob_rhomax_bits, pack_Ms);
				midgroup_count++;
			}
		}

		for (int rhomax=0; rhomax<LOGNMAX; rhomax++) {
			i=0;
			aligned_vector &pack_K_rhomax = pack_KK[rhomax];
			for (uint64_t C : tmp_KK[rhomax]) {
				i = insert_precise(i, C, glob_C_bits, pack_K_rhomax);
			}
		}

		for (int rho=0; rho<LOGNMAX; rho++) {
			i=0;
			aligned_vector &pack_B_rho = pack_BOT[rho];
			for (uint64_t ddl : tmp_BOT[rho]) {
				i = insert_precise(i, ddl, rho, pack_B_rho);
			}
		}

	}


	__attribute__((always_inline)) uint64_t static insert_precise(uint64_t i, uint64_t val, uint64_t size, aligned_vector &Data) {
		uint64_t w = i/64;
		uint64_t p = i%64;
		while (Data.size() < w+4) { //always keep the vector large enough
			Data.emplace_back(0);
		}
		Data[w] |= (val << p);
		if (p+size > 64) {
			Data[w+1] |= (val >> (64-p));
		}
		return i+size;
	}


	bool is_tree_too_large() const {
		return estimated_packed_size >  opt_packed_size * build_space_cutoff_factor;
	}


	uint64_t get_packed_size_in_bits() const {
		uint64_t words = 0;
		words += pack_TOP.size();
		for (auto &M : pack_MID) {words += M.size();}
		for (auto &K : pack_KK)  {words += K.size();}
		for (auto &B : pack_BOT) {words += B.size();}
		return words*64;
	}

	void save_treestats() const {
		treestats.a_exp = a_exp;
		treestats.b_exp = b_exp;
		treestats.topgroups = count_topgroups;
		treestats.midgroups = count_midgroups;
		treestats.botgroups = count_botgroups;
		treestats.final_tree_size_bits = opt_packed_size;
		treestats.final_tree_size_factor = opt_packed_size/(double)N_bits;
	}



	template<typename... Args>
	void LOG(Args&&... args) {
		if constexpr (SHOW_T3_BUILDING) {
			(std::cout << ... << std::forward<Args>(args)) << std::endl;
		}
	}

	static void print_TOP(const topgroup &TG) {
		// cout << "Printing TOP" << endl;
		// for (topgroup TG : vec_TOP) {
		std::cout << "l1 " << std::setw(10) << TG.l1 << "  j" << TG.j << "  k" << TG.kappa << "  g" << TG.gamma << "  s" << TG.s << std::endl;
		// }
	}


	void print_rho_count(const std::array<uint64_t, LOGNMAX>& rho_count, const uint64_t rho_max) const {
		for (uint64_t rho=rho_min; rho<=rho_max; rho++) {
			cout << setw(38) << "C_" << rho << ": " << rho_count[rho] << endl;
		}
	}

	void print_verbose_stats() {
		int val_width = 4;
		cout << endl;
		cout << "alpha = " << ALPHA << endl;
		cout << "min log(r) = log(alpha) = " << rho_min << endl;
		cout << "a " << setw(7) << a << "  (2^" << a_exp << ")" << endl;
		cout << "b " << setw(7) << b << "  (2^" << b_exp << ")" << endl;
		cout << "topgroups " << setw(9) << tmp_TOP.size() << endl;
		cout << "midgroups " << setw(9) << count_midgroups << endl;
		cout << "midentries " << setw(9) << count_midgroups * MIDENTRIES << endl;
		cout << "botgroups  " << setw(9) << count_botgroups << endl;
		cout << endl;
		for (uint64_t gamma=gamma_min; gamma<=gamma_max; gamma++) {
			cout << "    gamma " << setw(2) << gamma << " : " <<  setw(val_width) << num_MID_gamma[gamma] << endl;
		}
		cout << endl;
		for (uint64_t kappa=0; kappa<=kappa_max; kappa++) {
			cout << "    kappa " << setw(2) << kappa << " : " << setw(val_width) << num_MID_kappa[kappa] << endl;
		}
		for (uint64_t s=s_min; s<=s_max; s++) {
			cout << "    s " << setw(2) << s << " : " << setw(val_width) << num_MID_s[s] << endl;
		}
		for (uint64_t rho =0; rho <= rhomax_max; rho++) {
			if (rho > 0 && rho < rho_min) {
				continue;
			}
			cout << "    rho_max " << setw(2) << rho << " : " << setw(val_width) << num_K_rhomax[rho] << endl;
		}
		for (uint64_t rho = rho_min; rho <= rhomax_max; rho++) {
			cout << "    rho " << setw(2) << rho << " : " << setw(val_width) << num_BOT_rho[rho] << endl;
		}
		cout << "     l1 " << glob_l1_bits << " bits " << endl;
		cout << "      j " << glob_j_bits << " bits" << endl;
		cout << "  kappa " << glob_kappa_bits << " bits" << endl;
		cout << "top-all " << glob_topgroup_bits << " bits" << endl;
		cout << " -- " << endl;
		cout << "rhomax " << glob_rhomax_bits << " bits"<< endl;
		cout << "     h " << glob_h_bits << " bits"<< endl;
		cout << "     C " << glob_C_bits << " bits" << endl;
		cout << "estimated packed size " << estimated_packed_size << " bits" << endl;
	}


	double get_overhead_in_MB() const {
		uint64_t packed_bits = get_packed_size_in_bits();
		double MByte = packed_bits / (8.0*1024*1024);
		return MByte;
	}


	void print_ovhd_stats() const {
		uint64_t packed_bits = get_packed_size_in_bits();
		double tree_overhead = packed_bits /(double) N_bits;
		double MByte = get_overhead_in_MB();
		cout << "real      packed size " << packed_bits << " bits == " << MByte  << " MB" << endl;
		cout << "tree overhead factor  " << tree_overhead << endl;
	}


	void print_ovhd_result_line(uint64_t a_exp, uint64_t b_exp, double overhead) {
		string s = std::format("RESULT a_exp={} b_exp={} ovhd={:.12f} topg={} midb={} botb={} \n", a_exp, b_exp, overhead, count_topgroups, count_midgroups, count_botgroups);
		if (SHOW_T3_VERBOSE_STATISTICS) printf(s.c_str());
		result_lines.push_back(s);
	}

	void print_theo_ovhd_result_line(double a_exp, double b_exp, double range, double overhead) {
		string s1 = std::format("RESULT a_exp={:.3f} b_exp={:.3f} ovhd={} topg={} midb={} botb={} \n", a_exp, b_exp-range, overhead, count_topgroups, count_midgroups, count_botgroups);
		string s2 = std::format("RESULT a_exp={:.3f} b_exp={:.3f} ovhd={} topg={} midb={} botb={} \n", a_exp, b_exp+range, overhead, count_topgroups, count_midgroups, count_botgroups);
		if (SHOW_T3_VERBOSE_STATISTICS) {
			printf(s1.c_str());
			printf(s2.c_str());
		}
		result_lines.push_back(s1);
		result_lines.push_back(s2);
	}

    void write_ovhd_results_to_file() {
		std::ofstream file;
		string filename_out = "t3ovhds_" + BITGEN + "_r" + std::to_string(synth_01ratio);
		string foldername_out = "../../plotting/t3overheads/data/";
		if (BITGEN=="bimodal") {
			filename_out += "_bmf" + std::to_string(bimodal_factor);
		}
		filename_out += ".txt";
		file.open(foldername_out+filename_out);
		for(auto &s : result_lines) {
			file << s;
		}
		file.close();
	}


	struct theo_params {
		double a_star;
		double b_star;
		double a_star_exp;
		double b_star_exp;
	};


	theo_params get_strategy_params_a_b(bool show=false) const {
		switch (TREE3_STRATEGY) {
			case GET_CUBIC_THEORY_PARAMS:		return get_cubic_theory_params(show); break;
			case GET_LIFTED_CUBIC_THEORY_PARAMS:return get_lifted_cubic_theory_params(show); break;
			case GET_TFOCUS_THEORY_PARAMS:		return get_tfocus_theory_params(show); break;
			case GET_TOPHEAVY_PARAMS:			return get_top_level_avg_case_params(show); break;
			case GET_SMALLEST_TREE_PARAMS:		return get_lifted_cubic_theory_params(show); break; //as a reference value to scan around
			case GET_FASTEST_TREE_PARAMS:		return get_lifted_cubic_theory_params(show); break; //as a reference value to scan around
			default:
				cout << "Expected TREE3_PARAM_CHOICE={1,2,3,4,5,6}" << endl;
				return get_lifted_cubic_theory_params(show);
				break;
		}
	}

	//Clamp logarithms to one, to prevent negative values or division by zero in the theory estimations
	//Also models better that in practice we allocate at least one bit if there is something to allocate at all
	double clog2(double x) const {
		if (log2(x)>1)
			return log2(x);
		return 1;
	}

	theo_params get_cubic_theory_params(bool show=false) const {
		double L = 1<<L_exp;
		double s = N_bits /(ALPHA * L); //# max splits
		double alpha_over_logalpha = ALPHA / clog2(ALPHA);
		double n = N_bits;
		double m = N_ones;
		double b_star = std::pow(6*(m/n)* alpha_over_logalpha * L * log2(n/L), 0.333);
		double a_star = 0.5*b_star*b_star;
		double b_star_exp = log2(b_star);
		double a_star_exp = log2(a_star);
		if (show && SHOW_T3_VERBOSE_STATISTICS) {
			cout << "Cubic Theory Params" << endl;
			cout << "L=" << L << endl;
			cout << "L_exp" << log2(L) << endl;
			cout << "s=" << s << endl;
			cout << "alpha/log(alpha)=" << alpha_over_logalpha << endl;
			cout << "b*=" << b_star << endl;
			cout << "a*=" << a_star << endl;
			cout << "log(b*)=" << b_star_exp << endl;
			cout << "log(a*)=" << a_star_exp << endl;
		}
		return theo_params{a_star, b_star, a_star_exp, b_star_exp};
	}


	theo_params get_lifted_cubic_theory_params(bool show=false) const {
		double L = 1<<L_exp;
		double s = N_bits /(ALPHA * L); //# max splits
		double alpha_over_logalpha = ALPHA / clog2(ALPHA);
		double n = N_bits;
		double m = N_ones;
		double a_prefactor = std::pow((3/1.41) * (m/n) * alpha_over_logalpha, 0.666);
		double a = a_prefactor * L * std::pow(clog2(clog2(n)),0.333);
		double b = std::sqrt(2*a);
		double a_star_exp = log2(a);
		double b_star_exp = log2(b);
		if (show && SHOW_T3_VERBOSE_STATISTICS) {
			cout << "#\nLifted Cubic Theory Params" << endl;
			cout << "L=" << L << endl;
			cout << "L_exp" << log2(L) << endl;
			cout << "s=" << s << endl;
			cout << "alpha/log(alpha)=" << alpha_over_logalpha << endl;
			cout << "a*=" << a << endl;
			cout << "b*=" << b<< endl;
			cout << "log(a*)=" << a_star_exp << endl;
			cout << "log(b*)=" << b_star_exp << endl;
		}
		return theo_params{a, b, a_star_exp, b_star_exp};
	}


	theo_params get_tfocus_theory_params(bool show=false) const {
		double L = 1<<L_exp;
		double alpha_over_logalpha = ALPHA / clog2(ALPHA);
		double m = N_ones;
		double n = N_bits;
		double s = N_bits /(ALPHA * L);
		//First calculate a*
		double a_factor1 = std::sqrt(0.5 * (m/n) * alpha_over_logalpha * L) * log2(s);
		double a_factor2_core = std::sqrt(0.5*(n/m)) * log2(s) / std::sqrt(L);
		double a_factor2 = std::sqrt(clog2(a_factor2_core));
		double a_star = a_factor1 * a_factor2;
		//Use to calculate b*
		double b_factor1 = std::sqrt(2*(m/n)*alpha_over_logalpha * L);
		double b_factor2 = std::sqrt(clog2(a_star*n /(L*m)));
		double b_star = b_factor1 * b_factor2;
		double a_star_exp = log2(a_star);
		double b_star_exp = log2(b_star);
		double t = m / a_star;
		if (show && SHOW_T3_VERBOSE_STATISTICS) {
			cout << "t-focused Theory Params" << endl;
			cout << "L=" << L << endl;
			cout << "L_exp" << log2(L) << endl;
			cout << "s=" << s << endl;
			cout << "t=" << t << endl;
			cout << "alpha/log(alpha)=" << alpha_over_logalpha << endl;
			cout << "m/n=" << m/n << endl;
			cout << "a_factor_1=" << a_factor1 << endl;
			cout << "a_factor_2=" << a_factor2 << endl;
			cout << "b_factor_1=" << b_factor1 << endl;
			cout << "b_factor_2=" << b_factor2 << endl;
			cout << "a*=" << a_star << endl;
			cout << "b*=" << b_star << endl;
			cout << "log(a*)=" << a_star_exp << endl;
			cout << "log(b*)=" << b_star_exp << endl;
		}
		return theo_params{a_star, b_star, a_star_exp, b_star_exp};
	}

	theo_params get_top_level_avg_case_params(bool show=false) const {
		double L = 1<<L_exp;
		double alpha_over_logalpha = ALPHA / clog2(ALPHA);
		double m = N_ones;
		double n = N_bits;
		//Let a ~ (alpha L)/3, to have almost no expected splits
		double a_star = (m/n)*ALPHA*L/3;
		//And then take some simple choice for b*, like square root
		double b_star = std::sqrt(a_star);
		double a_star_exp = log2(a_star);
		double b_star_exp = log2(b_star);
		if (show && SHOW_T3_VERBOSE_STATISTICS) {
			cout << "#\nTop-level optimized" << endl;
			cout << "L      =" << L << endl;
			cout << "a*      " << a_star << endl;
			cout << "b*      " << b_star << endl;
			cout << "log(a*)=" << a_star_exp << endl;
			cout << "log(b*)=" << b_star_exp << endl;
		}
		return theo_params{a_star, b_star, a_star_exp, b_star_exp};

	}

	theo_params get_cycled_params(bool show=false) const {
		// auto theo = get_strategy_params_a_b();
		auto theo = get_lifted_cubic_theory_params();
		int start_shift_a_exp_cycle = -parameter_cycles/2;
		//With each benchmark run, curr_parameter_cycle is incremented by one, allowing us to test out different (a,b) values
		double a_exp_cycled = theo.a_star_exp + start_shift_a_exp_cycle + curr_parameter_cycle;
		double b_exp_cycled = a_exp_cycled / 2;
		a_exp_cycled = max(a_exp_cycled, 6.0); //Ensure some minimally useful values
		b_exp_cycled = max(b_exp_cycled, 3.0);
		double a_cycled = std::pow(2, a_exp_cycled);
		double b_cycled = std::pow(2, b_exp_cycled);
		if (show && SHOW_T3_VERBOSE_STATISTICS) {
			cout << "#\nCycled parameters" << endl;
			cout << "Total parmeter cycles="<<parameter_cycles<<endl;
			cout << "Curr cycle=" << curr_parameter_cycle << endl;
			cout << "log(a)=" << a_exp_cycled << endl;
			cout << "log(b)=" << b_exp_cycled << endl;
		}
		return theo_params{a_cycled, b_cycled, a_exp_cycled, b_exp_cycled};
	}

	void calc_theoretical_worstcase_overhead() {
		// auto params = get_strategy_params_a_b();
		auto params = get_lifted_cubic_theory_params();
		double L = 1<<L_exp;
		double s = N_bits /(ALPHA * L);
		double a = std::pow(2, params.a_star_exp);
		double b = std::pow(2, params.b_star_exp);
		double m = N_ones;
		double n = N_bits;
		double t = m / a;
		double G = std::min(s,t); //The number of midgroups is bounded both by the number of total splits and the number of topgroups, whichever is smaller
		double T = t *  (std::log2(n/ L) + std::log2(G) + clog2(clog2(a/b)));
		double M = G * (2*(a/b) * clog2(ALPHA*(s/G)) + clog2(clog2(ALPHA*(s/G))) + log2(G));
		double K = G * (std::log2(s/G)+1) * std::log2(s);
		double B = s * b * clog2(ALPHA);
		double S = T + M + K + B;
		double theo_worstcase_ovhd = S / double(N_bits);
		print_theo_ovhd_result_line(params.a_star_exp, params.b_star_exp, 0.1, theo_worstcase_ovhd);
		treestats.theo_a_exp = params.a_star_exp;
		treestats.theo_b_exp = params.b_star_exp;
		treestats.theo_worstcase_size_bits = S;
		treestats.theo_worstcase_size_factor = theo_worstcase_ovhd;
	}

	// uint64_t opt_tmp_size   = 1ULL<<60; //large.
	uint64_t opt_packed_size= 1ULL<<60;
	uint64_t opt_a_exp;
	uint64_t opt_b_exp;
	std::vector<string> result_lines = {};

	void test_specific_sampling_choice(uint64_t test_a_exp, uint64_t test_b_exp, bool are_theo_params=false) {
		a_exp = test_a_exp;
		b_exp = test_b_exp;
		reset_tree();
		build_tree();
		if (is_tree_too_large()) {
			string s = std::format("		a_exp {} b_exp {}: larger by factor >{}\n", a_exp, b_exp, build_space_cutoff_factor);
			if (SHOW_T3_VERBOSE_STATISTICS) printf(s.c_str());
			result_lines.push_back(s);
			return;
		}
		pack_tree();
		uint64_t pack_size = get_packed_size_in_bits();
		double overhead = pack_size/(double)N_bits;
		print_ovhd_result_line(a_exp, b_exp, overhead);
		if (are_theo_params) {
			treestats.theo_actual_size_bits = pack_size;
			treestats.theo_actual_size_factor = overhead;
		}
		//Remember these params if they give a new smallest tree
		if (pack_size < opt_packed_size) {
			opt_packed_size = pack_size;
			opt_a_exp = a_exp;
			opt_b_exp = b_exp;
		}
	}

	void find_best_sampling_choices() {
		result_lines.clear();
		opt_packed_size= 1ULL<<60; //large
		//get a*, b* from theory, as a first good guess
		auto strategy_params = get_strategy_params_a_b(true);
		uint64_t strat_a_exp = strategy_params.a_star_exp + 0.5; //round to nearest
		uint64_t strat_b_exp = strategy_params.b_star_exp + 0.5;
		//Do a first run with these theoretical values, to have some reference to quickly identify later test parameters that give way too large structures
		test_specific_sampling_choice(strat_a_exp, strat_b_exp, true);

		//scan different (a,b) parameters to find those which minimize the tree for the specific given bitvector
		if (TREE3_STRATEGY==GET_SMALLEST_TREE_PARAMS) {
			uint64_t const a_scan_min = max((uint64_t)3, strat_a_exp - 1);
			uint64_t const a_scan_max = max((uint64_t)6, strat_a_exp + 10);
			for (uint64_t test_a_exp = a_scan_min; test_a_exp <= a_scan_max; test_a_exp++) {
				uint64_t const b_scan_min = max((uint64_t)2, test_a_exp/2-1);
				uint64_t const b_scan_max = test_a_exp-3;
				for (uint64_t test_b_exp = b_scan_min; test_b_exp <= b_scan_max; test_b_exp++) {
					test_specific_sampling_choice(test_a_exp, test_b_exp);
				}
			}
		}
		//we are in a global optimization cycle, choose now some specific a,b parameters from this cycle to then measure their impact on query runtime
		else if (TREE3_STRATEGY==GET_FASTEST_TREE_PARAMS) {
			auto cycled_params = get_cycled_params(true);
			uint64_t cycle_a_exp = cycled_params.a_star_exp + 0.5; //round to nearest
			uint64_t cycle_b_exp = cycled_params.b_star_exp + 0.5;
			//Reset global knowledge, such that this cycled parameters are kept, even if they consume more space than with the theo parameters
			opt_packed_size= 1ULL<<60;
			test_specific_sampling_choice(cycle_a_exp, cycle_b_exp);
		}
		a_exp = opt_a_exp;
		b_exp = opt_b_exp;
		if constexpr (WRITE_RESULT_FILE) {
			write_ovhd_results_to_file();
		}
	}

	Tree3() {
		reset_tree();
	}

	void build_Tree3(uint64_t my_L_exp) {
		L_exp = my_L_exp;
		// cout << "Start build tree" << endl;
		calc_theoretical_worstcase_overhead();
		find_best_sampling_choices();
		if (TREE3_STRATEGY==GET_SMALLEST_TREE_PARAMS) {
			//in case we scanned multiple (a,b) parameters, now build the final tree with the parameters that led to the smallest tree
			reset_tree();
			build_tree();
			pack_tree();
		}
		save_treestats();
		if (SHOW_T3_VERBOSE_STATISTICS) {
			print_verbose_stats();
			print_ovhd_stats();
		}
	}
};

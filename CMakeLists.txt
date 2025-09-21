cmake_minimum_required(VERSION 3.28)
project(BV)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")

#Compiler paths
set(CLANG_PATH   "/usr/bin/clang"   CACHE STRING "path to clang compiler")
set(CLANGPP_PATH "/usr/bin/clang++" CACHE STRING "path to clang++ compiler")
set(GCC_PATH "/usr/bin/gcc" CACHE STRING "path to gcc compiler")
set(GPP_PATH "/usr/bin/g++" CACHE STRING "path to g++ compiler")
set(USE_CLANG_LIBC  0 CACHE BOOL "link with --stdlib=libc++, for pure clang run on server")

#Most general compile options
set(MYCOMPILER ""   CACHE STRING "specify the compiler: either gcc or clang")
set(NDEBUG      1   CACHE BOOL "Disable debug mode")
set(OPTIMIZE    3   CACHE STRING "Optimization level")
set(SHOW_WARNINGS 0 CACHE BOOL "Show warnings during compilation")
set(FORCE_CRIT_PATH 1 CACHE BOOL "Measure true(er) latencies by XOR-ing each result onto a global dummy, enforcing a critical path between all queries.
                                  Intended to prevent any parallelism that the CPU secretly might cook up by overlapping queries")
set(WARMUP_ROUNDS   0 CACHE STRING "Do a warmup with 5% of the last queries, to remove any lingering (cache) effects from the building process, and simulate a steady state")

if(NOT ${SHOW_WARNINGS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
endif ()

#Option to not hide old (by now unused) code. Not loading compiles faster.
set(LOAD_OLD_CODE 0 CACHE BOOL "load old code, which is generally unused")

#Competitors
set(LOAD_PASTA 1 CACHE BOOL "load pasta and tlx libraries")
set(LOAD_SDSL  1 CACHE BOOL "load SDLS library")
set(LOAD_SUX   1 CACHE BOOL "load sux library")
set(LOAD_POPPY 1 CACHE BOOL "load cs-poppy(-fs) library")

#Parameters for competitors
set(RRR_BLOCKSIZE       31 CACHE STRING "parameter for sdsl_rrr select")
#set(SIMPLE_SELECT_PARAM  2 CACHE STRING "tuning parameter for sux simple_select")
set(CMAKE_SAMPL_EXP 0 CACHE STRING "Sampling exponent, allows to inject custom sampling rate into pasta")



if(${LOAD_SDSL})
    #this "add_subdirectoy" needs to be added BEFORE target definition in "add_executable"
    add_subdirectory(external/sdsl-lite)
endif ()

if(${LOAD_SUX})
    add_library(sux_rank_select INTERFACE)
    target_include_directories(sux_rank_select INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/sux)
endif ()

if(${LOAD_POPPY})
    #original cs-poppy
    add_library(poppy_rank_select ${CMAKE_CURRENT_SOURCE_DIR}/external/rankselect/bitmap.cc)
    target_include_directories(poppy_rank_select PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/external/)

    #newer select version
    add_library(poppy_rank_select_improved ${CMAKE_CURRENT_SOURCE_DIR}/external/rankselect2/bitmap.cc)
    target_include_directories(poppy_rank_select PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/extlib/)
endif()


add_executable(BV bv.cpp)


#Structure of the output file
set(BENCHING_AXIS 1 CACHE STRING "Which parameter should be the default x-axis when plotting the results.
                                  0 = sampledist, 1 = 01ratio, 2=nbits.")


#Parameters for new compressed 3-level Tree
#set(SUMMARY_LEVELS 2 CACHE STRING "the number of summary levels, in {1,2}. Default is 1")
#set(ALPHA     4   CACHE STRING "the alpha-parameter for the 3* tree. If a group spans >= alpha superblocks, it is split into the next-lower level")
set(CUSTOM_L0   0     CACHE STRING "how many bits an L0-block should span in m3. Choices:{512,1024,2048,0}   (0 defaults to 2048)")
set(TREE_TYPE     3   CACHE STRING "select the tree type (3 == 3-level fully compressed delta-3*.  5==uncompressed 5-level original tree")
#set(TREE3_PARAM_CHOICE 6 CACHE STRING "Choose (a,b) params via 1:cubic theory, 2: t<s theory, 3: expected top-level, 4: experimentally scan for size")
#set(BV_COMPRESSION 1 CACHE BOOL "enable bitvector compression via delta-encoding")
set(CHECK_ALL_DENSE_SPARSE  0   CACHE BOOL "use alldense/allsparse shortcuts. Reduces the number of instructions but introduces new branches")
set(COUNT_HIT_TYPES 0 CACHE BOOL "Count the different ways the select queries resolves")
set(COUNT_CPU_CYCLES 0 CACHE BOOL "Count cpu cycles for different section of our query")
#set(Z_SEARCH_STRAT 1 CACHE STRING "How to find the correct z value in a non-matching sparse array. 0: linear scan. 1: binary search")

#Parameters for m2/m5 select within a superblock
set(SCAN        2   CACHE STRING "stepsize for scanning through L1s: m2/m5: {1,2,3,4}")
set(EXTRACT     4   CACHE STRING "method of extracting L2 value. m2:{4,5}, m5:{3,5,6}")
set(POP         2   CACHE STRING "number of popcnts per loop: m2/m5: {1,2}")
set(REVERSE     2   CACHE STRING "choosing if pop-direction can be reversed: 0=no, 1=yes via cmove, 2=yes via branch")
#m2: ex 4, scan 2
#m5: ex 5, scan 1
#m7: ex 5, scan 1


#Further parameters for m2/m3/m5/m7
set(L1MINISKIP  1   CACHE BOOL "Enable L1-Miniskip")
set(ALIGN       1   CACHE BOOL "align data to cacheline boundaries")
set(HUGEPAGES   1   CACHE BOOL "Set 2MB pages")
set(PREFETCH_NT 1   CACHE BOOL "Enable non-temporal prefetching hint")
set(PREFETCH_T0 1   CACHE BOOL "Enable prefetch into L1 and L2 cache")
set(INCL_LAST_L2 1   CACHE BOOL "Reconstruct the last L2 from neighbouring L1's (in m2)")
set(INV_POP_SELECT0 0   CACHE BOOL "Invert raw bits when popcounting in select 0 instead of subtracting (in m2)")
set(PACK_S0_S1  1   CACHE BOOL "pack the entries in the classical S0,S1 sampling")



#Parameters for old uncompressed 5-level Tree
set(PACK_TREE 0   CACHE BOOL "pack entries in the tree to log(n/L1) bits")
set(PACK_EXPL 1   CACHE BOOL "pack specifically the EXPL array")
set(TREE_EXP_LVL0      19  CACHE STRING "highest sampling exponent when using Tree, same role as SAMPLEDIST_EXP")
set(TREE_DIVISION_EXP  3   CACHE STRING "exponent of the number of subblocks per tree-division")
set(TREE_MAX_SUPERBLOCKS      50  CACHE STRING "how many superblocks can be spanned by one tree entry before it is split recursively" )
set(TREE_INV_EXPL_OVHD 1000 CACHE STRING "inverse of maximum explicit-array overhead")
set(BBRANCHLESS 0   CACHE BOOL "branchless B-Tree")

#Reduce terminal output noise
set(SILENT      1   CACHE BOOL "reduce terminal output")
set(EVEN_MORE_SILENT 1   CACHE BOOL "reduce terminal output even more")

#Deprecated or soon to be deprecated things
set(UNROLL      0   CACHE BOOL "Enable loop unrolling")
set(ASM_OUTPUT  0   CACHE BOOL "Output assembly")
#Testing how much pdep slows down on Zen2 architecture
set(AVOID_PDEP  0   CACHE BOOL "replace pdep in m2 and tree with equivalent bit-tricks")

#miscellaneous debugging helpers
set(LOAD_MATPLOTLIB 0 CACHE BOOL "load a python wrapper for matplotlib")
set(CHECK_DELTAW 0 CACHE BOOL "during validation, print the number of words that have been pop-counted")

#Legacy option, for algorithms before the L0/L1/L2 renaming, where L2 was still the lowest block level
set(CUSTOM_L2   0     CACHE STRING "how many bits an L2-block in should span. 0 == using default of each m{2,3,5,7} structure")

if(LOAD_MATPLOTLIB)
    find_package(Python3 COMPONENTS Interpreter Development NumPy)
    include_directories(${Python3_INCLUDE_DIRS})
    include_directories("/usr/lib/python3/dist-packages/numpy/core/include/")
    target_link_libraries(BV ${Python3_LIBRARIES})
endif ()

if(HUGEPAGES)
    set(ALIGNMENT 2097152)# 2^21, 2MiB
    set(ALIGN 1) #hugepages necessarily align
elseif (ALIGN)
    set(ALIGNMENT 64)
else ()
    set(ALIGNMENT 8)
endif ()

message(STATUS "\n")
message(STATUS "loading external libraries:")
if(${LOAD_PASTA})
    message(STATUS "---> pasta")
    target_include_directories(BV PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endif ()

if(${LOAD_SDSL})
    message(STATUS "---> sdsl")
    target_include_directories(BV PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/sdsl-lite/include)
    target_link_libraries(BV PRIVATE sdsl)
endif ()

if(${LOAD_SUX})
    message(STATUS "---> sux")
    target_link_libraries(BV PRIVATE sux_rank_select)
endif ()

if(${LOAD_POPPY})
    message(STATUS "---> poppy")
    target_link_libraries(BV PRIVATE poppy_rank_select)
    target_link_libraries(BV PRIVATE poppy_rank_select_improved)
endif ()
message(STATUS "")


target_compile_definitions(BV PUBLIC
        NDEBUG=$<BOOL:${NDEBUG}>
        OPTIMIZE=${OPTIMIZE}
        FORCE_CRIT_PATH=$<BOOL:${FORCE_CRIT_PATH}>
        WARMUP_ROUNDS=${WARMUP_ROUNDS}

        BENCHING_AXIS=${BENCHING_AXIS}

#        SUMMARY_LEVELS=${SUMMARY_LEVELS}
        TREE_TYPE=${TREE_TYPE}
#        ALPHA=${ALPHA}
        CUSTOM_L0=${CUSTOM_L0}
#        TREE3_PARAM_CHOICE=${TREE3_PARAM_CHOICE}
#        BV_COMPRESSION=$<BOOL:${BV_COMPRESSION}>
#        Z_SEARCH_STRAT=${Z_SEARCH_STRAT}
        COUNT_HIT_TYPES=$<BOOL:${COUNT_HIT_TYPES}>
        CHECK_ALL_DENSE_SPARSE=$<BOOL:${CHECK_ALL_DENSE_SPARSE}>

        SCAN=${SCAN}
        EXTRACT=${EXTRACT}
        REVERSE=${REVERSE}
        POP=${POP}

        L1MINISKIP=$<BOOL:${L1MINISKIP}>
        ALIGN=$<BOOL:${ALIGN}>
        ALIGNMENT=${ALIGNMENT}
        HUGEPAGES=$<BOOL:${HUGEPAGES}>
        PREFETCH_T0=$<BOOL:${PREFETCH_T0}>
        PREFETCH_NT=$<BOOL:${PREFETCH_NT}>
        INCL_LAST_L2=$<BOOL:${INCL_LAST_L2}>
        INV_POP_SELECT0=$<BOOL:${INV_POP_SELECT0}>

        PACK_TREE=$<BOOL:${PACK_TREE}>
        PACK_EXPL=$<BOOL:${PACK_EXPL}>
        TREE_EXP_LVL0=${TREE_EXP_LVL0}
        TREE_DIVISION_EXP=${TREE_DIVISION_EXP}
        TREE_MAX_SUPERBLOCKS=${TREE_MAX_SUPERBLOCKS}
        TREE_INV_EXPL_OVHD=${TREE_INV_EXPL_OVHD}
        BBRANCHLESS=$<BOOL:${BBRANCHLESS}>

        LOAD_OLD_CODE=$<BOOL:${LOAD_OLD_CODE}>
        LOAD_PASTA=$<BOOL:${LOAD_PASTA}>
        LOAD_SDSL=$<BOOL:${LOAD_SDSL}>
        LOAD_SUX=$<BOOL:${LOAD_SUX}>
        LOAD_POPPY=$<BOOL:${LOAD_POPPY}>
        RRR_BLOCKSIZE=${RRR_BLOCKSIZE}
#        SIMPLE_SELECT_PARAM=${SIMPLE_SELECT_PARAM}
        CMAKE_SAMPL_EXP=${CMAKE_SAMPL_EXP}
        CUSTOM_L2=${CUSTOM_L2} #not used by m3, but leave defined for legacy reasons

        JUMP_HEURISTIC=0
        SIMD_RANK=1
        STATIC_IN_FUNC=0
        UNROLL=$<BOOL:${UNROLL}>
        ASM_OUTPUT=$<BOOL:${ASM_OUTPUT}>
        AVOID_PDEP=$<BOOL:${AVOID_PDEP}>
        PACK_S0_S1=$<BOOL:${PACK_S0_S1}>

        SILENT=$<BOOL:${SILENT}>
        EVEN_MORE_SILENT=$<BOOL:${EVEN_MORE_SILENT}>

        LOAD_MATPLOTLIB=${LOAD_MATPLOTLIB}
        CHECK_DELTAW=$<BOOL:${CHECK_DELTAW}>

)

if(MYCOMPILER)
    if(MYCOMPILER MATCHES "gcc")
        set(CMAKE_C_COMPILER ${GCC_PATH})
        set(CMAKE_CXX_COMPILER ${GPP_PATH})
        #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
        message(STATUS "Compiler set to: gcc")

    elseif (MYCOMPILER MATCHES "clang")
        set(CMAKE_C_COMPILER ${CLANG_PATH})
        set(CMAKE_CXX_COMPILER ${CLANGPP_PATH})
        if(USE_CLANG_LIBC)
            #extra links,flags to compile with clang on the server,
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --stdlib=libc++")
            set(CMAKE_EXE_LINKER_FLAGS="-Wl,-rpath,${CLANG_FOLDER}/lib")
        endif()
        message(STATUS "Compiler set to: clang")
    else ()
        message(FATAL_ERROR "Invalid compiler shorthand. Either gcc or clang.")
    endif ()
else ()
    message(STATUS "Compiler: not set by user. Resorting to current IDE default")
endif ()
message(STATUS "${CMAKE_CXX_COMPILER}")


if(OPTIMIZE EQUAL 2)
    target_compile_options(BV PRIVATE -O2)
    message(STATUS "compiling with -O2")
elseif (OPTIMIZE EQUAL 3)
    target_compile_options(BV PRIVATE -O3)
    message(STATUS "compiling with -O3")
else ()
    message(STATUS "compiling without optimization, -O0")
endif ()

if(UNROLL)
    target_compile_options(BV PRIVATE -funroll-loops)
endif()

if(NDEBUG)
    target_compile_options(BV PRIVATE -DNDEBUG)
else ()
    target_compile_options(BV PRIVATE -g)
endif()

if(ASM_OUTPUT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -save-temps")
    message(STATUS "asm output")
endif ()


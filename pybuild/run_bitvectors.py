#!/usr/bin/python3
import copy
from typing import Dict, Union
import os, time, subprocess
from itertools import product

#Reference
my_names    = ["m2","m5","m7","m2-8","m3"]
pasta_names = ["pastaflat"]
sdsl_names  = ["sdsl_sd", "sdsl_rrr", "sdsl_mcl", "sdsl_v1","sdsl_v5"]
poppy_names = ["poppy", "poppy2"]
sux_names   = ["rank9sel", "simple_s"]

comment = ""
folder = ""

#Default Command line parameters
base_params = {
    "-benching_axis": "ratio",
    "-bvname": "m3",        #Algorithm-Selection. m3: the new 3-star structure. 
    "-job":    "multibench", 
    "-bitgen": "indep",     #Instance selection. Indep: independent bits.
    "-qtype":  "random",    #Query type selection.
    "-nbits":    "1",       #Instance size, in units of 10^9. Here small example with 1e9 bits.
    # "-nbits":    "1,8,64",  
    "-selects":  "0",
    "-select0s": "0",
    "-select1s": "2e7",     #Do 20 Mio select_1
    "-ranks":    "0",
    "-iterations": "1",
    "-instances":  "1",
    "-01ratio":"1,9,99", #Instance densities. Here small example 50%, 10%, 1%,
    # "-01ratio":"1,2,4,9,19,32,49,99", #50%, 33%, 20%, 10%, 5%, 3%, 2%, 1%
    "-folder_out":  "../multibench_out/",
    "-parameter_cycles":"18", 
    "-invert_all": "0",
    "-alpha" : "16",
    "-tree3-strat" : "4",#Choice of parameter "a" for 3* structure:  2=a*, 4=a_dense, 5=a_min, 6=a_fast
    "-summary-levels": "1",
    "-z-search-strat" : "1",
    "-bv-compression" : "0", #Experimental bitvector-compression
}


#Command-line parameters (starting with -):  Passed as one long string, values divided by comma
#Compile time parameters (UPPERCASE): Passed as individual entries in a list

sweep_choices1 : Dict[str, Union[str,list]] = {
    "-bvname" : "m3,sdsl_sd,sdsl_mcl,rank9sel,simple_s1,simple_s3",
}

sweep_choices2 : Dict[str, Union[str,list]] = {
    # INDEP
    # "CUSTOM_L0": [2048],
    # "-tree3-strat" : "2",#
    # "-alpha" : "16",
    # "-summary-levels": "1",
    # "-bv-compression" : "0,1",

    # "-bvname" : "m3",
    # "-bitgen": "gap",  # indep, bimodal, chunky, sin
    # "-gap-size" : "1000,10000,100000,1000000,10000000,100000000",
    # "-gap-size" : "300,3000,30000,300000,3000000,30000000",
    # "-01ratio":"1", 
    # "-instances":  "1",
    # "-nbits":  "8",  # 1,2,4,8,16,32,64 (in units of 10^9)
}

sweep_choices3 : Dict[str, Union[str,list]] = {
    #INDEP
    # "CUSTOM_L0": [2048],
    # "-tree3-strat" : "4",#
    # "-alpha" : "8",
    # "-summary-levels": "2",
    # "-bv-compression" : "0,1",

    # "-bvname" : "m3",
    # "-bitgen": "gap",  # indep, bimodal, chunky, sin
    # "-gap-size" : "1000,10000,100000,1000000,10000000,100000000",
    # "-gap-size" : "300,3000,30000,300000,3000000,30000000",
    # "-01ratio":"1", 
    # "-instances":  "1",
    # "-nbits":  "8",  # 1,2,4,8,16,32,64 (in units of 10^9)
}

sweep_choices4 : Dict[str, Union[str,list]] = {
    #INDEP
    # "CUSTOM_L0": [2048],
    # "-tree3-strat" : "5",#
    # "-alpha" : "32",
    # "-summary-levels": "1",
    # "-bv-compression" : "0,1",

    # "-bvname" : "m3",
    # "-bitgen": "gap",  # indep, bimodal, chunky, sin
    # "-gap-size" : "1000,10000,100000,1000000,10000000,100000000",
    # "-gap-size" : "300,3000,30000,300000,3000000,30000000",
    # "-01ratio":"1", 
    # "-instances":  "1",
    # "-nbits":  "8",  # 1,2,4,8,16,32,64 (in units of 10^9)

    # SINUS 20
    # "-bitgen": "sin",  # indep, bimodal, chunky, sin
    # "CUSTOM_L0": [512],
    # "-tree3-strat" : "2,4,5",
    # "-instances":  "20",
    # "-alpha" : "1,2,4,8,16,32,64",
    # "-01ratio":"1,2,4,6,9,13,19,32,49,65,99,199", 
    # "-bvname" : "m3,simple_s0,simple_s1,simple_s2,simple_s3,simple_s4,rank9sel,sdsl_sd,sdsl_mcl,pastaflat",
}


sweep_choices5 : Dict[str, Union[str,list]] = {
    # GAP GADGET

    #RANK

    # "-bvname" : "pastaflat,rank9sel,poppy2,sdsl_v1,sdsl_v5",
    # "-bvname" : "m3",
    # "-bitgen": "gap",  # indep, bimodal, chunky, sin
    # "-gap-size" : "1000,10000,100000,1000000,10000000,100000000",
    # "-01ratio":"1", 
    # "-instances":  "1",
    # "-nbits":  "8",  # 1,2,4,8,16,32,64 (in units of 10^9)
    # "-alpha" : "1,2,4,8,16,32,64",
    # "-tree3-strat" : "2,4,5",
    # "CUSTOM_L0": [512,1024,2048],
}


sweep_choices6 : Dict[str, Union[str,list]] = {
    # GAP GADGET
    # "-bvname" : "simple_s0,simple_s1,simple_s2,simple_s3,rank9sel,sdsl_sd,sdsl_mcl,pastaflat",
    # "-bitgen": "gap",  # indep, bimodal, chunky, sin
    # "-gap-size" : "1000,10000,100000,1000000,10000000,100000000",
    # "-gap-size" : "300,3000,30000,300000,3000000,30000000",
    # "-01ratio":"1", 
    # "-instances":  "1",
    # "-nbits":  "8",  # 1,2,4,8,16,32,64 (in units of 10^9)
    # "-alpha" : "8",
    # "-tree3-strat" : "2",
    # "CUSTOM_L0": [1024],
}


#Compile time cmake parameters
base_config: Dict[str, Union[str,int]] = {
    "MYCOMPILER": "clang",
    "OPTIMIZE" : 3, 
    "BENCHING_AXIS" : 1,  

    "TREE_TYPE": 3,
    "CUSTOM_L0": 2048,
    "CHECK_ALL_DENSE_SPARSE" : 0,
    "COUNT_HIT_TYPES" : 0,

    "ALIGN": 1,
    "HUGEPAGES": 1,
    "PREFETCH_NT": 1,
    "PREFETCH_T0": 1,

    "REVERSE": 2,
    "POP": 2,
    "SCAN": 2,
    "EXTRACT": 4,
    "L1MINISKIP": 1,
    "INCL_LAST_L2": 1,
    "INV_POP_SELECT0": 0,
    "AVOID_PDEP": 0,

    "PACK_TREE" : 0, #removed from m2 in favor of hardcoded 32-bit words
    "PACK_EXPL" : 0,#specifically pack the explicit array
    "TREE_EXP_LVL0" :19, #sampling dist when using tree
    "TREE_DIV" : 3, 
    "TREE_MAX_SUPERBLOCKS" : 50,   #[10,20,40,80,160,320]
    "DEFAULT_TREE_INV_EXPL_OVHD" : 1000, 
    "BBRANCHLESS" : 0,

    "LOAD_PASTA" : 0,
    "LOAD_SDSL" : 0,
    "LOAD_SUX" : 0,
    "LOAD_POPPY" : 0,

    "RRR_BLOCKSIZE" : 31,#sdsl_rrr's block size parameter
    "SAMPL_EXPOS" : 0, #pasta's sampling rate exponent. Leaving at 0 leaves it at default 2^13=8192
    # "DEFAULT_SIMPLE_SELECT_PARAM" : 3, #sux simple_select's tuning parameter

    "WARMUP_ROUNDS" : 0, 
    "FORCE_CRIT_PATH" : 1, 
    "ASM_OUTPUT" : 0,
    "CHECK_DELTAW" : 0,
    "CMAKE_SAMPL_EXP": 0,
    "SHOW_WARNINGS": 0,
    "NDEBUG": 1,
    "SILENT" : 1, 
    "EVEN_MORE_SILENT" : 1,

    "CLANG_PATH": "/usr/bin/clang",
    "CLANGPP_PATH": "/usr/bin/clang++",
    "GCC_PATH": "/usr/bin/gcc",
    "GPP_PATH": "/usr/bin/g++",
    "USE_CLANG_LIBC": 0,
}

sweep_choice_list = [sweep_choices1, sweep_choices2, sweep_choices3,
                     sweep_choices4, sweep_choices5, sweep_choices6]
automatic_compiler_paths = 1
clean_build  = True
show_warnings= False
RUN=1

################## end of user inputs ################################









############## start of automatic script #####################

status_track = {}

def get_cpu_alias():
    alias_from_full= {"":"vb",
                     "feigenbaum":"fg",
                     "hellman"   :"hl",
                     "blum"      :"bl",
                     "yao"       :"yo",
                     "dijkstra"  :"dij",
                     "cook"      :"ck",
                     "karp"      :"kp",
                     "hamming"   :"hg",
                     "hoare"     :"ho",
                     "diffie"    :"dif"}
    fullname = subprocess.run("echo $SLURM_NODELIST",   shell=True, capture_output=True, text=True).stdout.strip()
    alias = alias_from_full.get(fullname)
    if alias is None:
        alias = fullname
    print("fullname",fullname)
    print("alias", alias)
    return alias, fullname


def get_query_alias():
    if base_params["-selects"] != "0": return "s01"
    if base_params["-select1s"] != "0": return "s1"
    if base_params["-select0s"] != "0": return "s0"
    if base_params["-ranks"] != "0": return "r"
    return "mixed_queries"

cpu_alias,fullname   = get_cpu_alias()
query_alias = get_query_alias() 

if folder!="":
    if folder=="+":
        folder = cpu_alias + "-" + query_alias + "e" + base_params["-nbits"]+"-"+base_params["-bitgen"]
    path = "../multibench_out/"+folder+"/"
    base_params["-folder_out"] = path
    try:
        os.mkdir(path)
    except:
        print("Folder already there")

def cmake_flags_from_dict(config: dict) -> list[str]:
    return [f"-D{key}={value}" for key, value in config.items()]


def run_config(config: Dict[str, Union[str,int]], params : Dict[str, str], tag, progress, sweep_overview : Dict[str, Union[str,list]] = {}):
    params["-progress"] = str(progress)
    print("received")
    # print(params)
    # print(config)
    bvnames = params["-bvname"]
    print("bvnames", bvnames)
    bvnames_set = set(bvnames.split(","))
    print("bvnames set:", bvnames_set)
    if len(bvnames_set.intersection(pasta_names))>0 : config["LOAD_PASTA"]=1
    if len(bvnames_set.intersection(sdsl_names))>0:  config["LOAD_SDSL"]=1
    if len(bvnames_set.intersection(poppy_names))>0:  config["LOAD_POPPY"]=1
    if len(bvnames_set.intersection(sux_names))>0:   config["LOAD_SUX"]=1

    if automatic_compiler_paths:
        config["CLANG_PATH"]   = subprocess.run("which clang",   shell=True, capture_output=True, text=True).stdout.strip()
        config["CLANGPP_PATH"] = subprocess.run("which clang++", shell=True, capture_output=True, text=True).stdout.strip()
        config["GCC_PATH"] = subprocess.run("which gcc", shell=True, capture_output=True, text=True).stdout.strip()
        config["GPP_PATH"] = subprocess.run("which g++", shell=True, capture_output=True, text=True).stdout.strip()

    cpu_alias,fullname   = get_cpu_alias()
    if(fullname!=""): #i.e. we are on the server
        config["USE_CLANG_LIBC"]=1 #for the server

    # bv_alias    = bv_alias_from_full[bvname]
    # if bv_alias != "m3": bv_alias="others"
    bv_alias="m3"
    if(bvnames!="m3"):
        bv_alias="others"

    query_alias = get_query_alias() 
    bitgen_alias = params["-bitgen"]
    filename = cpu_alias+"_"+bv_alias+"_"+query_alias+"_n"+params["-nbits"]+"_"+bitgen_alias
    if comment != "":
        filename +="-"+comment
    if folder=="":
        filename = tag+"_"+filename
    params["-filename_out"] = filename
    
    if status_track['sweep_overview_written_to_file']==False:
        write_sweep_overview(filename, sweep_overview)
        status_track['sweep_overview_written_to_file']=True

    nproc = int(subprocess.run("nproc", shell=True, capture_output=True, text=True).stdout.strip())
    useproc = int(nproc/2)

    cmake_flags = cmake_flags_from_dict(config)
    cmake_cmd = "cmake " + " ".join(cmake_flags) + " .."
    build_cmd = f"cmake --build . --clean-first -- -j{useproc}"

    #Show each cmake_flag with aligned values
    max_key_length = max(len(flag.split('=')[0]) for flag in cmake_flags)
    for flag in cmake_flags:
        key, value = flag.split('=')
        offs = int((max_key_length - len(key))/2)
        offs = offs + offs%2
        print(f"{key.ljust(max_key_length - offs)} {value}")
    
    # if(ASM_OUTPUT):
        # llvm = "llvm-mca -timeline -bottleneck-analysis -output-asm-variant=1 bv.s > llvm.txt"
        # os.system(llvm)

    os.system(cmake_cmd)
    os.system(build_cmd)

    run_cl_command="./BV"
    for flag,param in params.items():
         run_cl_command+=" "+flag+" "+param   

    print(run_cl_command)
    if RUN:
        os.system(run_cl_command)


def write_sweep_overview(filename, sweep_overview):
    folder = base_params['-folder_out']
    path = folder + filename + ".txt"
    pad=18
    with open(path, "a") as f:
        f.write("Sweep \n")
        for key, val in sweep_overview.items():
            f.write(key.ljust(pad)+"="+str(val)+"\n")
        f.write("\nBase Command-Line Params\n")
        for key, val in base_params.items():
            f.write(key.ljust(pad)+"="+str(val)+"\n")
        f.write("\nBase Compile Config \n")
        for key, val in base_config.items():
            f.write(key.ljust(pad)+"="+str(val)+"\n")
        f.write("\n \n")


def run_sweep():
    tag = str(int(time.time()))
    # for bvname in bvname_list:
    #     params = cl_params_dict.copy()
    #     params["-bvtype"] = bvname
    #     if bvname in my_names:
    for sweep_choices in sweep_choice_list:
        status_track['sweep_overview_written_to_file'] =  False
        for key, value in sweep_choices.items(): 
            if not isinstance(value, list): 
                sweep_choices[key]=[value]

        keys = list(sweep_choices.keys())
        value_lists = list(sweep_choices.values())
        if len(keys)==0:
            continue
        print("Sweeping", keys)
        print("Sweeping", value_lists)
        # Create all combinations
        combinations = [dict(zip(keys, values)) for values in product(*value_lists)]
        for i,combo in enumerate(combinations):
            config = base_config.copy()
            params = base_params.copy()
            for key,value in combo.items():
                print("Now updating", key, value)
                if(key[0]=="-"):
                    params.update({key : value})
                else:
                    config.update({key : value})
            run_config(config, params, tag, i, sweep_choices)

run_sweep()


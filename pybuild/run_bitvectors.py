#!/usr/bin/python3
import copy
from typing import Dict, Union
import os, time, subprocess
from itertools import product


my_names    = ["m3"]
pasta_names = ["pastaflat"]
sdsl_names  = ["sdsl_sd", "sdsl_rrr", "sdsl_mcl", "sdsl_v1","sdsl_v5"]
poppy_names = ["poppy", "poppy2"]
sux_names   = ["rank9sel", "simple_s0", "simple_s1", "simple_s2", "simple_s3"]

base_params = {
    "-benching_axis": "ratio",
    "-bvname": "m3",
    "-job":    "multibench",
    "-bitgen": "indep",      #Choose independen bit instance (alternative is "sin", for sinusoidally based interference instances)
    "-qtype":  "random",     #Do Random queries 
    "-nbits":    "1",        #Chose number of bits in units 10^9
    "-selects":  "0",
    "-select0s": "0",
    "-select1s": "2e7",
    "-ranks":    "0",
    "-iterations": "1",      #Choose number of iterations, for statistical security. For quick test run set to 1 
    "-instances":  "1",
    "-01ratio":"1,4,49", #Choose ratio of Zeros to Ones. A ratio 3 corresponds to a density of 25%. 
    "-folder_out":  "../multibench_out/",
    "-invert_all": "0",
    "-alpha" : "16",            # Choose The 3* Alpha parameter
    "-tree3-strat" : "4",       # Choose strategy for calculating a and b. 2=a*, 4=a_fast, 5=a_min
    "-summary-levels": "1",     # Choose a value in [1,2], corresponds to k=2 and k=3 summary levels
    "-bv-compression" : "0",    # Choose whether to activate the preliminiary bitvector compression
}


#Command-line options: Single string, with commata between entries
#Compile-time options: List of values

sweep_choices1 : Dict[str, Union[str,list]] = {
    # Also Test some competitors on the same instance
    "-bvname" : "m3,rank9sel,simple_s1,pastaflat,sdsl_mcl",
    # "-01ratio":"1,9,99", 
}

sweep_choices2 : Dict[str, Union[str,list]] = {
    "-bvname" : "m3",
    "CUSTOM_L0": [512,2048],   #Test L0=512 and L0=2048
    "-tree3-strat" : "4", #Test with a = a_fast strategy
    "-alpha" : "8",  # Test alpha=8 
    # "-alpha" : "2,8,32",  # Test different alpha
    # "-summary-levels": "1,2",
    # "-bv-compression" : "0,1",
    # "-01ratio":"1", 
    # "-instances":  "1",
    # "-nbits":  "1,2,4,8",  
}



#Compile time cmake parameters. Typicall not needed to change them
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

    #Legacy for old uncompressed tree
    "PACK_TREE" : 0, 
    "PACK_EXPL" : 0,
    "TREE_EXP_LVL0" :19, 
    "TREE_DIV" : 3, 
    "TREE_MAX_SUPERBLOCKS" : 50,   
    "DEFAULT_TREE_INV_EXPL_OVHD" : 1000, 
    "BBRANCHLESS" : 0,

    "LOAD_PASTA" : 0,
    "LOAD_SDSL" : 0,
    "LOAD_SUX" : 0,
    "LOAD_POPPY" : 0,

    "RRR_BLOCKSIZE" : 31,#sdsl_rrr's block size parameter
    "SAMPL_EXPOS" : 0, #pasta's sampling rate exponent. Leaving at 0 leaves it at default 2^13=8192

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

sweep_choice_list = [sweep_choices1, sweep_choices2]
automatic_compiler_paths = 1
clean_build  = True
show_warnings= False
RUN=1

################## end of user inputs ################################









############## start of automatic script #####################
folder =""
status_track = {}


def get_query_alias():
    if base_params["-selects"] != "0": return "s01"
    if base_params["-select1s"] != "0": return "s1"
    if base_params["-select0s"] != "0": return "s0"
    if base_params["-ranks"] != "0": return "r"
    return "mixed_queries"

        
cpu_alias ="mypc"
fullname=""
query_alias = get_query_alias() 

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


repos = {
    "sdsl-lite": "https://github.com/simongog/sdsl-lite.git",
    "sux": "https://github.com/vigna/sux.git",
    "rankselect": "https://github.com/efficient/rankselect.git",
    "rankselect2": "https://github.com/pasta-toolbox-forks/rankselect2.git"
}

def clone_repos(target_dir):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    for folder_name, repo_url in repos.items():
        repo_path = os.path.join(target_dir, folder_name)
        if os.path.exists(repo_path):
            print(f"[OK] Repository '{folder_name}' already exists at {repo_path}")
        else:
            print(f"[CLONE] Cloning {repo_url} into {repo_path}")
            subprocess.run(["git", "clone", repo_url, repo_path], check=True)

def run_sweep():
    tag = str(int(time.time()))
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

clone_repos("../external/")
run_sweep()


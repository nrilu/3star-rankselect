# 3star-rankselect
This is the accompanying source-code to **Theory Meets Practice for Bit Vectors Supporting Rank and Select** (arXiv link)

The new 3* structure is implemented in ```m3.hpp```  (Query and Summary-Tree construction) and ```Tree3.hpp``` (Sample-Tree construction).

# Run
To start a simple select_1 benchmark: run 
```./run_bitvectors.py``` in ```pybuild/```.
Two example settings are provided.


Live results are shown after starting ```./run_bitvectors.py```. In particular, the column "query" contains the average query time and the column "overhead" the space overhead of all introduced structures.
 
Logged results are found in the folder ```multibench_out/```. Runtime is stored in the field "t(ns)", and overhead in "overhead_percent" (among many other logged statistics).



  

# 3star-rankselect
This is the accompanying source-code to **Theory Meets Practice for Bit Vectors Supporting Rank and Select** (arXiv link)

Our new 3* structure is implemented in ```m3.hpp``` and ```Tree3.hpp```. 

# Run
To start an exemplary benchmark of *select_1* queries, run 

```./run_bitvectors.py``` within the ```pybuild/``` folder.

Two example settings are provided, the full example run should take around a minute.


In the live results, the query time and space overhead are directly shown. The results are also logged to ```multibench_out/```. Query time is stored in the field **t(ns)**, and overhead in **overhead_percent**.

At first run, the script auto-downloads the competitor codes into the folder ```external/```.

  

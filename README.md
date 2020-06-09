# NetExtractor

![GitHub Logo](workflow1.png)

## USAGE:

## 1) Pre-processing

Pre-processing includes normalizing of the input Gene expression matrix (GEM). 

Details can be found within the NetExtractor paper, and the normalizing code can be accessed at: https://github.com/SystemsGenetics/GEMprep


## 2) NetExtractor Algorithm [Modules A-D]
```
python NetExtractor.py
```

This code reads in the GEM (GTEx_v7_brain_subGEM-log-no.txt) and outputs a file with with GeneA_name, GeneB_name, MI value, Inter-cluster score value.

The code is for multiprocess and runs on 20 threads. Modify based on resources available.

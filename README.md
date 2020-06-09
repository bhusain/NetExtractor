# NetExtractor

![GitHub Logo](NetEx_workflow.png)

## USAGE:

## 1) Pre-processing

Pre-processing includes normalizing of the input Gene expression matrix (GEM). 

Details can be found within the NetExtractor paper, and the normalizing code can be accessed at: https://github.com/SystemsGenetics/GEMprep


## 2) NetExtractor Algorithm [Modules A-D]
```
python NetExtractor.py
```

This code reads in the GEM (GTEx_v7_brain_subGEM-log-no.txt) and outputs a file (NetExtractor_output.txt) which contains GeneA_name, GeneB_name, MI value, Inter-cluster score value.

The code is parallelized for multiprocesses and runs on 20 threads. Modify based on resources available.

### Important Notes:

1) NetExtractor_output.txt will include unfiltered MI and S values. Thresholds are determined in postprocessing.


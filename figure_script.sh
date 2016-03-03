#!/bin/bash
for ((i=16; i<33; i++)); do echo ========; echo k = $i; quartz/ref_distance_variable $i hg19.fa 2>&1 | tee -a kmer_analysis.log ; done

#! /usr/bin/bash

## merge together the output of the summary results clumping - only select clumps where the lead SNP <= 1e-8
## merge across all chromosomes for each trait
## 1) file stub to search for, i.e. Noise-<protein>:<cell_type>

# this is the directory containing the .clumped files from plink
CLUMP_DIR="GCTA_clump"

# change this to the output directory for the relevant file
MERGE_DIR="GCTA_merged"

clump_files=$(eval "find $CLUMPDIR/. ! -name '.' -prune -name '*clumped' | grep '$1'")

for cfile in $clump_files;
do
    # subset to association with p< 1E-07 - this will be more stringently filtered later
    cat $cfile | awk '$5 <= 1e-7 {print $0}' | sed '/^$/D' >> $MERGE_DIR/$1.merged
done

#! /usr/bin/bash

# for each results file, run the plink report clumping

# the location of plink binary bed files
PLINK_FILES="plink_files"

# the path to a list of SNPs to exclude from analysis - one SNP per line
EXCLUDE_FILES="plink_files/Exclusion_snps.txt"

# results directory containing output files from LIMIX QTL-mapping - with .txt suffix
QTL_DIR="qtl_results"

# scripts directory
SRCDIR="src"

for result in `find $QTL_DIR/. ! -name . -prune -name '*txt'`;
do
    #echo $result
    bash $SRCDIR/submit_plink_clump.sh $PLINK_FILES $result $EXCLUDE_FILES
done

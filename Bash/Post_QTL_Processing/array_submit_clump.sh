#! /usr/bin/bash

# Change the paths below to the location of binary plink files, list of SNPs to exclude and the 
# directory containing the output files from GCTA, i.e. .mlma files, as well as the location of the 
# scripts directory on your machine

# for each results file, run the plink report clumping
PLINK_FILES="plink_files/"
EXCLUDE_FILES="plink_files/Exclusion_snps.txt"
ASSOC_DIR="GCTA_assoc"
SRCDIR="src"

for result in `find $ASSOC_DIR/. ! -name . -prune -name '*mlma'`;
do
    #echo $result
    bash $SRCDIR/submit_plink_clump.sh $PLINK_FILES $result $EXCLUDE_FILES
done

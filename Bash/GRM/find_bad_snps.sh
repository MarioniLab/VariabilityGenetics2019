#! /usr/bin/bash

## Find duplicate, overlapping and tri allelic SNPs on each chromosome
## 1) output file prefix
## 2) input .bim file

# change this to the location of the CGAT geno2geno scripts
CGAT_SRCDIR="CGAT/scripts"


python $CGAT_SRCDIR/geno2geno.py  --task=detect_duplicates  --outfile-pattern=$1  $2

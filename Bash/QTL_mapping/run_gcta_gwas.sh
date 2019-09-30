#! /usr/bin/bash

# 1) input plink files
# 2) input GRM file
# 3) input phenotype file
# 4) output file name
# 5) covariates file
# 6) bad SNPs to exclude, i.e. tri-allelic, duplicate, etc

gcta64 --mlma  --bfile $1  --grm $2  --pheno $3  --out $4 --qcovar $5  --maf 0.02  --threads 10 --exclude $6

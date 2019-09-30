#! /usr/bin/bash

## loop over traits, extract 1Mb cis window and run cis-QTL mapping 
## with LIMIX

# need a table of protein, chr, TSS positions - this can be generated manually or by using an online database such as Ensembl
TSS_TABLE="TSS_table.txt"

# location of file that contains the filepath to each plink file - one per line
PLINK_FILES="plink_files/autosome_plink_files.txt"

# directory containing all Chr and Not-Chr GRM files 
GRM_FILES="GCTA_out/"

# files of covariates to adjust for in the LMM - filenames are the same as triat files with .covar suffix
COVAR_FILES="covar_files"

# file of SNPs to exclude - 1 SNP per line
EXCLUDE_FILES="plink_files/Exclusion_snps.txt"

# directory containing trait files group in to sub-directories for the chromsome on which the protein is encoded
# phenotype files have a .pheno suffix
TRAIT_DIR="trait_files"

# scripts directory
SRCDIR="src"

for pheno in `find $TRAIT_DIR/chr*/*pheno`;
do
 bash $SRCDIR/submit_limix_qtl.sh $PLINK_FILES $GRM_FILES $pheno $COVAR_FILES $EXCLUDE_FILES $TSS_TABLE

done

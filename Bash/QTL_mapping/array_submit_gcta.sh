#! /usr/bin/bash

# change the following variables to the locations of the appropriate files
# for each trait file run a full GWAS using each chromosome

# this is a text file that contains the path to each plink binary bed file - one entry per line
PLINK_FILES="/nfs/research1/marioni/mdmorgan/Noise_genetics/Milieu_Interior/plink_files/autosome_plink_files.txt"

# the directory containing the GRM files
GRM_FILES="/nfs/research1/marioni/mdmorgan/Noise_genetics/Milieu_Interior/GCTA_out/"

# the directory for all covariates files - in this case with same name as the trait to test with the .covar suffix
COVAR_FILES="/nfs/research1/marioni/mdmorgan/Noise_genetics/Milieu_Interior/FCGR2A_covar_files/"

# a text file containing one SNP per line to remove from analyses
EXCLUDE_FILES="/nfs/research1/marioni/mdmorgan/Noise_genetics/Milieu_Interior/plink_files/Exclusion_snps.txt"

# the directory containing the trait phenotype files across individuals - sorted in to sub-directories based on
# which chromosome the protein is encoded on - used to find the correct GRM for the linear mixed model
TRAIT_DIR="trait_files"

# directory containing scripts
SRCDIR="src"

for pheno in `find $TRAIT_DIR/chr*/ -name '*pheno' | grep Mean`;
do
    #echo $pheno
    bash $SRCDIR/submit_gcta_gwas.sh $PLINK_FILES $GRM_FILES $pheno $COVAR_FILES $EXCLUDE_FILES
done

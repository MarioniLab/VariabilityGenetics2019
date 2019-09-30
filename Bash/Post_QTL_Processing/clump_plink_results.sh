#! /usr/bin/bash

## Use LD clumping on the association results to extract the roughly independent regions
## 3) input results file - per chromosome
## 2) SNP exclusions file
## 1) input binary bed file prefix
## 4) output file stub

plink --bfile $1 --exclude $2  --clump $3 --clump-field p --out $4 --memory 5000


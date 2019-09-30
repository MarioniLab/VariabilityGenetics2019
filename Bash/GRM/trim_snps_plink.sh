#! /usr/bin/bash

# This script assumes that plink v1.9 is in your PATH variable - if not then change the front call to plink to the 
# PATH to the location of the plink v1.9 binary file

## Extract a set of ~independent SNPs on each chromsome based on the LD between them
## Remove SNPs in a 5Mb window, shifting 100kb at a time, removing SNP pairs with an r^2 <= 0.01
## 1) input plink file prefix
## 2) output file stub

plink --bfile $1 --extract $2 --make-bed  --exclude $4 --out $3

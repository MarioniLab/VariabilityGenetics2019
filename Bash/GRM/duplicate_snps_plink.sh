#! /usr/bin/bash

## 1) input plink file
## 2) output file

plink --bfile $1  --list-duplicate-vars  --out $2

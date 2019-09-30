#! /usr/bin/bash

# files need to be block gzipped to make use of tabix

# 1) input file name

cat $1 | bgzip > $1.bgz
tabix -f -b 3 -e 3 -S 1 $1.bgz

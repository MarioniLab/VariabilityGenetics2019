#! /usr/bin/bash

# for each results file tabix index the file based on chromosome and variant position
# the tabix indexed files are automatically written to the directory containing the original results files

# change this to the location of your plink binary files
PLINK_FILES="plink_files/"

# this should be filepath for the association results files from GCTA, i.e. .mlma files
ASSOC_DIR="GCTA_assoc"

# change this to the directory containing your scripts
SRCDIR="src"

for result in `find $ASSOC_DIR/ -name '*mlma'`;
do
    JOB="bsub -R "rusage[mem=5000]" -n 1 -T 5 -M 10000 -J $result -q research-rh7 -e $result-tabix.err -o $result-tabix.out $SRCDIR/make_tabix.sh $result"
    #echo $JOB
    eval $JOB

done

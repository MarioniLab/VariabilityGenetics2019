#! /usr/bin/bash


## set the two variables below to the locations of the plink binary bed files, 
## output file location and directory containing the appropriate scripts
## Change the job submission command to the appropriate one for your cluster queue manager

# loop over plink chromosome files and create a GRM per chromosome
# used to estimate h2 with GREML, partitioned by chromosome

PLINK_FILES="plink.dir"
OUTDIR="GCTA_out"
SRCDIR="src"

for x in `find $PLINK_FILES/ -name '*bed'`;
do
    # this assumes a particular file name strcture to extract the total name and chromosome ID
    PLINK=$(eval "echo $x | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev")
    
    CHROME=$(eval "echo $PLINK | cut -d '_' -f 1")
    #oecho $PLINK
    #echo $CHROME
    JOB="bsub -R "rusage[mem=8000]" -n 1 -M 12000 -T 15 -q research-rh74 -J $CHROME.grm -e $OUTDIR/$CHROME.grm.err -o $OUTDIR/$CHROME.grm.out  $SRCDIR/make_grm.sh $PLINK_FILES/$PLINK $OUTDIR/$CHROME"

    eval $JOB
    #echo $JOB
done

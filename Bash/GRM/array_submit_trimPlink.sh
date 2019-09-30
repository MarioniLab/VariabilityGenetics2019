#! /usr/bin/bash

######################################################################################################
## Provide the input directory containing the binary plink files, the intended output directory and 
## directory containing bash scripts  as positional arguments to this script
## It is assumed that you use LSF as your queue manager - therefore change the JOB variable to the appropriate
## script submission commands for your cluster
######################################################################################################                                                                               

## Trim SNPs in LD on each chromosome separately
INDIR=$1
OUTDIR=$2
SRCDIR=$3

for plink in `find $INDIR -name '*bed' | rev | cut -d '.' -f 2- | rev`;
do
    # this assumes a particular structure of the file names to extract the chromosome from the plink bed file name - change according to your filenaming
    CHROME=$(eval "echo $plink | cut -d '/' -f 11 | cut -d '-' -f 1")
    JOBNAME=$CHROME-LDindependent
    #echo $CHROME
    #echo $JOBNAME
    JOB="bsub -R "rusage[mem=5000]" -n 10 -T 16 -M 12000 -J $JOBNAME -q research-rh74 -e $OUTDIR/$JOBNAME.err -o $OUTDIR/$JOBNAME.out   $SRCDIR/trim_snps_plink.sh $plink  $OUTDIR/keep_independent.snps   $OUTDIR/$CHROME-independent"

    eval $JOB
    #echo $JOB

done

#! /usr/bin/bash

## submit a bsub job
## 1) input plink files
## 2) input results file
## 3) SNP exclusions file

# change this to the location of your scripts
SRCDIR="src"

# change this to the output directory on your filesystem
OUTDIR="GCTA_clump"

# These commands assume a particular filename and directory structure - change to suit your own data and filesystem
## remove PATH and .mlma
infile_stub=$(eval "echo $2 | cut -d '/' -f 10 |  rev | cut -d '.' -f 2- | rev ")
summary=$(eval "echo $infile_stub | cut -d '-' -f 1 | cut -d '_' -f1")
#echo $infile_stub
trait=$(eval "echo $infile_stub | rev | cut -d '-' -f 2 | rev | tr -s ':' '_' ")
#echo $trait

JOBNAME=$(eval "echo $infile_stub")
ERROR_FILE=$JOBNAME.err
CLUSTER_OUT=$JOBNAME.out

# loop over the plink chromosome files
# append chromosome name
chrome=$(eval "echo $2 | cut -d '/' -f 10 | rev | cut -d '.' -f 2- | rev | cut -d '-' -f3")

# change this to match the structure of your filenames
plink=$(eval "find $1$chrome-LabExMI* | grep bed | rev | cut -d '.' -f 2- | rev")


# change the job submission commands to those for your queue manager
JOB="bsub -R "rusage[mem=5000]" -n 1 -T 16 -M 12000 -J $JOBNAME -q research-rh74 -e $OUTDIR/$ERROR_FILE -o $OUTDIR/$CLUSTER_OUT   $SRCDIR/clump_plink_results.sh $plink $3 $2 $OUTDIR/$infile_stub"

#echo $chrome
#echo $JOB
eval $JOB

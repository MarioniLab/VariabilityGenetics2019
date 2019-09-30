#! /usr/bin/bash

## submit LIMIX qtl job
## 1) input plink files
## 2) input GRM file
## 3) input phenotype file
## 4) covariates directory
## 5) SNP exclusion file
## 6) table of protein/gene TSS


# these commands are based on a particular filename and directory structure - adjust according to your own filesystem
infile_stub=$(eval "echo $3 | cut -d '/' -f 10 | rev | cut -d '.' -f 2- | rev")
summary=$(eval "echo $infile_stub | cut -d '-' -f 1 | cut -d '_' -f1")

trait=$(eval "echo $infile_stub | rev | cut -d '-' -f 1 | rev | tr -s ':' '_' ")

JOBNAME=$(eval "echo $trait")
ERROR_FILE=$JOBNAME.err
CLUSTER_OUT=$JOBNAME.out

COVAR_FILE=$(eval "echo $infile_stub")

# find the TSS matching the protein trait
protein=$(echo $trait | cut -d '_' -f 1 | awk '{printf("%s\t", $1)}')

TSS=$(eval 'cat $6 | grep "$protein" | cut -f 3')
START=$((TSS - 500000))
END=$((TSS + 500000))

# find the correct GRM, i.e. Not-chrN
CHROME=$(eval 'cat $6 | grep "$protein" | cut -f 2 | sed "s/chr/Chr/g"')
#echo $CHROME

GRM=$(eval "find $2/ -name '*grm.bin' | grep Not-$CHROME.grm")
#GRM=$(echo $GRM | rev | cut -d '.' -f 3- | rev)
#echo $GRM

# output file directoy
OUTDIR="qtl_results"

# directory to re-route script logging to
LOGDIR="logs"

for plink in `cat $1 | grep "$CHROME-"`;
do
    chrome=$(eval "echo $plink | cut -d '/' -f 9 | cut -d '-' -f 1")

    # change the following job submissions command for your relevant queue manager
    JOB="bsub -R "rusage[mem=22000]" -n 10 -T 150 -M 28000 -J $JOBNAME -q research-rh74 -e $OUTDIR/$ERROR_FILE -o $OUTDIR/$CLUSTER_OUT python3 $SRCDIR/cis_limix.py --genotypes $plink --phenotype $3 --covariates $4/$infile_stub.covar --grm $GRM --maf 0.05 --start $START --end $END  --log $LOGDIR/$JOBNAME.log --output $OUTDIR/$summary-$trait-$chrome-cisqtl.txt --permutation --iters 1000"
 
    echo $JOB
    eval $JOB
done

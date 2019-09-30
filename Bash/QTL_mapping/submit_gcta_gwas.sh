#! /usr/bin/bash

# submit a bsub job
# 1) input plink files
# 2) input GRM file
# 3) input phenotype file
# 4) covariates directory
# derive output file name from input pheno file
# derive jobname, output and error files from input pheno file

# these commands assume a particular directory and filename structure - please change these according to your own filesystem and filenames
# remove PATH and .pheno
infile_stub=$(eval "echo $3 | cut -d '/' -f 10 |  rev | cut -d '.' -f 2- | rev ")
summary=$(eval "echo $infile_stub | cut -d '-' -f 1 | cut -d '_' -f1")

trait=$(eval "echo $infile_stub | rev | cut -d '-' -f 1 | rev | tr -s ':' '_' ")

JOBNAME=$(eval "echo $infile_stub")
ERROR_FILE=$JOBNAME.err
CLUSTER_OUT=$JOBNAME.out

COVAR_FILE=$(eval "echo $infile_stub")

# the location of scripts
SRCDIR="src"

# the location of the output files
OUTDIR="GCTA_assoc"


# loop over the plink chromosome files
for plink in `cat $1`;
do
    # append chromosome name
    #echo $plink
    chrome=$(eval "echo $plink | cut -d '/' -f 9 | cut -d '-' -f 1")
    #echo $chrome
    grm=$(eval 'find $2 -name "*grm.bin" | grep  "${chrome}.grm.bin" | grep "Not" | rev | cut -d "." -f 3- | rev')
    #echo $grm

    # change the following job submission to use the appropriate commands for your queue manager
    JOB="bsub -R "rusage[mem=5000]" -n 10 -T 16 -M 12000 -J $JOBNAME -q research-rh74 -e $OUTDIR/$ERROR_FILE -o $OUTDIR/$CLUSTER_OUT $SRCDIR/run_gcta_gwas.sh $plink $grm $3 $OUTDIR/$summary-$trait-$chrome $4/$infile_stub.covar $5"

    echo $JOB
    eval $JOB
done

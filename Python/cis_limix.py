#! /usr/bin/python3

import limix
import pandas_plink
import os
import sys
import numpy as np
import pandas as pd
import argparse
import logging
import xarray
import scipy.stats as stats
import functools

def get_min_p(p_series, min_series):
    '''
    Given a set of minimum p-values, assess whether any of the new p-values 
    are smaller. If they are then store these and return, otherwise retain
    the original input p-values

    p_series: pandas.Series

    min_series: pandas.Series
    '''

    # compare p-values and store if new p-values are smaller than old
    lt_p = p_series < min_series
    #print(p_series)
    #print(min_series)
    #print(lt_p)
    min_series.loc[lt_p] = p_series.loc[lt_p]

    return(min_series)    


def perm_qtl(genos, phenos, covars, grm, bim, method="beta", n_iters=100):
    '''
    Perform gene-wise permutation-based multiple testing procedure. The null distribution is either
    generated empirically from a fixed number of iterations, or approximated using a beta-distribution 
    as per FastQTL. Only genotypes are permuted, all other values are fixed, i.e. phenotypes, covariates, GRM

    genos: numpy.ndarray of genotypes (0, 1, 2, -9)

    phenos: numpy.ndarray of phenotyps

    grm: xarray of realized genetic relationships

    bim: SNP position and ID information, i.e. Plink BIM file

    covars: pandas.DataFrame of covariates to adjust in model

    method: either beta or fixed

    n_iters: int of number of iterations to perform for either beta-approximation or fixed permutations.
    '''
    BETA_SHAPE1_MIN = 0.05
    BETA_SHAPE1_MAX = 10
    BETA_SHAPE2_MIN = 1
    BETA_SHAPE2_MAX = 1000000

    # need to approximate beta-distribution parameters using permutations to then calculate p-values per gene
    # first step is to randomly permute the genotypes
    # merge dataframe on the fly

    null_p = []
    for i in range(n_iters):
        logging.info("Performing iteration {} of {}".format(i, n_iters))
        perm_phenos = np.random.permutation(phenos)
        perm_lmm = limix.qtl.scan(G=genos.T, Y=perm_phenos,
                                  lik="normal", K=grm, M=covars,
                                  verbose=False)

        res_df = perm_lmm.stats
        res_df = pd.merge(res_df, bim, left_index=True, right_index=True)
        # extract the SNP effect sizes
        coef_df = perm_lmm.effsizes['h2']

        cov_eff_df = coef_df.loc[coef_df.effect_name.isin(covars.columns[2:]), :]
        coef_df = coef_df.loc[~coef_df.effect_name.isin(cov_names), :]
        coef_df.set_index('test', drop=False, inplace=True)
        res_df = pd.merge(res_df, coef_df, left_index=True, right_index=True)
        res_df = res_df.loc[:, ["snp", "pv20"]]

        if i == 0:
            min_p = res_df.loc[:, "pv20"]
        else:
            min_p = get_min_p(p_series=res_df.loc[:, "pv20"],
                              min_series=min_p)

        # merging within the permutation loop makes it really slow, and doesn't
        # necessarily get around the memory issues
        # I only actually need to store the smallest null p-value for each SNP
        # rather than all of them

    top_null_p = min_p
    #print(top_null_p)

    # MLE of a(k), b (n) parameters of the beta distribution
    # for each SNP
    try:
        k, n, floc, fscale = stats.beta.fit(top_null_p, floc=0, fscale=1)
    except stats._continuous_distns.FitSolverError:
        logging.warning("MLE convergence failed - trying restart")
        k, n, floc, fscale = stats.beta.fit(top_null_p, floc=0, fscale=1)

    # need to check convergence of MLEs
    # this can be done by taking the partial dreivatites with respect to the shape parameter is negative,
    # or that the variances of the log-transformed variables are positive, i.e. > 0
    # equivalently I can check that the shape parameters are within expected limits as per Marc's implementation
    
    if(k < BETA_SHAPE1_MIN or k > BETA_SHAPE1_MAX or n < BETA_SHAPE2_MIN or n > BETA_SHAPE2_MAX):
        logging.warning("Beta-parameter estimates are outside expected bounds. This is most likely caused by "
                        "a lack of convergence of MLEs. Setting P-value to NaN")
        k = np.nan
        n = np.nan
    else: pass

    # calculate a beta-distribution based on these parameter estimates
    bdist = stats.beta(k, n)
    correction_function = lambda X: bdist.cdf(X)
    
    return(correction_function)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="A script for performing cis-eQTL association mapping using linear mixed models")
    parser.add_argument("--genotypes", dest="geno_path", type=str,
                        help="The file prefix to the plink-format genotypes in fam, bim and bed format")

    parser.add_argument("--phenotype", dest="pheno_file", type=str,
                        help="The file containing a vector of phenotypes and sample IDs. Missing values should be NA")

    parser.add_argument("--covariates", dest="covar_file", type=str,
                        help="A table of covariates to include in the LMM, with FID and IID columns to identify individuals")

    parser.add_argument("--grm", dest="grm_file", type=str,
                        help="A genetic relatedness matrix constructed from SNPs that are not on the same chromosome as the "
                        "target gene/protein - used to model family-relationships and population structure")

    parser.add_argument("--grm-id", dest="grm_id_file", type=str, default=None,
                        help="File of smaple IDs accompanying the GRM file")

    parser.add_argument("--start", dest="start_pos", type=int,
                        help="Starting position along the chromsome of the cis-window")

    parser.add_argument("--end", dest="end_pos", type=int,
                        help="Ending position along the chromosome of the cis-window")

    parser.add_argument("--maf", dest="maf_thresh", type=float,
                        help="Minimal MAF at which to test SNPs")

    parser.add_argument("--exclude-snps", dest="exclude_file", type=str,
                        help="File of SNPs to exclude, one SNP per line")

    parser.add_argument("--log", dest="logfile", type=str, default=sys.stdout,
                        help="PATH to logging file")

    parser.add_argument("--output", dest="outfile", type=str, default=sys.stdout,
                        help="File PATH to write results to")

    parser.add_argument("--permutation", dest="perms", action='store_true',
                        help="Perform permutation-based per-gene multiple testing correction")

    parser.add_argument("--perm-method", dest="perm_method", type=str,
                        choices=["beta", "fixed"], default="beta",
                        help="Permutation method to use, either fixed iterations, or beta-approximation")

    parser.add_argument("--iters", dest="perm_n", type=int,
                        default=100, help="Number of permutations to use for gene-wise multiple testing procedure")

    args = parser.parse_args()
    # setup the logger
    if type(args.logfile) == str:
        logging.basicConfig(level=logging.INFO,
                            filename=args.logfile)
    else:
        # logging defaults to stdout
        logging.basicConfig(level=logging.INFO)

    # read in the covariates file
    logging.info("Reading in covariates file {}".format(args.covar_file))
    covar_df = pd.read_table(args.covar_file, sep="\t", header=None)
    cov_names = ["COV{}".format(i) for i in range(covar_df.shape[1]-2)]
    covar_df.columns = [x for sub in [["FID", "IID"], cov_names] for x in sub]
    logging.info("Found {} covariates".format(covar_df.shape[1]-2))

    # read in the plink genotype files - assumed to be a single chromosome
    logging.info("Reading plink files: {}".format(args.geno_path))
    (bim, fam, bed) = limix.io.plink.read(args.geno_path, verbose=False)

    # read in PHENOTYPE and merge with FAM file - important for subsetting the 
    # correct individual-level genotypes - merge with covariates
    logging.info("Reading phenotype file: {}".format(args.pheno_file))
    pheno_df = pd.read_table(args.pheno_file, sep="\t", header=None)
    pheno_df.columns = ["FID", "IID", "PHENO"]
    pheno_df = pd.merge(pheno_df, covar_df, on=["FID", "IID"])

    # remove infinite and missing values individuals
    missing_cov = pheno_df.loc[:, cov_names].apply(lambda X: not X.isnull().values.any(), axis=1)
    pheno_df = pheno_df.loc[missing_cov, :]

    # I need to cast the IID to the same type
    fam_iid_type = fam.iid.dtype
    pheno_iid_type = pheno_df.IID.dtype

    if fam_iid_type != pheno_iid_type:
        logging.warning("Casting FAM IID {} to type: {}".format(fam_iid_type, pheno_iid_type))
        fam.loc[:, "iid"] = fam.iid.astype(pheno_iid_type, copy=True)
    else:
        pass

    logging.info("Matching individuals across files")
    fam_merge = pd.merge(fam, pheno_df, left_on=["iid"], right_on=["IID"])
    fam_merge.trait = fam_merge.PHENO

    logging.info("{} individuals matching between genotypes and phenotype".format(fam_merge.shape[0]))
    covar_merge = fam_merge.loc[:, cov_names]

    # have to set sample IDs as the indices
    covar_merge.index = fam_merge.iid
    fam_merge = fam_merge.loc[:, ["fid", "iid", "father", "mother", "gender", "trait", "i"]]
    fam_merge.index = fam_merge.iid

    logging.info("{} individuals matching between phenotypes and covariates".format(covar_merge.shape[0]))
    # need to subset the genotypes to just the cis-window
    # the bed and bim file are in the same order - therefore indexing the bim will 
    # also index the appropriate positions in the bed file

    test_chr = args.geno_path.split("/")[-1].split("-")[0]
    logging.info("Extracting genotypes from {}:{}-{}".format(test_chr, args.start_pos, args.end_pos))
    bim_cis = bim.query("(pos >= {}) & (pos <= {})".format(args.start_pos, args.end_pos))
    bed_cis = bed[bim_cis.i.values, ].compute()

    # calculate MAFs
    logging.info("Calculating MAFs")
    cis_maf = limix.qc.compute_maf(bed_cis.T)
    keep_geno = cis_maf >= args.maf_thresh
    logging.info("Retaining {} SNPs with MAF >= {}".format(sum(keep_geno), args.maf_thresh))

    logging.info("Subsetting genotypes to {} indviduals".format(fam_merge.shape[0]))
    # subset to individuals
    bed_cis = bed_cis[keep_geno, :]
    bed_cis = bed_cis[:, fam_merge.i.values]
    bim_cis = bim_cis.loc[keep_geno, :]
    bim_cis.reset_index(inplace=True, drop=False)
    logging.info("Extracted genotypes for {} variants".format(bed_cis.shape[0]))

    # remove genotypes with non-finite values
    # set missing values to -9?
    nan_genos = sum(np.apply_along_axis(func1d=lambda X: np.any(np.isnan(X)), axis=1, arr=bed_cis))
    if nan_genos > 0:
        logging.warning("{} non-finite genotype values found, setting to -9".format(nan_genos))
        bed_cis[np.isnan(bed_cis)] = -9
    else: pass

    logging.info("Reading GRM: {}".format(args.grm_file))
    (grm, nsnps) = pandas_plink.read_grm(filepath=args.grm_file, id_filepath=args.grm_id_file)

    logging.info("Subsetting individauls from GRM")
    grm_cis = grm[fam_merge.i.values, fam_merge.i.values]

    # check arrays are concordant
    logging.info("Check matrices are concordant")
    nonzero = bed_cis.shape[1]
    concord_mats = (bed_cis.shape[1] == fam_merge.shape[0]) & (bim_cis.shape[0] == bed_cis.shape[0]) & (bed_cis.shape[1] == grm_cis.shape[0])

    if concord_mats:
        pass
    else:
        raise ValueError("Matrices are not concordant - check overlapping individuals")

    if nonzero > 0:
        pass
    else:
        raise Exception("No individuals are overlapping between files - exiting")
        exit

    logging.info("Matrices are concordant: {}".format(concord_mats))

    if args.perms:
        logging.info("Permforming null distribution estimation by phenotype permutation")
        # this function returns a CDF for all SNPs
        correct_func = perm_qtl(genos=bed_cis, phenos=fam_merge.trait.values,
                                covars=covar_merge, grm=grm_cis, bim=bim_cis,
                                method="beta", n_iters=args.perm_n)
    else: pass

    ## Setup the linear mixed model
    # the G matrix needs to be transposed, i.e. nxp
    logging.info("Setting up LMM")
    lmm = limix.qtl.scan(G=bed_cis.T, Y=fam_merge.trait.values, 
                         lik="normal", K=grm_cis, M=covar_merge,
                         verbose=False)

    # add the SNP names from the bim file
    res_df = lmm.stats
    res_df = pd.merge(res_df, bim_cis, left_index=True, right_index=True)
    # extract the SNP effect sizes
    coef_df = lmm.effsizes['h2']
    #print(coef_df.loc[coef_df.effect_type == 'candidate'].shape)
    cov_eff_df = coef_df.loc[coef_df.effect_name.isin(cov_names), :]
    coef_df = coef_df.loc[~coef_df.effect_name.isin(cov_names), :]
    coef_df.set_index('test', drop=False, inplace=True)

    res_df = pd.merge(res_df, coef_df, left_index=True, right_index=True)
    res_df = res_df.loc[:, ["chrom", "pos", "a0", "a1", "snp", "effsize", "effsize_se", "pv20"]]
    res_df.columns = ["CHR", "BP", "A1", "A2", "SNP", "BETA", "SE", "P"]

    if args.perms:
        logging.info("Getting multiple-testing corrected p-values")
        corrected_pvals = np.array([correct_func(x) for x in res_df.P.values])
        res_df["PermP"] = corrected_pvals
    else: pass

    # output the relevant information: SNP, A1, A2, DF, PVAL, CHR, SNP, BP, BETA, SE, PCorrect
    res_df.to_csv(args.outfile, sep="\t", index=False, index_label=False)

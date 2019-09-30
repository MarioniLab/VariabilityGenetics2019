#! /usr/bin/python3

import limix
import pandas_plink
import sys
import numpy as np
import pandas as pd
import argparse
import logging
import re
import itertools


def file_len(fname):
    '''
    Check the number of lines in file quickly
    '''

    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return(i + 1)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="A script for decomposing the variance of a trait into multiple genetic components")

    parser.add_argument("--phenotype", dest="pheno_file", type=str,
                        help="The file containing a vector of phenotypes and sample IDs. Missing values should be NA")

    parser.add_argument("--covariates", dest="covar_file", type=str, default=None,
                        help="A table of covariate fixed effects to include in the LMM, with FID and IID columns to identify individuals")

    parser.add_argument("--grm", dest="grm_file", type=str,
                        help="The file path to the GRM file for estimating the genetic variance component")

    parser.add_argument("--grm-id", dest="grm_id_file", type=str, default=None,
                        help="File of sample IDs accompanying the GRM file. Leave blank to derive ID file name from parent GRM file")

    parser.add_argument("--multi-grm", dest="grm_list", type=str,
                        default=None,
                        help="Path to a text file that contains two columns and one GRM per row. The first is a label for the GRM, the second is the path of the GRM file. This is used for decomposing the trait variance into multiple components, e.g. cis & trans. If there is only 1 column then the label is taken from the filename.")

    parser.add_argument("--log", dest="logfile", type=str, default=sys.stdout,
                        help="Path to logging file")

    parser.add_argument("--output", dest="outfile", type=str, default=sys.stdout,
                        help="File PATH to write results to")

    parser.add_argument("--model-common-environment", dest="model_env", action="store_true",
                        help="Include a variance component in the model that takes into account the shared familial environment between individuals")

    args = parser.parse_args()
    # setup the logger
    if type(args.logfile) == str:
        logging.basicConfig(level=logging.INFO,
                            filename=args.logfile)
    else:
        logging.basicConfig(level=logging.INFO)

    if args.covar_file:
        logging.info("Reading covariates file {}".format(args.covar_file))
        covar_df = pd.read_table(args.covar_file, sep="\t", header=None)
        cov_names = ["COV{}".format(i) for i in range(covar_df.shape[1]-2)]
        covar_df.columns = [x for sub in [["FID", "IID"], cov_names] for x in sub]
        logging.info("Found {} covariates".format(covar_df.shape[1]-2))
    else: pass

    # read in PHENOTYPE and merge with FAM file - important for subsetting the 
    # correct individual-level genotypes - merge with covariates
    logging.info("Reading phenotype file: {}".format(args.pheno_file))
    pheno_df = pd.read_table(args.pheno_file, sep="\t", header=None)
    pheno_df.columns = ["FID", "IID", "PHENO"]

    if args.covar_file:
        pheno_df = pd.merge(pheno_df, covar_df, on=["FID", "IID"])

        # remove infinite and missing values individuals
        missing_cov = pheno_df.loc[:, cov_names].apply(lambda X: not X.isnull().values.any(), axis=1)
        pheno_df = pheno_df.loc[missing_cov, :]

    if args.grm_file and not args.grm_list:
        logging.info("Reading GRM: {}".format(args.grm_file))
        (grm, nsnps) = pandas_plink.read_grm(filepath=args.grm_file, id_filepath=args.grm_id_file)

        logging.info("Subsetting individauls from GRM")
        common_iid = [x for x in set(pheno_df.IID.values).intersection(grm['sample_0'].values)]
        grm_cis = grm.loc[common_iid, common_iid]
        #print(grm_cis[:5, :5])
        pheno_df.set_index("IID", drop=False, inplace=True)
        final_pheno = pheno_df.loc[common_iid, :]
    elif not args.grm_file and not args.grm_list:
        raise IOError("Missing GRM files. Please provide at least one GRM for genetic variance decomposition")
    elif not args.grm_file and args.grm_list:
        # check a label first, then a valid path to the GRMs
        # if there is no label then use the file name
        # may need to append '.grm.bin' as well if not present
        grm_list = {}
        with open(args.grm_list, "r") as gfile:
            file_size = file_len(args.grm_list)
            logging.info("Found {} GRM files in {}".format(file_size, args.grm_list))
            for line in gfile.readlines():
                if len(line.rstrip("\n").split("\t")) > 1:
                    label = line.split("\t")[0]
                    grm_file = line.split("\t")[-1].rstip("\n") + ".grm.bin"
                    logging.info("Reading GRM: {}".format(grm_file))
                    (grm, nsnps) = pandas_plink.read_grm(filepath=grm_file, id_filepath=None)

                    logging.info("Subsetting individauls from GRM")
                    common_iid = [x for x in set(pheno_df.IID.values).intersection(grm['sample_0'].values)]
                    grm_cis = grm.loc[common_iid, common_iid]

                elif len(line.rstrip("\n").split("\t")) == 1:
                    grm_file = line.rstrip("\n") + ".grm.bin"
                    label = grm_file.split("/")[-1]
                    (grm, nsnps) = pandas_plink.read_grm(filepath=grm_file, id_filepath=None)
                    logging.info("Subsetting individauls from GRM")
                    common_iid = [x for x in set(pheno_df.IID.values).intersection(grm['sample_0'].values)]
                    grm_cis = grm.loc[common_iid, common_iid]

                grm_list[label] = grm_cis
        pheno_df.set_index("IID", drop=False, inplace=True)
        final_pheno = pheno_df.loc[common_iid, :]       
    
    # check arrays are concordant
    logging.info("Check matrices are concordant")
    concord_mats = (final_pheno.shape[0] == grm_cis.shape[0])
    nonzero = final_pheno.shape[0]

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
    logging.info("Setting up LMM for variance decomposition")
    #print(final_pheno.head())
    if args.covar_file:
        logging.info("Adjusting for {} covariates".format(len(cov_names)))
        vardec = limix.vardec.VarDec(final_pheno.PHENO, "normal", final_pheno.loc[:, cov_names])
    else:
        vardec = limix.vardec.VarDec(final_pheno.PHENO, "normal")

    # append each grm
    trans_grm = 0
    for k in grm_list.keys():
        logging.info("Adding {} to LMM".format(k))
        # vardec doesn't like the xarray.DataArray format
        # cast to numpy array with .values
        vardec.append(grm_list[k].values, k)
        if re.search("Not", k):
            trans_grm = grm_list[k].values

    if args.model_env:
        logging.info("Generating a common environmental matrix - adding to model")
        # this should be obvious from the genotype data
        # calculate from the trans GRM, the cis is too variable
        env_matrix = np.absolute(np.around(trans_grm, 1))
        #print(env_matrix[:25, :25])
        # convert all 1's and 0.5's to 1's
        env_matrix[env_matrix > 0] = 1
        #print(env_matrix[:25, :25])

        logging.info("Modelling environment for {} families".format(int(np.sum(env_matrix == 1)/2.0)))
        vardec.append(env_matrix, "ENV")

    logging.info("Adding residual covariance as IID noise")
    vardec.append_iid("noise")
    
    logging.info("Decomposing variance")
    vardec.fit(verbose=False)

    #print(vardec)
    # extract variance estimates from .GivenCov objects
    # how do I get the standard errors from the variance estimates?
    # I assume this an MLE - so it should be possible?
    cov_dict = {}
    # check the variances sum to 1
    trait_var = np.var(final_pheno.PHENO)
    variances = [v.scale for v in vardec._covariance]
    sumvar = sum(variances)
    # the sum of variances should be the same as the total trait variance
    #print(sumvar)
    #print(trait_var)
    for x in vardec._covariance:
        if re.search("Not", x._name):
            cov_dict["trans"] = x.scale/sumvar
        elif re.search("noise", x._name):
            cov_dict["noise"] = x.scale/sumvar
        elif re.search("ENV", x._name):
            cov_dict["env"] = x.scale/sumvar
        else:
            cov_dict["cis"] = x.scale/sumvar
        #print(x.__dict__)
    #print(sum(cov_dict.values()))
    try:
        np.testing.assert_approx_equal(sumvar, trait_var, significant=2)
        logging.info("Total variance sum {} is approximately equal to trait variance {}".format(sumvar, trait_var))
    except AssertionError:
        logging.error("Total variance not approximately equal to {}: {}".format(trait_var, sumvar))

    # information to write out:
    # phenotype name
    # cis variance estimate
    # trans variance estimate
    # residual variance estimate
    # sample size

    trait = args.pheno_file.split("/")[-1].rstrip(".pheno")
    trait_var = np.var(final_pheno.PHENO)
    # variance compnents are in absolute terms
    #print(trait)
    if type(args.outfile) == str:
        with open(args.outfile, "w") as ofile:
            if args.model_env:
                ofile.write("Trait\tcis.H2\ttrans.H2\tCommonEnvironment\tResidual\tNSize\tTraitVar\n")
                ofile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t\n".format(trait, cov_dict["cis"],
                                                                    cov_dict["trans"], cov_dict["env"], 
                                                                    cov_dict["noise"], final_pheno.shape[0],
                                                                    trait_var))
            else:
                ofile.write("Trait\tcis.H2\ttrans.H2\tResidual\tNSize\tTraitVar\n")
                ofile.write("{}\t{}\t{}\t{}\t{}\t{}\t\n".format(trait, cov_dict["cis"],
                                                                cov_dict["trans"],
                                                                cov_dict["noise"],
                                                                final_pheno.shape[0], trait_var))

    else:
        ofile = args.outfile
        if args.model_env:
            ofile.write("Trait\tcis.H2\ttrans.H2\tCommonEnvironment\tResidual\tNSize\tTraitVar\n")
            ofile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(trait, cov_dict["cis"],
                                                              cov_dict["trans"], cov_dict["env"],
                                                              cov_dict["noise"],
                                                              final_pheno.shape[0], trait_var))
        else:
            ofile.write("Trait\tcis.H2\ttrans.H2\tResidual\tNSize\tTraitVar\n")
            ofile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(trait, cov_dict["cis"],
                                                          cov_dict["trans"],
                                                          cov_dict["noise"],
                                                          final_pheno.shape[0], trait_var))
            

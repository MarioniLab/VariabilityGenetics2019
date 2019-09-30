import os
import re
import sys
import subprocess
import argparse
import logging
import pandas as pd

parser = argparse.ArgumentParser(description="Extract lead SNP results from .mlma files. Additionally add rsIDs.")
parser.add_argument("--output", dest="output", type=str,
                    default=sys.stdout, help="Output file path")

parser.add_argument("--results-directory", dest="input_dir", type=str,
                    help="The directory containing the .mlma results files")

parser.add_argument("--clump-file", dest="clump_file", type=str,
                    help="The clumping results for a trait/phenotype")

parser.add_argument("--trait", dest="trait", type=str,
                    choices=["Mean", "Noise", "Joint"],
                    help="The specific trait phenotype")

parser.add_argument("--celltype", dest="celltype", type=str,
                    help="The cell type for the expression trait")

parser.add_argument("--protein", dest="protein", type=str,
                    help="The protein for the expression trait")

parser.add_argument("--panel", dest="panel", type=str,
                    help="The relevant antibody panel for the trait")

parser.add_argument("--log", dest="logfile", type=str,
                    help="Logging file destination",
                    default=sys.stdout)

parser.add_argument("--SNP-files", dest="id_file", type=str,
                    default=False,
                    help="PATH to files containing chr:position information and SNP "
                    "IDs - split by chromosome")

args = parser.parse_args()

# setup the logger

if type(args.logfile) == str:
    logging.basicConfig(level=logging.INFO,
                    filename=args.logfile)
else:
    logging.basicConfig(level=logging.INFO)


logging.info("Reading clumped results file")
try:
    # the input clumped file should be able to fit in memory - need a variable length whitespace delimiter
    clump_df = pd.read_table(args.clump_file, header=None, index_col=None, sep=r"\s*")
except pd.errors.EmptyDataError:
    logging.error("Clump file is empty - exiting")
    sys.exit()

if clump_df.shape[0] == 0:
    raise IOError("Clump file is empty")
else:
    pass

clump_df.columns = ["CHR", "F", "SNP", "BP", "P", "TOTAL", "NSIG", "S05", "S01", "S001", "S0001", "SP2"]
clump_df.set_index(keys="SNP", inplace=True, drop=False)
# remove duplicate entries based on the SNP Id
clump_df.drop_duplicates(keep='first', subset=["CHR", "F", "SNP", "BP"], inplace=True)
#print(clump_df.head())

# iterate over the lead SNPs and extract the relevant results from the .mlma file
search_trait = "-".join([args.trait, "_".join([args.protein, args.celltype, args.panel])])

logging.info("Searching for matching mlma results files")
match_files = [px for px in os.listdir(args.input_dir) if re.search(search_trait, px)]
mlma_files = [bz for bz in match_files if re.search("bgz$", bz)]
#print(mlma_files[:5])

if len(mlma_files) == 0:
    raise ValueError("No matching mlma results files were found - check file path and contents")
else:
    pass

logging.info("Found {} matching results files".format(len(mlma_files)))

clump_list = []
logging.info("Extracting results for {} lead variants".format(clump_df.shape[0]))
#print(clump_df.head())
for x in clump_df["SNP"].values:
    bp = clump_df.loc[x, "BP"].astype(int)
    chr_int = clump_df.loc[x, "CHR"]
    chrome = "Chr" + str(clump_df.loc[x, "CHR"])
    #print(clump_df.loc[x, ])

    # search the results directory for the appropriate chromosome
    if chrome == "chr1":
        chrome = "chr01"
    chrome_file = [cx for cx in mlma_files if re.search(chrome + ".mlma", cx)]
    #print(chrome)
    #print(chrome_file)
    # tabix query based on the bp position
    tab_query = "tabix {}{} {}:{}-{}".format(args.input_dir, chrome_file[0], chr_int,
                                               bp, bp)
    logging.info("Running tabix query {}".format(tab_query))
    proc = subprocess.Popen(tab_query, shell=True, stdout=subprocess.PIPE)
    proc_channel = proc.communicate()
    proc_rc = proc.returncode

    if proc_rc:
        logging.info("Tabix query exit code is {}".format(proc_rc))
        raise SystemError("tabix query executed with non-zero status - check file format")
    else:
        # I need to check that the exit status is 0
        logging.info("Tabix query exit code is {}".format(proc_rc))
        
        result_dict = {}
        count = 1
        # executing proc.communicate() returns all of the output in a tuple
        # the tuple is (stdout, stderr)
        for line in proc_channel[0].decode("utf-8").split("\n"):
            if line == '':
                pass
            else:
                snp_dict = {}
                parse = line.split("\t")
                #print(parse)
                try:
                    if int(parse[2]) > bp + 1:
                        break
                    else:
                        snp_dict["CHR"] = int(parse[0])
                        snp_dict["SNP"] = parse[1]
                        snp_dict["BP"] = int(parse[2])
                        snp_dict["A1"] = parse[3]
                        snp_dict["A2"] = parse[4]
                        snp_dict["MAF"] = float(parse[5])
                        snp_dict["BETA"] = float(parse[6])
                        snp_dict["SE"] = float(parse[7])
                        snp_dict["P"] = float(parse[8])
                        count += 1
    
                        result_dict[count] = snp_dict
                    #print(result_dict)
                    #print(args.id_file)
                    if args.id_file:
                        # need to run tabix query to extract the appropriate rsIDs too
                        # find the correct chromosome BED file file first
                        try:
                            chr_bed = [bx for bx in os.listdir(args.id_file) if re.search("chr{}.bed.gz$".format(int(parse[0])), bx)][0]
                        except IndexError:
                            logging.warn("No matching BED file detected for chromosome {}".format(parse[0]))

                        bed_tabix_query = "tabix {}/{} chr{}:{}-{}".format(args.id_file, chr_bed, int(parse[0]),
                                                                           int(parse[2])-1, int(parse[2]))
                
                        logging.info("Running BED tabix query {}".format(bed_tabix_query))

                        bed_proc = subprocess.Popen(bed_tabix_query, shell=True, stdout=subprocess.PIPE)
                        bed_proc_channel = bed_proc.communicate() 
                        bed_proc_rc = bed_proc.returncode                                                                                                               

                        if bed_proc_rc:
                            logging.info("Tabix BED query exit code is {}".format(bed_proc_rc))
                            raise SystemError("tabix query executed with non-zero status - check file format")                                                          
                        else:                                                                                                                                                        
                            # I need to check that the exit status is 0
                            logging.info("Tabix query exit code is {}".format(bed_proc_rc))

                            # parse the BED tabix query output
                            bed_result_dict = {}
                            count = 1                                                                                                                                              
                            # executing proc.communicate() returns all of the output in a tuple
                            # the tuple is (stdout, stderr)
                            for line in bed_proc_channel[0].decode("utf-8").split("\n"):
                                bed_snp_dict = {}
                                bed_parse = line.split("\t")

                                try:
                                    if int(bed_parse[2]) > bp + 1:
                                        break
                                    else:
                                        bed_snp_dict["CHR"] = int(bed_parse[0].strip("chr"))
                                        bed_snp_dict["SNP"] = bed_parse[3]
                                        bed_snp_dict["BP"] = int(bed_parse[1]) + 1
                                        bed_snp_dict["BPEND"] = int(bed_parse[2])
                                        bed_snp_dict["STRAND"] = bed_parse[5]
                                        count += 1
                                        bed_result_dict[count] = bed_snp_dict
                            
                                        if count > 1:
                                            break
                                except IndexError:
                                    logging.warn("No matching SNP id found for {}:{}".format(parse[0], int(parse[2])-1))
                except IndexError:
                    logging.warn("No SNP result found")
            
            lead_df = pd.DataFrame(result_dict).T

            if args.id_file:
                rsid_df = pd.DataFrame(bed_result_dict).T
        
                # merge rsID information
            if args.id_file:
                try:
                    lead_df = pd.merge(lead_df, rsid_df, left_on=['CHR', 'BP'], right_on=['CHR', 'BP'], how='outer')
                    clump_list.append(lead_df)
                except:
                    logging.warn("Matching rsID not found at {}".format(lead_df["CHR:POS"].values))
                    # need to add the empty columns anyway
                    lead_df["SNP"] = "NA"
                    lead_df["BPEND"] = lead_df["BP"]
                    lead_df["STRAND"] = "NA"
            else:
                clump_list.append(lead_df)

if len(clump_list) > 0:
    clumped_lead_df = pd.concat(clump_list)
    clumped_lead_df.drop_duplicates(inplace=True, subset=["SNP", "BETA"], keep='first')
    # remove all SNPs with MAF <= 4%
    clumped_lead_df = clumped_lead_df.loc[clumped_lead_df["MAF"] >= 0.04, :]
    #print(clumped_lead_df)
    # merge this with the input clump file
    # the p-value column from the results file is more accruate than from the clumped
    logging.info("Merging results and clumps")
    select_cols = ["SNP", "CHR", "BP", "F", "NSIG", "TOTAL", "S05", "S01",
                   "S001", "S0001", "SP2"]
    clump_df = clump_df.loc[:, select_cols]

    # later versions of pandas struggle if a column and index are the same
    clump_df.reset_index(inplace=True, drop=True)
    #print(clump_df.head())
    #print(clumped_lead_df.head())
    if args.id_file:
        clump_df.columns = ["CHR:POS", "CHR", "BP", "F", "NSIG", "TOTAL", "S05", "S01",
                            "S001", "S0001", "SP2"]
        merged_clump = pd.merge(clumped_lead_df, 
                                clump_df,
                                how='left',
                                left_on=["BP", "CHR", "CHR:POS"],
                                right_on=["BP", "CHR", "CHR:POS"])
    else:
        merged_clump = pd.merge(clumped_lead_df,
                                clump_df,
                                how='left',
                                left_on=["BP", "CHR", "SNP"],
                                right_on=["BP", "CHR", "SNP"])
        
    merged_clump.to_csv(args.output, sep="\t", index=None)

else:
    logging.warn("No lead SNP results found - exiting")

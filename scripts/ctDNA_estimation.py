# Estimate cfDNA using mutations
import pandas as pd
import numpy as np 
import os
from scipy.stats import binom

DIR_working = "/home/vpc/datasets/plasma_chipseq/pipeline"

# Inputs
PATH_mutations_curated = os.path.join(DIR_working, "analysis/all_mutations_curated.tsv")
PATH_gene_CN = os.path.join(DIR_working, "analysis/CN.csv")
PATH_sample_dict = os.path.join(DIR_working, "sample_lists/naming_dictionary.tsv") # to help translate takeda ids into GUBB IDs

sheet_id = "1W3Xw_t56dxOY3iHJ3FXW_L2wT8DVQ8GoLbapwhzRjFo"
sheet_name = "WGD_[ARCHIVED]"
URL_CN_ctDNA = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"

# outputs 
PATH_mutation_ctfractions = os.path.join(DIR_working, "analysis/ct_fractions.tsv")

muts = pd.read_csv(PATH_mutations_curated, "\t")
cna = pd.read_csv(PATH_gene_CN)
sample_dict = pd.read_csv(PATH_sample_dict, "\t")
cn_ct = pd.read_csv(URL_CN_ctDNA)[["Takeda_Patient", "CN_ct", "Diploid level", "WGD"]]

def get_adj_vaf(vaf, depth):
    """
    Adjust the mutation VAF using binomial distribution.
    """
    alt = vaf*depth
    for p in np.arange(0,1,0.005):
        dist = binom(depth, p)
        if dist.cdf(alt) < 0.95:
            vaf = p; break;
    return(vaf)

muts = muts.merge(cna[["Patient", "Sample", "Gene", "Coverage logratio"]], "left")
muts = (muts[(muts['Coverage logratio'] < 0.3) &
             (muts["Characteristic"] != "Germline") &
             (muts['Total_reads'] > 40) &
             (~muts['Chrom'].str.contains("X")) &
             (~muts['Chrom'].str.contains("Y"))])

muts['Adj_allele_frequency'] = muts.apply(lambda row: get_adj_vaf(row['VAF'], row['Total_reads']), axis = 1)
muts["Max_Sample_adj_VAF"] = muts.groupby('Sample')["Adj_allele_frequency"].transform('max')
muts = muts[muts["Adj_allele_frequency"] == muts["Max_Sample_adj_VAF"]]
muts['row_max'] = muts[['Total_reads']].max(axis=1) #To break any ties for max VAF!
muts = muts.loc[muts.groupby('Sample')['row_max'].idxmax()]

muts['Mutation_ctDNA_fraction'] = 2/(1+1/(muts['Max_Sample_adj_VAF']))

muts = muts.merge(sample_dict, how = "outer") # some pts have no mutations available to call CT fractions, this retains them and allows a helpful comparison between CN based and mut based estimates across the whole cohort.
# reorder columns after the merge
columns = muts.columns.tolist()
new_columns = [columns[0]] + columns[-3:] + columns[1:-3] # brings the last two columns into second and third positions
muts = muts[new_columns]
muts = muts.merge(cn_ct, "left") # add the CN based ct estimates too
muts["CN_ct"] = muts["CN_ct"]/100 # mut based ct is fraction, so this one should be as well
muts["Mut ct-CN ct"] = muts["Mutation_ctDNA_fraction"] - muts["CN_ct"]
muts["Mut ct-CN ct (PERCENT)"] = muts["Mut ct-CN ct"]*100
muts['Final_ct_fraction'] = muts.apply(lambda row: row['CN_ct'] if row['CN_ct'] < 0.2 and pd.isnull(row['Mutation_ctDNA_fraction']) else row['Mutation_ctDNA_fraction'], axis=1) # if CN based estimate is less than 20, use the mutation estimate
muts['Final_ct_fraction (PERCENT)'] = muts['Final_ct_fraction']*100
muts.to_csv(PATH_mutation_ctfractions, "\t", index = False) # this is the file that shows the mutation ctDNA fraction estimates, along with the mutation used to estimate it
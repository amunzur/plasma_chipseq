# This scr

import os 
import pandas as pd
import numpy as np
import re

cn_dict_category = {
    "Deep_deletion": (-np.inf, -1),
    "Deletion": (-1, -0.3),
    "Neutral": (-0.3, 0.3),
    "Gain": (0.3, 0.7),
    "Amplification": (0.7, np.inf)
}

def subset_mutations_to_panel(PATH_mutation_vcf, PATH_target_vcf):
    """
    Given a mutations file and a bed file with the target regions, only keep mutations falling within the bed file.
    """
    muts = pd.read_csv(PATH_mutation_vcf, "\t")
    target = pd.read_csv(PATH_target_vcf, "\t", header = None)
    target.columns = ["Chrom", "Start", "Stop", "Gene"]
    mutations_within_target = []
    for _, mutation in muts.iterrows():
        chrom = mutation['CHROM']
        position = mutation['POSITION']
        target_matches = target[(target['Chrom'] == chrom) & (target['Start'] <= position) & (target['Stop'] >= position)]
        if not target_matches.empty:
            mutations_within_target.append(mutation)
    mutations_within_target_df = pd.DataFrame(mutations_within_target)
    return(mutations_within_target_df)    

DIR_project="/home/vpc/datasets/plasma_chipseq/pipeline"
PATH_altering=os.path.join(DIR_project, "mutations", "somatic_protein_altering.vcf")
PATH_silent=os.path.join(DIR_project, "mutations", "somatic_silent.vcf")
PATH_germline=os.path.join(DIR_project, "mutations", "rare_germline.vcf")
PATH_ct_fraction=os.path.join(DIR_project, "ctdna_fractions.tsv")
PATH_gene_CN=os.path.join(DIR_project, "gene_cna_prostate.tsv")
PATH_target=os.path.join(DIR_project, "target_probes", "PROSTATE_target.bed")

PATH_mutations=os.path.join(DIR_project, "analysis", "mutations.csv")
PATH_save_CN=os.path.join(DIR_project, "analysis", "CN.csv")
PATH_save_ct=os.path.join(DIR_project, "analysis", "ct.csv")

PATH_mutations_silent_altering_PROCESSED=os.path.join(DIR_project, "analysis", "silent_and_altering_mutations.csv")
PATH_mutations_germline_PROCESSED=os.path.join(DIR_project, "analysis", "germline_mutations.csv")

# silent and altering somatic only to begin with
mutations_list = []
for path in [PATH_altering, PATH_silent, PATH_germline]: 
    muts = subset_mutations_to_panel(path, PATH_target)
    muts = pd.melt(muts, id_vars=['CHROM', 'POSITION', 'REF', 'ALT', 'GENE', 'EFFECT', 'NOTES'], var_name='SAMPLES', value_name='MUTATION')
    muts = muts.loc[muts['MUTATION'].str.contains('\*', regex=True)].reset_index(drop=True)
    muts[["EFFECT", "PROTEIN_ALTERATION"]] = muts["EFFECT"].str.split(" ", 1, expand = True)
    muts["PATIENT"] = muts["SAMPLES"].apply(lambda x: re.split("-cfDNA|_cfDNA|-WBC|_WBC", x)[0]) # get patient ID from sample name
    if path != PATH_germline:
        muts[["Alt_reads", "Total_reads", "MQ", "Strand_bias", "Sidedness"]] = muts["MUTATION"].str.split(":", 4, expand = True)
        muts["Sidedness"] = muts["Sidedness"].str.replace(":\*", "") # remove extra characters
        muts["VAF"] = muts["Alt_reads"].astype(int)/muts["Total_reads"].astype(int)
        muts = muts[["PATIENT", "SAMPLES", "CHROM", "POSITION", "REF", "ALT", "GENE", "EFFECT", "PROTEIN_ALTERATION", "NOTES", "Alt_reads", "Total_reads", "VAF", "MQ", "Strand_bias", "Sidedness"]]
        muts.columns = ["Patient", "Sample", "Chrom", "Position", "Ref", "Alt", "Gene", "Effect", "Protein_alteration", "Notes", "Alt_reads", "Total_reads", "VAF", "MQ", "Strand_bias", "Sidedness"]
        if path == PATH_silent: 
            muts = muts[muts["MQ"].astype(int) >= 30].reset_index(drop = True)
        elif path == PATH_altering: 
            muts = muts[muts["MQ"].astype(int) >= 10].reset_index(drop = True) 
    else:
        muts = muts[["PATIENT", "SAMPLES", "CHROM", "POSITION", "REF", "ALT", "GENE", "EFFECT", "PROTEIN_ALTERATION", "MUTATION", "NOTES"]]
        muts.columns = ["Patient", "Sample", "Chrom", "Position", "Ref", "Alt", "Gene", "Effect", "Protein_alteration", "Mutation", "Notes"]
    if path == PATH_altering:
        muts["Characteristic"] = "Altering"
    elif path == PATH_germline: 
        muts["Characteristic"] = "Germline"
        muts = muts[muts["Notes"].str.contains("pathogen", na=False)].reset_index(drop = True) #only keep pathogenic germline mutations
    elif path == PATH_silent: 
        muts["Characteristic"] = "Silent"
    mutations_list.append(muts)

# These are for google sheets
mutations_silent_altering = pd.concat(mutations_list[0:2]).reset_index(drop = True)
mutations_silent_altering.to_csv(PATH_mutations_silent_altering_PROCESSED, index = False)
mutations_list[2].to_csv(PATH_mutations_germline_PROCESSED, index = False)

# COPY NUMBER
cn = pd.read_csv(PATH_gene_CN, "\t")[['Sample', 'Gene', 'Coverage logratio']]
cn['CN event'] = cn['Coverage logratio'].apply(lambda x: next((category for category, (start, end) in cn_dict_category.items() if start <= x <= end), "Unknown"))
cn["Patient"] = cn["Sample"].apply(lambda x: re.split("-cfDNA|_cfDNA|-WBC|_WBC", x)[0]) # get patient ID from sample name
cn = cn[["Patient", "Sample", "Gene", "Coverage logratio", "CN event"]]

# mutations["Gene alteration percentage"] = (mutations.groupby(["Gene", "Effect"])["Sample"].transform("count")/cn["Sample"].unique().size)*100

# TUMOR FRACTION
ct = pd.read_csv(PATH_ct_fraction, sep="\t", header=None)
ct.columns = ["Sample", "ct percentage"]

# Save to path
ct.to_csv(PATH_save_ct, index = False)
cn.to_csv(PATH_save_CN, index = False)
# mutations.to_csv(PATH_mutations, index = False)


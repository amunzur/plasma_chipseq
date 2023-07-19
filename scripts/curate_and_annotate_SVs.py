# This script 
import os
import pandas as pd
import numpy as np
import re

DIR_SVs = "/home/vpc/datasets/plasma_chipseq/pipeline/rearrangements/annotated"
PATH_bed_prostate = "/home/vpc/datasets/plasma_chipseq/pipeline/target_probes/PROSTATE_target.bed" # need to have gene names
PATH_bed_bladder = "/home/vpc/datasets/bladder_uc/bcap21/probes_hg38.bed" # need to have gene names
PATH_prostate_patients = "/home/vpc/datasets/plasma_chipseq/pipeline/sample_lists/prostate_patients.txt"
PATH_bladder_patients = "/home/vpc/datasets/plasma_chipseq/pipeline/sample_lists/bladder_patients.txt"

def annotate_sv_type(sv):
    """
    Given an sv file, annotate the type of SV. 
    """
    sv['CATEGORY'] = "SV" #Dummy entry to manually check that all SVs have received a proper annotation
    for index, row in sv.iterrows():
        if row['CHROM'] == row['CHROM.1']:
            if row['STRAND'] == row['STRAND.1']:
                if (row['STRAND'] == '+') & (row['POSITION'] < row['POSITION.1']):
                    sv.at[index,'CATEGORY'] = "Deletion"
                elif (row['STRAND'] == '-') & (row['POSITION'] > row['POSITION.1']):
                    sv.at[index,'CATEGORY'] = "Deletion"
                else:
                    sv.at[index,'CATEGORY'] = "Tandem_duplication"
            elif np.abs(row['POSITION'] - row['POSITION.1']) < 5e3:
                sv.at[index,'CATEGORY'] = "Inversion"
            else:
                sv.at[index,'CATEGORY'] = "Intrachromosomal_translocation"
        else:
            sv.at[index,'CATEGORY']  = "Interchromosomal_translocation"
    return(sv)

def subset_svs_to_panel(sv, PATH_bed):
    """
    Given an sv file and a path to a bed file, only keep the genes that are affected by the SV.
    """
    bed = pd.read_csv(PATH_bed, "\t", header = None)
    bed.columns = ["CHROM", "START", "END", "GENE"]
    bed_genes = bed["GENE"].unique()
    rows_keep = []
    for idx, row in sv.iterrows():
        if (not pd.isnull(row["NEARBY FEATURES"])) and (not pd.isnull(row["NEARBY FEATURES.1"])): # if the there are nearby features and the column isn't empty
            nearby_features = row["NEARBY FEATURES"].split(" ") + row["NEARBY FEATURES.1"].split(" ")
            for i, feature in enumerate(nearby_features):
                if feature == "AR": 
                    rows_keep.append(row)
                elif (feature in bed_genes) and (nearby_features[i + 1] == "(0),"): 
                    rows_keep.append(row)
    return pd.DataFrame(rows_keep).reset_index(drop = True)

def filter_genes(row, genes):
    nearby_features = row['NEARBY FEATURES'].split(', ')
    nearby_features_1 = row['NEARBY FEATURES.1'].split(', ')
    
    filtered_nearby_features = [gene for gene in nearby_features if gene.split(' ')[0] in genes]
    filtered_nearby_features_1 = [gene for gene in nearby_features_1 if gene.split(' ')[0] in genes]
    
    row['NEARBY FEATURES'] = ', '.join(filtered_nearby_features)
    row['NEARBY FEATURES.1'] = ', '.join(filtered_nearby_features_1)
    
    return row

prostate_patients = pd.DataFrame({"PATIENT": open(PATH_prostate_patients, "r").read().strip().split(" "), "DISEASE": "PROSTATE"})
bladder_patients = pd.DataFrame({"PATIENT": open(PATH_bladder_patients, "r").read().strip().split(" "), "DISEASE": "BLADDER"})
disease_patients = pd.concat([prostate_patients, bladder_patients], ignore_index = True)

prostate_genes = pd.read_csv(PATH_bed_prostate, "\t", header=None).rename(columns={3: 'GENE'})['GENE'].unique().tolist()
bladder_genes = pd.read_csv(PATH_bed_bladder, "\t", header=None).rename(columns={3: 'GENE'})['GENE'].unique().tolist()
genes_set = set(prostate_genes + bladder_genes)

sv_list_panel = []
sv_files = sorted([file for file in os.listdir(DIR_SVs) if file.endswith(".sv")])

for file in sv_files: 
    print(file)
    if not file == "annotations.txt":
        sv = pd.read_csv(os.path.join(DIR_SVs, file), "\t")
        sv_list_test.append(sv)
        if not sv.empty:
            sv = annotate_sv_type(sv)
            sv["SAMPLE"] = os.path.basename(file).replace(".sv", "") # add the sample name
            sv["PATIENT"] = sv["SAMPLE"].apply(lambda x: re.split("-cfDNA|_cfDNA|-WBC|_WBC", x)[0]) # get patient ID from sample name
            sv = sv.merge(disease_patients, how = "left") # indicate diagnosis
            sv["READ SUPPORT"] = sv['SUPPORTING READS'].apply(lambda x: x.count(";") + 1) # GET READ SUPPORT
            sv = sv[['PATIENT', 'SAMPLE', 'DISEASE'] + list(sv.columns.difference(['PATIENT', 'SAMPLE', 'DISEASE']))] 
            sv_list_all.append(sv) 
            if sv["DISEASE"].unique() == "PROSTATE":
                sv = subset_svs_to_panel(sv, PATH_bed_prostate)
            else: 
                sv = subset_svs_to_panel(sv, PATH_bed_bladder)
            if not sv.empty:
                sv[['NEARBY FEATURES', 'NEARBY FEATURES.1']] = sv[['NEARBY FEATURES', 'NEARBY FEATURES.1']].apply(filter_genes, axis=1, genes = genes_set)
                sv_list_panel.append(sv)

sv_panel = pd.concat(sv_list_panel)
annots = pd.read_csv("/home/vpc/datasets/plasma_chipseq/pipeline/analysis/rearrangements_annotation.tsv", "\t")
sv_panel = sv_panel.merge(annots, "left") # add the annotations in
sv_panel.to_csv("/home/vpc/datasets/plasma_chipseq/pipeline/analysis/Structural_variants_panel.csv", index = False)


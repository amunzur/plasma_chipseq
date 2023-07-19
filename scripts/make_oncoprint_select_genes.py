import os 
import pandas as pd
import numpy as np
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

DIR_project="/home/vpc/datasets/plasma_chipseq/pipeline"
PATH_cn=os.path.join(DIR_project, "analysis", "CN.csv")
PATH_mutations=os.path.join(DIR_project, "analysis", "mutations.csv")
PATH_ct=os.path.join(DIR_project, "analysis", "ct.csv")
PATH_figure_png=os.path.join("/home/vpc/asli/plasma_chip_seq/figures", "oncoprint_8_genes.png")
PATH_figure_pdf=os.path.join("/home/vpc/asli/plasma_chip_seq/figures", "oncoprint_8_genes.pdf")


# COLORS FOR MUTATIONS and CN 
mut_dict = {
    "Frameshift": '#FFC907', 
    "Missense": '#79B443', 
    "Stopgain": '#BD4398', 
    "Non-frameshift": '#a9a9a9', 
    "Promoter": "Red", 
    "Splice": "Black"}

mut_dict_shape = {
    "Altering": 's', 
    "Germline": '^'}

cn_dict_color = {
    "Deletion" : '#3F60AC', 
    "Deep_deletion" : '#9CC5E9',
    "Gain" : '#F59496',
    "Amplification" : '#EE2D24',
    "Neutral" : '#f2f2f2'
}

gene_groups = {
    'DNA Repair Genes': ['ARID1A', 'MSH2', 'MSH6', 'FANCD2', 'MLH1', 'CHD1', 'PMS2', 'RAD51D', 'BRCA2', 'BRCA1', 'ERCC2', 'BRIP1', 'PARP1', 'BARD1'],
    'Oncogenes and Tumor Suppressor Genes': ['MYC', 'PIK3R1', 'MET', 'BRAF', 'NKX3-1', 'PPP2R2A', 'NCOA2', 'CDKN2A', 
                                             'PTEN', 'CCND1', 'ZBTB16', 'KRAS', 'RB1', 'FOXA1', 
                                             'TP53', 'CDK12', 'SPOP', 'ASXL1', 'CHEK1', 'APC', 'MDM4'],
    'Transcription Factors and Chromatin Regulators': ['FOXP1', 'RYBP', 'ZFHX3'],
    'Cell Cycle Regulators': ['CDK6', 'CDKN1B', 'CDK4'],
    'Signaling Pathway Genes': ['PIK3CB', 'PIK3CA', 'KRAS', 'AKT1'],
    'Ubiquitin Ligases and E3 Ligases': ['MDM2', 'FANCL', 'FBXW7', 'RNF43'],
    'DNA Damage Response and Checkpoint Genes': ['ATM', 'ATR', 'CHEK1', 'CHEK2'],
    'Epigenetic Regulators': ['KMT2C', 'KMT2D', 'KDM6A', 'SMARCA1', 'TET2', 'DNMT3A'],
    'Hormone Receptor Genes': ['AR'],
    'DNA Polymerase Genes': ['POLE'],
    'Receptor Tyrosine Kinase Genes': ['NTRK1', 'NTRK3'],
    'Other Genes': ['CTNNB1', 'IDH1', 'FANCC', 'TSC1', 'FOLH1', 'RAD51B', 'TP53BP1', 'PALB2', 'FANCA', 'RAD51C']
}

cn=pd.read_csv(PATH_cn)
mut=pd.read_csv(PATH_mutations)
ct=pd.read_csv(PATH_ct)

# to make for select genes
mut = mut[mut['Gene'].isin(["AR", "TP53", "PTEN", "RB1", "SPOP", "FOXA1", "BRCA2", "CDK12"])]
cn = cn[cn['Gene'].isin(["AR", "TP53", "PTEN", "RB1", "SPOP", "FOXA1", "BRCA2", "CDK12"])]


# COPY NUMBER
cn['Pathway'] = cn['Gene'].map({gene: group for group, genes in gene_groups.items() for gene in genes})
cn = cn.sort_values(by = "Pathway").reset_index(drop=True)
cn["Genes_enumerated"] = cn['Gene'].map({sample: i / (len(cn['Gene'].unique()) - 1) for i, sample in enumerate(cn['Gene'].unique())})
cn['CN_Color'] = cn["CN event"].map(cn_dict_color)

# NON SILENT MUTATIONS
nonsilent=mut[mut["Characteristic"] != "Silent"]
nonsilent['Mutation_order'] = nonsilent.groupby(['Sample', 'Gene']).cumcount() + 1
nonsilent['TMB_nonsilent'] = nonsilent.loc[nonsilent['Characteristic'] != 'Germline'].groupby('Sample')['Sample'].transform('count') # not including germline muts in TMB calculation
nonsilent['Total_mutations_nonsilent_per_sample_per_gene'] = nonsilent.groupby(['Sample', 'Gene'])['Gene'].transform('count') # this is to help plotting

silent=mut[mut["Characteristic"] == "Silent"]
silent['TMB_silent'] = silent.groupby('Sample')['Sample'].transform('count')
silent = silent[["Sample", "Patient", "TMB_silent"]]
nonsilent = nonsilent.merge(silent, "left")
nonsilent['TMB_silent'] = nonsilent['TMB_silent'].fillna(0) # For the samples that have no silent mutations
nonsilent['TMB_nonsilent'] = nonsilent['TMB_nonsilent'].fillna(0) # For the samples that have no altering mutations
nonsilent["Mutation_color"] = nonsilent["Effect"].map(mut_dict)
nonsilent["Shape"] = nonsilent["Characteristic"].map(mut_dict_shape)

combined = cn.merge(nonsilent, "left")
combined = combined.merge(ct, "left")
combined = combined.sort_values(by = "ct percentage", ascending = False).reset_index(drop = True)
combined["Samples_enumerated"] = combined['Sample'].map({sample: i / (len(combined['Sample'].unique()) - 1) for i, sample in enumerate(combined['Sample'].unique())})

# per gene information
genes = nonsilent[["Gene", "Gene alteration percentage", "Effect"]].drop_duplicates()
genes = genes.merge(combined[["Gene", "Genes_enumerated"]], "left")
genes["Color"] = genes["Effect"].map(mut_dict)
genes = genes.drop_duplicates().reset_index(drop = True)
genes["Gene alteration percentage"] = genes["Gene alteration percentage"]

# ONCOPRINT
fig_width = 12
fig_height = 15

fig = plt.figure(figsize = (fig_width, fig_height))
gs = mpl.gridspec.GridSpec(ncols = 2, nrows = 3, height_ratios = [0.13, 0.13, 1], width_ratios = [1, 0.15], hspace = 0.05, wspace =0)

#Add and space out subplots
main = fig.add_subplot(gs[2,0])
tc = fig.add_subplot(gs[0,0], sharex=main)
tmb = fig.add_subplot(gs[1,0], sharex=main)
alt_percentage = fig.add_subplot(gs[2,1], sharey=main)
fig.subplots_adjust(bottom=0.3)  # Increase the bottom padding

# X and Y limits 
main.set_xlim(-0.01, 1.01)
main.set_ylim(-0.01, 1.1)

bar_height = main.get_position().height/combined["Gene"].unique().size + 0.006
bar_width = main.get_position().width/combined["Sample"].unique().size + 0.0058

# ct percentage barplot
tc.bar(x = combined["Samples_enumerated"], height = combined["ct percentage"], width=bar_width, color = "black", zorder = 10)

# tmb bar plot
tmb.bar(x = combined["Samples_enumerated"], height = combined["TMB_nonsilent"], width=bar_width, color = "skyblue", zorder = 10)
tmb.bar(x = combined["Samples_enumerated"], height = combined["TMB_silent"], bottom = combined["TMB_nonsilent"], width=bar_width, color = "mediumpurple", zorder = 10)

# Plot copy numbers 
for index, row in combined.iterrows():
    main.bar(x = row["Samples_enumerated"], height=bar_height, width=bar_width, bottom = row["Genes_enumerated"], color = row["CN_Color"])
    if row["Characteristic"] == "Altering" or row["Characteristic"] == "Germline":
        if row["Total_mutations_nonsilent_per_sample_per_gene"] == 1: # if the patient has only 1 mutation in a given gene
            main.scatter(x = row["Samples_enumerated"], y = row["Genes_enumerated"] + bar_height/2 , color = row["Mutation_color"], marker = row["Shape"], zorder = 90, s = 35, lw = 0)
        if row["Total_mutations_nonsilent_per_sample_per_gene"] == 2:
            if row["Mutation_order"] == 1: # first mutation in the pair
                main.scatter(x = row["Samples_enumerated"], y = row["Genes_enumerated"] + 0.01 + bar_height/2 , color = row["Mutation_color"], marker = row["Shape"], zorder = 90, s = 35, lw = 0)
            if row["Mutation_order"] == 2: # second mutation in the pair
                main.scatter(x = row["Samples_enumerated"], y = row["Genes_enumerated"] - 0.01 + bar_height/2 , color = row["Mutation_color"], marker = row["Shape"], zorder = 90, s = 35, lw = 0)
        if row["Total_mutations_nonsilent_per_sample_per_gene"] >= 3:
            main.scatter(x = row["Samples_enumerated"], y = row["Genes_enumerated"] + bar_height/2 , color = "red", marker = "*", zorder = 90, s = 35, lw = 0)

# PLOT THE GENE ALTERATION PERCENTAGE AS A STACKED HORIZONTAL BAR PLOT
grouped_genes = genes.groupby(by = "Gene")
for gene_name, gene_info in grouped_genes:
    gene_info = gene_info.sort_values(by="Gene alteration percentage").reset_index(drop = True)
    if gene_info.shape[0] == 1:  # if there is only 1 mutation type in this gene
        alt_percentage.barh(y = gene_info["Genes_enumerated"] + bar_height/2, left = 0, width = gene_info["Gene alteration percentage"], align="center", color=gene_info["Color"], height = bar_height, zorder = 10)
    else:
        left_position = 0  # Initialize the left position for the first bar
        for idx, row in gene_info.iterrows():
            alt_percentage.barh(y = row["Genes_enumerated"] + bar_height/2, left = left_position, width = row["Gene alteration percentage"], align="center", color=row["Color"], height = bar_height, zorder = 10)
            left_position += row["Gene alteration percentage"]


# X and Y ticks
main.set_yticks(cn["Genes_enumerated"].unique().tolist() + bar_height/2)
main.set_yticklabels(cn["Gene"].unique().tolist())
main.set_xticks(combined["Samples_enumerated"].unique().tolist()) # so they align to the center of the squares
main.set_xticklabels(combined["Sample"].unique().tolist(), rotation = 90)

alt_percentage.set_xticks([0, 25, 50, 75, 100])
alt_percentage.set_xticklabels(["0", "25", "50", "75", "100"])
alt_percentage.tick_params(axis='y', left=False, labelleft=False)

main.tick_params(axis='x', bottom=False, pad=-2)
main.tick_params(axis='y', left=False, pad=-2)

tc.tick_params(axis='x', bottom=False, labelbottom=False)
tmb.tick_params(axis='x', bottom=False, labelbottom=False)

tmb.set_yticks([0, 10, 20, 30])
tmb.set_yticklabels(["0", "10", "20", "30"])


# Spines
main.spines["top"].set_visible(False)
main.spines["bottom"].set_visible(False)
main.spines["right"].set_visible(False)
main.spines["left"].set_visible(False)

tmb.spines["top"].set_visible(False)
tmb.spines["right"].set_visible(False)

tc.spines["top"].set_visible(False)
tc.spines["right"].set_visible(False)

alt_percentage.spines["top"].set_visible(False)
alt_percentage.spines["right"].set_visible(False)
alt_percentage.spines["left"].set_visible(False)

# X and Y labels
tc.set_ylabel("ctDNA (%)")
tmb.set_ylabel("TMB")
main.set_ylabel("Genes")
main.set_xlabel("Samples")
alt_percentage.set_xlabel("Gene \nalteration (%)")

# Dashed lines across ticks
for tick in tc.get_yticks():
    tc.axhline(tick, color='lightgray', linestyle='dashed', linewidth=0.5)

for tick in tmb.get_yticks():
    tmb.axhline(tick, color='lightgray', linestyle='dashed', linewidth=0.5)

for tick in alt_percentage.get_xticks():
    alt_percentage.axvline(tick, color='lightgray', linestyle='dashed', linewidth=0.5)

fig.tight_layout()
fig.savefig(PATH_figure_png)
fig.savefig(PATH_figure_pdf)
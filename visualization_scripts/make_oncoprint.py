import os 
import pandas as pd
import numpy as np
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

DIR_project="/home/vpc/datasets/plasma_chipseq/pipeline"
PATH_cn=os.path.join(DIR_project, "analysis", "CN.csv")
PATH_mutations=os.path.join(DIR_project, "analysis", "all_mutations_curated.tsv")
PATH_ct=os.path.join(DIR_project, "analysis", "ct_fractions.tsv")
PATH_figure_png=os.path.join(DIR_project, "figures", "oncoprint_8_genes.png")
PATH_figure_pdf=os.path.join(DIR_project, "figures", "oncoprint_8_genes.pdf")
PATH_prostate_pairs=os.path.join(DIR_project, "sample_lists", "prostate_paired.tsv")
PATH_naming_dictionary=os.path.join(DIR_project, "sample_lists", "naming_dictionary.tsv")

PATH_mutations_silent_altering_PROCESSED=os.path.join(DIR_project, "analysis", "silent_and_altering_mutations.csv")
PATH_mutations_germline_PROCESSED=os.path.join(DIR_project, "analysis", "germline_mutations.csv")

# Whole genome duplication status
sheet_id = "1W3Xw_t56dxOY3iHJ3FXW_L2wT8DVQ8GoLbapwhzRjFo"
sheet_name = "WGD_[ARCHIVED]"
URL_CN_ctDNA = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"

# Imaging data
sheet_id = "1qyar6XDXhH7SBJ3uAVTCJlFXsXrav4DIoadIi_n8kfo"
sheet_name = "CT"
URL_imaging = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"

sheet_name_scantiming = "Scan_timing"
URL_scan_timing = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name_scantiming}"

# clinical data
sheet_id = "1W3Xw_t56dxOY3iHJ3FXW_L2wT8DVQ8GoLbapwhzRjFo"
sheet_name = "Clinical_Data"
URL_clinical = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"

# SV data
sheet_id = "1W3Xw_t56dxOY3iHJ3FXW_L2wT8DVQ8GoLbapwhzRjFo"
sheet_name = "Structural_variants"
URL_SV = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"


# COLORS FOR MUTATIONS and CN 
mut_dict = {
    "Frameshift": '#FFC907', 
    "Missense": '#79B443', 
    "Stopgain": '#BD4398', 
    "Non-frameshift": '#a9a9a9', 
    "Splice": "Black"}

mut_dict_shape = {
    "Altering": 's', 
    "Germline": '*'}

mut_dict_size = {
    "Altering": 25, 
    "Germline": 30}

mut_dict_zorder = {
    "Altering": 30, 
    "Germline": 40}

cn_dict_color = {
    "Deep_deletion" : '#3F60AC', 
    "Deletion" : '#9CC5E9',
    "Gain" : '#F59496',
    "Amplification" : '#EE2D24',
    "Neutral" : '#f2f2f2'
}

sv_color = {
    "Deletion" : '#33b8af', 
    "Tandem_duplication" : '#ffd681',
    "Interchromosomal_translocation" : '#916e53',
    "Intrachromosomal_translocation" : '#c87fc6',
    "Inversion" : '#61a165'
}

ar_sv_color = {
    "LBD truncating": "black",
    "Enhancer": "#585858",
    "Neighborhood": "#9a9a9a", 
    "Other": "#dadada"    
}

wgd_color = {
    True: "Black", 
    False: "Gray",
    np.nan: "Whitesmoke"
}

imaging_dict = {
    "Yes": True,
    "No": False
}


# LOADING DFs
wgd = pd.read_csv(URL_CN_ctDNA)[["Takeda_Patient", "Patient", "WGD"]]
wgd["WGD_colors"] = wgd["WGD"].map(wgd_color)

cn = pd.read_csv(PATH_cn)
cn = cn[cn['Gene'].isin(["TP53", "SPOP", "CDK12", "FOXA1", "PTEN", "RB1", "AR", "BRCA2"])].reset_index(drop = True)
cn = cn.merge(wgd, "left")

mut=pd.read_csv(PATH_mutations, "\t")
mut = mut[mut['Gene'].isin(["TP53", "SPOP", "CDK12", "FOXA1", "PTEN", "RB1", "AR", "BRCA2"])].reset_index(drop = True)

ct=pd.read_csv(PATH_ct, "\t")[["Patient", "Sample", "Final_ct_fraction"]]
ct.columns = ["Patient", "Sample", "ct percentage"]
ct["ct percentage"] = ct["ct percentage"]*100

pairs=pd.read_csv(PATH_prostate_pairs, "\t")
pairs_dict = pairs.set_index('WBC')['cfDNA'].to_dict()

# COPY NUMBER
cn["Genes_enumerated"] = cn['Gene'].map({sample: i / (len(cn['Gene'].unique()) - 1) for i, sample in enumerate(cn['Gene'].unique())})
cn['CN_Color'] = cn["CN event"].map(cn_dict_color)

# NON SILENT MUTATIONS
nonsilent = mut[mut["Characteristic"] != "Silent"]
keys = nonsilent[nonsilent["Characteristic"] == "Germline"]["Sample"]
nonsilent.loc[nonsilent["Characteristic"] == "Germline", "Sample"] = [pairs_dict[key] for key in keys] # replace the wbc samples with their cfDNA counterpart
nonsilent = nonsilent[nonsilent["Gene"] != "TERT"]

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
nonsilent["Size"] = nonsilent["Characteristic"].map(mut_dict_size)
nonsilent["Zorder"] = nonsilent["Characteristic"].map(mut_dict_zorder)

combined = cn.merge(nonsilent, on=['Patient', 'Sample', 'Gene'], how='left')
combined = combined.merge(ct, "left")
combined = combined.sort_values(by = "ct percentage", ascending = False).reset_index(drop = True)
combined["Samples_enumerated"] = combined['Sample'].map({sample: i / (len(combined['Sample'].unique()) - 1) for i, sample in enumerate(combined['Sample'].unique())})
combined = combined.merge(pd.read_csv(PATH_naming_dict, "\t")[["Patient", "Takeda_Patient", "Sample", "Timepoint"]])
combined["X ticks"] = combined["Patient"].str.replace("GU-", "") + " " + combined["Takeda_Patient"] + " " + combined["Timepoint"].astype(str)

# per gene information
nonsilent["Gene alteration percentage"] = (nonsilent.groupby(["Gene", "Effect"])["Sample"].transform("count")/cn["Sample"].unique().size)*100
genes = nonsilent[["Gene", "Gene alteration percentage", "Effect"]].drop_duplicates()
genes = genes.merge(combined[["Gene", "Genes_enumerated"]], "left")
genes["Color"] = genes["Effect"].map(mut_dict)
genes = genes.drop_duplicates().reset_index(drop = True)
genes["Gene alteration percentage"] = genes["Gene alteration percentage"]

# germline mutations are tackled separately
germline = nonsilent[nonsilent["Characteristic"] == "Germline"]
germline = germline.merge(combined[["Patient", "Gene", "Genes_enumerated", "Samples_enumerated"]], on = ["Patient", "Gene"], how = "left")

#===========================================================================================================
#IMAGING DATA
naming_dict = pd.read_csv(PATH_naming_dictionary, "\t")[["Patient", "Diagnosis", "Sample"]]
naming_dict["Collection_date"] = naming_dict['Sample'].str.extract(r'(\d{4}[A-Za-z]{3}\d{2})')

imaging = pd.read_csv(URL_imaging)[["GUBB ID", "Date", "Bone mets", "Any lytic component", "Liver mets", "Lung mets", "Other \nvisceral mets"]]
imaging = imaging.rename(columns={'GUBB ID': 'Patient', 'Other \nvisceral mets': 'Other visceral mets', 'Date': 'CT_Date'})
imaging["Patient"] = "GU-" + imaging["Patient"]

scantiming = pd.read_csv(URL_scan_timing)[["GUBB_ID", "Collection_date", "Nearest CT"]]
scantiming["GUBB_ID"] = "GU-" + scantiming["GUBB_ID"]
scantiming = scantiming.rename(columns={'GUBB_ID': 'Patient', 'Nearest CT':'CT_Date'})
imaging = imaging.merge(scantiming, "left") # this merge adds the sample collection date to the imaging sheet
imaging["Collection_date"] = pd.to_datetime(imaging["Collection_date"]).dt.strftime('%Y%b%d') # putting it in the same format as the naming dict

imaging = imaging.merge(naming_dict, "left") # adding the sample id to the imaging data
for column in ["Bone mets", "Any lytic component", "Lung mets", "Liver mets"]: # Turn to boolean
    imaging[column] = imaging[column].map(imaging_dict)

imaging["Other visceral mets"] = imaging["Other visceral mets"].apply(lambda x: False if pd.isnull(x) else True)

imaging = imaging.merge(combined[["Samples_enumerated", "Sample"]].drop_duplicates().reset_index(drop= True), "left") # add the samples enumarated column to the imaging data
imaging = imaging[~pd.isnull(imaging["Samples_enumerated"])] # this helps drop samples that wont be included in the analysis
imaging = imaging[["Samples_enumerated", "Collection_date", "Bone mets", "Any lytic component", "Liver mets", "Lung mets", "Other visceral mets"]] # only keeping what is necessary
imaging = pd.melt(imaging, id_vars=["Samples_enumerated", "Collection_date"], var_name='Mets', value_name='Value') # reformatting will help with formatting
imaging["Mets enumerated"] = imaging["Mets"].map({sample: i / (len(imaging['Mets'].unique()) - 1) for i, sample in enumerate(imaging['Mets'].unique())})
imaging["Mets color"] = imaging["Value"].map(wgd_color)

#===========================================================================================================
#CLINICAL DATA (PSA, LDH, ALP)
clinical = pd.read_csv(URL_clinical)
clinical = clinical[clinical["Exclude"] != "Yes"].reset_index(drop = True)[['GUBB ID', 'Collection date', 'PSA', 'LDH', 'ALP']]

colorbars = []  # List to store the colorbar legend objects
for variable in ['PSA', 'LDH', 'ALP']:
    clinical[variable] = clinical[variable].replace(0, 1)  # Before log transformation, replace 0 with 1
    clinical[variable] = np.log10(clinical[variable])
    cmap = plt.cm.get_cmap("coolwarm")
    normalize = plt.Normalize(clinical[variable].min(), clinical[variable].max())
    clinical[variable + "_color"] = clinical[variable].map(lambda x: cmap(normalize(x)))
    # Create the colorbar
    # sm = plt.cm.ScalarMappable(cmap=cmap, norm=normalize)
    # colorbar = plt.colorbar(sm, label='Colorbar Label')
    # colorbars.append(colorbar)



clinical["GUBB ID"] = "GU-" + clinical["GUBB ID"]
clinical = clinical.rename(columns={'GUBB ID': 'Patient', "Collection date": "Collection_date"})
clinical["Collection_date"] = pd.to_datetime(clinical["Collection_date"]).dt.strftime('%Y%b%d') # putting it in the same format as the naming dict
clinical = clinical.merge(imaging[["Samples_enumerated", "Collection_date"]], "left")
clinical = clinical[["Samples_enumerated", "Collection_date", "PSA_color", "ALP_color", "LDH_color"]]
clinical = pd.melt(clinical, id_vars=["Collection_date", "Samples_enumerated"], var_name='Color', value_name='Value')
clinical["Color"] = clinical["Color"].str.replace("_color", "") # remove the color prefix
clinical = clinical.rename(columns={"Color": "Blood counts"})
clinical["Clinical_enumerated"] = clinical["Blood counts"].map({sample: i / (len(clinical['Blood counts'].unique()) - 1) for i, sample in enumerate(clinical['Blood counts'].unique())})
#===========================================================================================================

#===========================================================================================================
# STRUCTURAL VARIANTS
# to plot a bar chart of SV frequency
sv_all = pd.read_csv(URL_SV)
sv_all = sv_all.merge(combined[["Sample", "Samples_enumerated"]].drop_duplicates().reset_index(drop = True), how = "left", left_on = "SAMPLE", right_on = "Sample")
sv_all["SV color"] = sv_all["CATEGORY"].map(sv_color)
sv_all = sv_all.loc[(sv_all['PATHOGENIC'] != False)].reset_index(drop = True)
region_counts_ALL_SVs = sv_all.groupby(['Samples_enumerated', 'CATEGORY']).size().unstack(fill_value=0)

genes = ["TP53", "RB1", "PTEN"]
pattern = r"\b(" + "|".join(genes) + r")\b"
sv_select_genes = sv_all[sv_all["NEARBY FEATURES"].str.contains(pattern, flags=re.IGNORECASE, na=False) | sv_all["NEARBY FEATURES.1"].str.contains(pattern, flags=re.IGNORECASE, na=False)].reset_index(drop = True)
sv_select_genes = sv_select_genes.loc[(sv_select_genes['PATHOGENIC'] != False)].reset_index(drop = True)
sv_select_genes['GENES'] = sv_select_genes['NEARBY FEATURES'].combine_first(sv_select_genes['NEARBY FEATURES.1']).str.extract(pattern) + " " + "SVs"
sv_select_genes['Genes_enumerated'] = sv_select_genes['GENES'].map({sample: i / (len(sv_select_genes['GENES'].unique()) - 1) for i, sample in enumerate(sv_select_genes['GENES'].unique())})
#===========================================================================================================

#===========================================================================================================
# AR STRUCTURAL VARIANTS
sv_ar = sv_all.loc[(sv_all['NEARBY FEATURES'].str.contains("AR")) | (sv_all['NEARBY FEATURES.1'].str.contains("AR"))].reset_index(drop = True)
sv_ar["REGION"] = "Dummy"
sv_ar.loc[(sv_ar['POSITION'] >= 66890000) & (sv_ar['POSITION'] <= 66920000), 'REGION'] = 'Enhancer'
sv_ar.loc[(sv_ar['POSITION.1'] >= 66890000) & (sv_ar['POSITION.1'] <= 66920000), 'REGION'] = 'Enhancer'
sv_ar.loc[sv_ar['ANNOTATION'].str.contains("LBD truncating"), 'REGION'] = 'LBD truncating' # the category becomes LBD truncating if the annotation states that it is
sv_ar.loc[sv_ar['ANNOTATION'].str.contains("AR neighborhood rearrangement"), 'REGION'] = 'Neighborhood' # the category becomes LBD truncating if the annotation states that it is
sv_ar['REGION'] = sv_ar['REGION'].replace('Dummy', 'Other')
region_counts = sv_ar.groupby(['Samples_enumerated', 'REGION']).size().unstack(fill_value=0)
#===========================================================================================================






#===========================================================================================================

# ONCOPRINT
fig_width = 15
fig_height = 17

fig = plt.figure(figsize = (fig_width, fig_height))
gs = mpl.gridspec.GridSpec(ncols = 1, nrows = 8, height_ratios = [0.17, 0.17, 0.02, 0.15, 0.10, 0.10, 0.10, 1], hspace = 0.08)

#Add and space out subplots
tc = fig.add_subplot(gs[0], sharex=main) # tumor content
tmb = fig.add_subplot(gs[1], sharex=main) # tumor mutation burden
wgd = fig.add_subplot(gs[2], sharex=main) # whole genome duplication status
scan = fig.add_subplot(gs[3], sharex=main) # imaging (scan) results
clinic = fig.add_subplot(gs[4], sharex=main) # clinical data
svar_select_genes = fig.add_subplot(gs[5], sharex=main) # TP53, RB1, PTEN
svar_AR = fig.add_subplot(gs[7], sharex=main) # barplot showing the frequency of each type of AR SV
main = fig.add_subplot(gs[8]) # gene CN and mutation status

fig.subplots_adjust(bottom=0.3)  # Increase the bottom padding
fig.subplots_adjust(left=0.2)  # Increase the bottom padding

# X and Y limits 
main.set_xlim(-0.01, 1.01)
main.set_ylim(-0.01, 1.1)

# bar_height = main.get_position().height/combined["Gene"].unique().size + 0.0075
bar_height = main.get_position().height/combined["Gene"].unique().size + 0.05
bar_width = main.get_position().width/combined["Sample"].unique().size + 0.0058

#===========================================================================================================
# PLOTTING
tc.bar(x = combined["Samples_enumerated"], height = combined["ct percentage"], width=bar_width, color = "black", zorder = 10)
tmb.bar(x = combined["Samples_enumerated"], height = combined["TMB_nonsilent"], width=bar_width, color = "limegreen", zorder = 10)
tmb.bar(x = combined["Samples_enumerated"], height = combined["TMB_silent"], bottom = combined["TMB_nonsilent"], width=bar_width, color = "mediumpurple", zorder = 10)
wgd.bar(x = combined["Samples_enumerated"], height = 1, width=bar_width, color = combined["WGD_colors"])
scan.bar(x=imaging["Samples_enumerated"], bottom=imaging["Mets enumerated"], height=0.2, width=bar_width, color=imaging["Mets color"])
clinic.bar(x=clinical["Samples_enumerated"], bottom=clinical["Clinical_enumerated"], height=0.45, width=bar_width, color=clinical["Value"])
svar_select_genes.bar(x=sv_select_genes["Samples_enumerated"], bottom=sv_select_genes["Genes_enumerated"], height=0.45, width=bar_width, color="black")

bar_position = [0] * len(region_counts.index)
for region in region_counts.columns:
    svar_AR.bar(x=region_counts.index, height=region_counts[region], bottom=bar_position, label=region, width=bar_width, color=ar_sv_color.get(region, 'gray'))  # Set the color for the region from the dictionary
    bar_position += region_counts[region] # Update the starting position for the next region

#===========================================================================================================


# Plot copy numbers and somatic mutations
for index, row in combined.iterrows():
    main.bar(x = row["Samples_enumerated"], height=bar_height, width=bar_width, bottom = row["Genes_enumerated"], color = row["CN_Color"])
    if row["Characteristic"] == "Altering" or row["Characteristic"] == "Germline":
        if row["Total_mutations_nonsilent_per_sample_per_gene"] == 1: # if the patient has only 1 mutation in a given gene
            main.scatter(x = row["Samples_enumerated"], y = row["Genes_enumerated"] + bar_height/2 , color = row["Mutation_color"], marker = row["Shape"], zorder = row["Zorder"], s = row["Size"], lw = 0)
        if row["Total_mutations_nonsilent_per_sample_per_gene"] == 2:
            if row["Mutation_order"] == 1: # first mutation in the pair
                main.scatter(x = row["Samples_enumerated"], y = row["Genes_enumerated"] + 0.01 + bar_height/2 , color = row["Mutation_color"], marker = row["Shape"], zorder = row["Zorder"], s = row["Size"], lw = 0)
            if row["Mutation_order"] == 2: # second mutation in the pair
                main.scatter(x = row["Samples_enumerated"], y = row["Genes_enumerated"] - 0.01 + bar_height/2 , color = row["Mutation_color"], marker = row["Shape"], zorder = row["Zorder"], s = row["Size"], lw = 0)
        if row["Total_mutations_nonsilent_per_sample_per_gene"] > 2:
            main.scatter(x=row["Samples_enumerated"], y=row["Genes_enumerated"] + bar_height/2, marker="^", zorder=row["Zorder"], s=row["Size"], lw=1, facecolor="white", edgecolor="black")

# X and Y ticks
main.set_yticks(cn["Genes_enumerated"].unique().tolist() + bar_height/2)
main.set_yticklabels(cn["Gene"].unique().tolist())
main.set_xticks(combined["Samples_enumerated"].unique().tolist()) # so they align to the center of the squares
main.set_xticklabels(combined["X ticks"].unique().tolist(), rotation = 90)

main.tick_params(axis='x', bottom=False, pad=-2)
main.tick_params(axis='y', left=False, pad=-2)

wgd.set_yticks([0.55])
wgd.set_yticklabels(["WGD"])

svar_AR.set_ylim(0, 20)
svar_AR.set_yticks([0, 10, 20])
svar_AR.set_yticklabels(["0", "10", "20"])
#===========================================================================================================
# y ticks for the scan axis
scan.set_yticklabels([]) # first remove the default
scan.set_yticks([])
scan.set_yticks(np.sort(imaging["Mets enumerated"].unique()) + 0.12) # add new ones
scan.set_xlim(-0.01, 1.01) # to remove the horizontal gap after setting y ticks
mets_dict = {mets: enumerated for mets, enumerated in zip(imaging["Mets"], imaging["Mets enumerated"])} # a dict matching the mets with the mets enumerated, helps label the y ticks
ytick_labels = []
for mets in scan.get_yticks() - 0.12:
    for key, value in mets_dict.items():
        if np.isclose(mets, value, atol=0.01):
            ytick_labels.append(key)
            break

scan.set_yticklabels(ytick_labels)
#===========================================================================================================

#===========================================================================================================
# ticks for the clinic axis
clinic.set_yticklabels([]) # first remove the default
clinic.set_yticks([])
clinic.set_yticks(np.sort(clinical["Clinical_enumerated"].unique()) + 0.2) # add new ones
clinic.set_xlim(-0.01, 1.01) # to remove the horizontal gap after setting y ticks
blood_counts_dict = {Blood_counts: enumerated for Blood_counts, enumerated in zip(clinical["Blood counts"], clinical["Clinical_enumerated"])} # a dict matching the mets with the mets enumerated, helps label the y ticks
ytick_labels = []
for y in clinic.get_yticks() - 0.2:
    for key, value in blood_counts_dict.items():
        if np.isclose(y, value, atol=0.01):
            ytick_labels.append(key)
            break

clinic.set_yticklabels(ytick_labels)
#===========================================================================================================

#===========================================================================================================
# ticks for the SVs in select genes axis (TP53, RB1, PTEN)
svar_select_genes.set_yticklabels([]) # first remove the default
svar_select_genes.set_yticks([])
svar_select_genes.set_yticks(np.sort(sv_select_genes["Genes_enumerated"].unique()) + 0.2) # add new ones
# svar_select_genes.set_xlim(-0.01, 1.01) # to remove the horizontal gap after setting y ticks
select_genes_dict = {SV_status: enumerated for SV_status, enumerated in zip(sv_select_genes["GENES"], sv_select_genes["Genes_enumerated"])} # a dict matching the mets with the mets enumerated, helps label the y ticks
ytick_labels = []
for y in svar_select_genes.get_yticks() - 0.2:
    for key, value in select_genes_dict.items():
        if np.isclose(y, value, atol=0.01):
            ytick_labels.append(key)
            break

svar_select_genes.set_yticklabels(ytick_labels)
#===========================================================================================================


tmb.set_yticks([0, 5, 10])
tmb.set_yticklabels(["0", "5", "10"])

for ax in [tc, tmb, wgd, scan, clinic, svar_AR, svar_select_genes]:
    ax.tick_params(axis='x', bottom=False, labelbottom=False)

for ax in [wgd, scan, clinic, svar_select_genes]:
    ax.tick_params(axis='y', left=False, labelbottom=False)

# this blank space above is to stay
#===========================================================================================================
# Spines
main.spines["top"].set_visible(False)
main.spines["bottom"].set_visible(False)
main.spines["right"].set_visible(False)
main.spines["left"].set_visible(False)

wgd.spines["top"].set_visible(False)
wgd.spines["bottom"].set_visible(False)
wgd.spines["right"].set_visible(False)
wgd.spines["left"].set_visible(False)

scan.spines["top"].set_visible(False)
scan.spines["bottom"].set_visible(False)
scan.spines["right"].set_visible(False)
scan.spines["left"].set_visible(False)

clinic.spines["top"].set_visible(False)
clinic.spines["bottom"].set_visible(False)
clinic.spines["right"].set_visible(False)
clinic.spines["left"].set_visible(False)

svar_select_genes.spines["top"].set_visible(False)
svar_select_genes.spines["bottom"].set_visible(False)
svar_select_genes.spines["right"].set_visible(False)
svar_select_genes.spines["left"].set_visible(False)

svar_AR.spines["top"].set_visible(False)
svar_AR.spines["right"].set_visible(False)

tmb.spines["top"].set_visible(False)
tmb.spines["right"].set_visible(False)

tc.spines["top"].set_visible(False)
tc.spines["right"].set_visible(False)
#===========================================================================================================

# X and Y labels
tc.set_ylabel("ctDNA (%)")
tmb.set_ylabel("TMB")
svar_AR.set_ylabel("AR SVs")
main.set_ylabel("Genes")
main.set_xlabel("Samples")

for ax in [main, svar_AR, svar_select_genes, tmb, tc]:
    ax.yaxis.set_label_coords(-0.1, 0.5)  # Move the label slightly to the left
    ax.yaxis.label.set_rotation(0)  

# Dashed lines across ticks
for ax in [tc, tmb, svar_AR]:
    for tick in ax.get_yticks():
        ax.axhline(tick, color='lightgray', linestyle='dashed', linewidth=0.5)

########################
# LEGEND
########################

mut_dict_shape_legend = {
    "Somatic": 's', 
    "Germline": '*', 
    ">2 mutations": "^"}

tmb_legend = {
    "Nonsilent": "limegreen",
    "Silent": "mediumpurple"
}

wgd_legend = {
    "WGD": "Black",
    "No WGD": "Gray", 
    "Unevaluable": "Whitesmoke"  
}

imaging_legend = {
    "Metastasis": "Black",
    "No metastasis": "Gray", 
    "No information": "Whitesmoke"  
}

sv_legend = {
    "Deletion" : '#33b8af', 
    "Tandem dup." : '#ffd681',
    "Interch. transl." : '#916e53',
    "Intrach. transl." : '#c87fc6',
    "Inversion" : '#61a165'
}

handles_cn = []
for key in cn_dict_color:
    handle = mpatches.Patch(color = cn_dict_color.get(key), label = key)
    handles_cn.append(handle)

handles_muts = []
for key in mut_dict:   
    handle = mpatches.Patch(color = mut_dict.get(key), label = key)
    handles_muts.append(handle)

handles_mut_shapes = []
for key in mut_dict_shape_legend:    
    handle = Line2D([0], [0], linestyle = "none", marker = mut_dict_shape_legend.get(key), label = key, markersize=10, markeredgecolor='black', markerfacecolor='white')
    handles_mut_shapes.append(handle)

handles_tmb = []
for key in tmb_legend:   
    handle = mpatches.Patch(color = tmb_legend.get(key), label = key)
    handles_tmb.append(handle)

handles_wgd = []
for key in wgd_legend:   
    handle = mpatches.Patch(color = wgd_legend.get(key), label = key)
    handles_wgd.append(handle)

handles_imaging = []
for key in imaging_legend:   
    handle = mpatches.Patch(color = imaging_legend.get(key), label = key)
    handles_imaging.append(handle)

handles_AR_SV = []
for key in ar_sv_color:
    handle = mpatches.Patch(color = ar_sv_color.get(key), label = key)
    handles_AR_SV.append(handle)

handles_SV = []
for key in sv_legend:
    handle = mpatches.Patch(color = sv_legend.get(key), label = key)
    handles_SV.append(handle)

legend1 = fig.legend(handles=handles_cn, bbox_to_anchor=(0.18, 0.17), frameon=False, title = "Copy number variants", title_fontsize = 12)
legend2 = fig.legend(handles=handles_mut_shapes, bbox_to_anchor=(0.28, 0.17), frameon=False, title = "Shapes", title_fontsize = 12)
legend3 = fig.legend(handles=handles_muts, bbox_to_anchor=(0.39, 0.17), frameon=False, title = "Mutations", title_fontsize = 12)
legend4 = fig.legend(handles=handles_tmb, bbox_to_anchor=(0.47, 0.17), frameon=False, title = "TMB colors", title_fontsize = 12)
legend5 = fig.legend(handles=handles_wgd, bbox_to_anchor=(0.57, 0.17), frameon=False, title = "WGD status", title_fontsize = 12)
legend6 = fig.legend(handles=handles_imaging, bbox_to_anchor=(0.69, 0.17), frameon=False, title = "Metastasis", title_fontsize = 12)
legend7 = fig.legend(handles=handles_AR_SV, bbox_to_anchor=(0.80, 0.17), frameon=False, title = "AR SVs", title_fontsize = 12)
legend8 = fig.legend(handles=handles_SV, bbox_to_anchor=(0.90, 0.17), frameon=False, title = "SVs", title_fontsize = 12)

# align the legend titlesshapes_dict = {}
legend1._legend_box.align = "left"
legend2._legend_box.align = "left"
legend3._legend_box.align = "left"
legend4._legend_box.align = "left"
legend5._legend_box.align = "left"
legend6._legend_box.align = "left"
legend7._legend_box.align = "left"
legend8._legend_box.align = "left"

fig.tight_layout()
fig.savefig(PATH_figure_png)
fig.savefig(PATH_figure_pdf)


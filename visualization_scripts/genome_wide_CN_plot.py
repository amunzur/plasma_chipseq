"""
This script make a genome wide CN status plot.
"""
import numpy as np 
import pandas as pd
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
import statsmodels.api as sm
import matplotlib.patches as mpatches

DIR_igv_tracks = "/home/vpc/datasets/plasma_chipseq/pipeline/igv_tracks_admixture_prostate/"
PATH_naming_dict = "/home/vpc/datasets/plasma_chipseq/pipeline/sample_lists/naming_dictionary.tsv"
PATH_chrom_sizes = "/home/vpc/asli/resources/chrom_sizes.tsv"

sheet_id = "1W3Xw_t56dxOY3iHJ3FXW_L2wT8DVQ8GoLbapwhzRjFo"
sheet_name = "WGD_[ARCHIVED]"
URL_CN_ctDNA = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"

sheet_id = "1W3Xw_t56dxOY3iHJ3FXW_L2wT8DVQ8GoLbapwhzRjFo"
sheet_name = "Tumor_fraction"
URL_TF = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"

tf = pd.read_csv(URL_TF)
tf = tf[tf["Final_ct_fraction"] > 0.20]

pts = pd.read_csv(PATH_naming_dict, delimiter="\t")
pts = pts[pts["Diagnosis"] == "Prostate"]
pts = pts[pts["Sample"].isin(tf["Sample"])] # limit this analysis to high TF patients

def bin_chrom(group, bin_size):
    """
    Given a df grouped by chroms, bin the chroms. This function is applied to each groups.
    """
    start = 0 
    end = group["START"].max()
    bin_edges = np.arange(start, end + bin_size, bin_size)
    group['Bin'] = pd.cut(df['START'], bins=bin_edges)
    return(group)

def determine_CN_events(patient, median_values_df, wgd):
    """
    Given a sample name, the binned logratio df and the wgd sheet that indicates if a sample has wgd, this function annotates the CN events 
    in the log ratio as amplification and deletion.
    """
    if wgd[wgd["Patient"] == patient]["WGD"]: 
        min_threshold_amp = 2 
        min_threshold_del = -1
        filtered_bins = median_values.loc[median_values > 2]

def raw_cn(lr,tf,dl,sex_chr):
    """
    Considering logratio (lr), tumor fraction (tf), diploid level (dl), return the absolute CN.
    """
    if sex_chr == 'chrX':
        return ((2**(lr-dl)-1+tf) / tf)*2 # just double the number of copies for the X chr
    else:
        return (2**(lr-dl+1)-2*(1-tf)) / tf 

bin_size = 6 * 10**6  # 6 megabases

df_list = []
for sample in pts["Sample"]:
    print(sample)
    # some variables about the sample that will be useful
    pt = sample.split("-cfDNA")[0]
    WGD_status = tf[tf["Sample"] == sample]["WGD"]
    tf_value = tf[tf["Sample"] == sample]["Final_ct_fraction"]
    dl = tf[tf["Sample"] == sample]["Diploid level"]
    # 
    igv_path = os.path.join(DIR_igv_tracks, sample + "_logratio.igv")
    df = pd.read_csv(igv_path, delimiter = "\t")
    df.columns = ["CHROM", "START", "END", "FEATURE", "LOGRATIO"]
    df = df.groupby("CHROM").apply(lambda group: bin_chrom(group, bin_size)) # bin the chromosomes 
    df = df.rename(columns={"CHROM": "Chrom"})
    df["Absolute_CN"] = df.apply(lambda row: raw_cn(lr = row["LOGRATIO"], tf = tf_value, dl = dl, sex_chr = row["Chrom"]), axis=1)
    median_values_df = df.groupby(['CHROM', 'Bin'])['Absolute_CN'].median().to_frame() # mean value within each bin within each chromosome
    if all(WGD_status) == True:  
        median_values_df["Amplification"] = median_values_df["Absolute_CN"] > 5
        median_values_df["Loss"] = median_values_df["Absolute_CN"] < 3
    else: 
        median_values_df["Amplification"] = median_values_df["Absolute_CN"] > 3
        median_values_df["Loss"] = median_values_df["Absolute_CN"] < 1
    df_list.append(median_values_df)

# combine and get percentages for each bin
combined_df = pd.concat(df_list)
bin_percentages = combined_df.groupby(['CHROM', 'Bin'])[['Amplification', 'Loss']].mean() * 100
bin_percentages["Loss"] = bin_percentages["Loss"]*-1
bin_percentages.reset_index(inplace=True)
bin_percentages["Position"] = bin_percentages["CHROM"].astype(str) + " " + bin_percentages["Bin"].astype(str)
bin_percentages["Position_enumerated"] = bin_percentages['Position'].map({sample: i / (len(bin_percentages['Position'].unique()) - 1) for i, sample in enumerate(bin_percentages['Position'].unique())})

#==============================================================================
# Prepare the genes to annotate on the plot
genes = pd.DataFrame({"GENE": ["BRCA2", "PTEN", "RB1", "TP53", "MYC", "AR", "NKX3-1"], "CHROM": ["chr13", "chr10", "chr13", "chr17", "chr8", "chrX", "chr8"], "POSITION": [32315086, 87862563, 48303744, 7661779, 127735434, 67544021, 23678693]})
# genes['Bin'] = pd.cut(genes['POSITION'], bins=bin_percentages['Bin'], labels=bin_percentages['Bin'], right=False)

# Convert the gene positions to bins using pandas merge function
genes = pd.merge(genes, bin_percentages, on=['CHROM'], how='left')
genes["Position in Bin"] = genes.apply(lambda row: row["POSITION"] in row["Bin"], axis=1) # find the bin the gene position belongs to 
genes = genes[genes["Position in Bin"]]
genes["Subplot"] = genes["CHROM"].str.replace("chr", "") # this will help figure out which subplot the gene belongs to
genes.loc[genes["CHROM"] == "chrX", "Subplot"] = 23
genes.loc[genes["CHROM"] == "chrY", "Subplot"] = 24
genes["Subplot"] = genes["Subplot"].astype(int) - 1 
#==============================================================================



#==============================================================================
# Plotting
chrom_sizes = pd.read_csv(PATH_chrom_sizes, delimiter = "\t")
proportions = [size / sum(chrom_sizes["SIZE"]) for size in chrom_sizes["SIZE"]]

# Create the subplots using GridSpec with specified width ratios
fig = plt.figure(figsize=(15, 5))
gs = GridSpec(1, len(chrom_sizes), width_ratios=proportions, wspace = 0)

# to annotate on the plot

ax_dict = {}
for chrom in range(23): # dont include the Y chromosome
    ax = fig.add_subplot(gs[0,chrom])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(-100, 100)
    # annotate some genes on the plot
    ax_dict[chrom] = ax # add the ax to the dict   
    # only the very left plot has ticks
    if chrom == 0:
        ax.set_yticks([-100, -50, 0, 50, 100])
        ax.set_yticklabels(["100", "50", "0", "50", "100"])
        trans = ax.get_xaxis_transform() # x in data untis, y in axes fraction
        ax.set_ylabel("Patients with                        Patiens with\n        loss(%)                     amplification (%)")
    else:
        ax.set_yticks([])  # Remove y-axis ticks for other subplots
        ax.spines["left"].set_visible(False)
    ax.set_xticks([])
    if chrom in range(22): # if not sex chromosome
        chrom_percentages = bin_percentages[bin_percentages["CHROM"] == "chr" + str(chrom + 1)] # the if loop is 0 based but our chr start with 1
        ax.set_xlabel(str(chrom + 1))
    elif chrom == 22: 
        chrom_percentages = bin_percentages[bin_percentages["CHROM"] == "chrX"]
        ax.set_xlabel("X")
    # amplification
    # ax.bar(x = chrom_percentages["Position_enumerated"], height = chrom_percentages["Amplification"], color = "red", zorder = 10, width = 0.001)
    lowess_amplification = sm.nonparametric.lowess(chrom_percentages["Amplification"], chrom_percentages["Position_enumerated"], frac=0.1)
    ax.plot(lowess_amplification[:, 0], lowess_amplification[:, 1], color="red", linewidth=1, label="LOESS (Amplification)") # plot the loess curve
    ax.fill_between(chrom_percentages["Position_enumerated"], lowess_amplification[:, 1], 0, color="red", alpha=0.3) # fill between the curve and the axis
    # loss
    # ax.bar(x = chrom_percentages["Position_enumerated"], height = chrom_percentages["Loss"], color = "blue", zorder = 10, width = 0.001)
    lowess_loss = sm.nonparametric.lowess(chrom_percentages["Loss"], chrom_percentages["Position_enumerated"], frac=0.1)
    ax.plot(lowess_loss[:, 0], lowess_loss[:, 1], color="blue", linewidth=1, label="LOESS (Loss)")
    ax.fill_between(chrom_percentages["Position_enumerated"], lowess_loss[:, 1], 0, color="blue", alpha=0.3)
    # annotate the genes on the plot
    genes_to_annotate = genes[genes["Subplot"] == chrom]
    if not genes_to_annotate.empty: # if there are genes to annotate in a given subplot
        for i, row in genes_to_annotate.iterrows():
            x_position = row["Position_enumerated"]
            # Get the y-values from the LOESS curve
            if abs(row["Amplification"]) > abs(row["Loss"]):
                # Use the LOESS curve for amplification
                y_values = lowess_amplification[:, 1]
                label = "Amplification"
            else:
                # Use the LOESS curve for loss
                y_values = lowess_loss[:, 1]
                label = "Loss"
            # Find the nearest y-value to the x_position
            idx = (np.abs(lowess_amplification[:, 0] - x_position)).argmin()
            # Annotate the gene in the correct subplot
            y_values[idx]
            if chrom != 22: # X chromosome
                ax.annotate(row["GENE"], xy=(x_position, y_values[idx]), xytext=(x_position - 0.02, y_values[idx] + 10))
            else: # for chrX, the x values of the text get shifted a bit less.
                ax.annotate(row["GENE"], xy=(x_position, y_values[idx]), xytext=(x_position - 0.01, y_values[idx] + 10))

legend_dict = {"Amplification": "Red", 
               "Loss": "Blue"}       
  
handles = []
for key in legend_dict:
    handle = mpatches.Patch(color = legend_dict.get(key), label = key, alpha = 0.3, edgecolor = legend_dict.get(key))
    handles.append(handle)

legend = fig.legend(handles=handles, bbox_to_anchor=(0.2, 0.98), frameon=False, title = "Copy number changes", title_fontsize = 12)
legend._legend_box.align = "left"


plt.tight_layout()
fig.savefig("/home/vpc/datasets/plasma_chipseq/pipeline/figures/genome_wide_CN.png")

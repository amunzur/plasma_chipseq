# Script to make IGV snapshots with tumor and WBC bams shown together. 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 13:45:38 2021
@author: amunzur
"""
import os
import pandas as pd
import shutil
from shutil import copyfile

def make_igv_batch_script(
	variants_df, 
	igv_script, 
	DIR_snapshots, 
	add_prefix, 
	add_suffix, 
	given_range):

	for index, row in variants_df.iterrows(): # iterate through each snv in the maf

		BAM_tumor = str(add_prefix + str(row["Path_bam_t"]) + add_suffix)
		BAM_wbc = str(add_prefix + str(row["Path_bam_n"])+ add_suffix)
				
		position = int(row["Position"])
		if given_range: # if the user wants to see a given range, consider the end position of the range as well
			start_position = str(int(position - given_range/2))
			end_position = str(int(position + given_range/2))
			
		# os.remove(IGV_script)
		os.makedirs(DIR_snapshots, exist_ok = True) # make the snapshost dir if it doesnt exist already
		# output_file_name = str(row["Protein_annotation"]) + "_" + row["Gene"] + "_" + row["Sample_name_t"] + ".png" # one snapshot for each variant
		output_file_name = row["Gene"] + "_" + str(row["Position"]) + "_" + row["Ref"] + "_" + row["Alt"] + "_" + row["Patient_ID"] + ".png" # one snapshot for each variant
		
		# Begin compiling the batch file
		with open(igv_script, 'a') as the_file:
			the_file.write('new\n')
			if index == 0: the_file.write('genome hg38\n') # only load the genome if it is the first time we are starting IGV
			the_file.write(str("load " + BAM_tumor + '\n'))
			the_file.write(str("load " + BAM_wbc + '\n'))
			the_file.write(str("snapshotDirectory " + DIR_snapshots + '\n'))
				
			chrom = str(row["Chrom"])
			if given_range: the_file.write(str('goto ' + chrom + ":" + start_position + "-" + end_position + "\n"))
			else: the_file.write(str('goto ' + chrom + ":" + str(position) + "\n"))
			the_file.write('sort start location\n')
			the_file.write('sort base\n')
			the_file.write(str('snapshot ' + output_file_name + '\n'))
			the_file.write("\n")

	with open(igv_script, 'a') as the_file: the_file.write('exit') # append an exist statement at the end so that IGV closes on its own
		
	print("/home/amunzur/IGV_Linux_2.11.3/igv.sh --batch ", igv_script) # print the exact command needed to turn IGV to terminal


# PATH_variants_df = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined/tumor_wbc.csv"
PATH_variants_df = "/groups/wyattgrp/users/jbacon/projects/CHIP/variant_calling/chip/new_intersected_list.tsv"
igv_script = "/groups/wyattgrp/users/jbacon/projects/CHIP/variant_calling/chip/new_intersected_list_IGV_BATCHFILE.txt" # batch script
DIR_snapshots = "/groups/wyattgrp/users/jbacon/projects/CHIP/variant_calling/chip/snapshots"
DIR_bams = "/groups/wyattgrp/users/jbacon/projects/CHIP/alignments"
given_range = 100
add_prefix = ""
add_suffix = ""

# put the variants df in correct format
df = pd.read_csv(PATH_variants_df, "\t")[['Sample_ID', 'CHROM', 'POS', 'Gene', 'REF', 'ALT']]
df.columns = ["Patient_ID", "Chrom", "Position", "Gene", "Ref", "Alt"]
df['Patient_ID'] = df['Patient_ID'].str.replace('.table', '')

# add the path to the cfdna and wbc bams
cfDNA_bams=[]
wbc_bams=[]
for filename in os.listdir(DIR_bams):
    if filename.endswith('.bam') and 'cfDNA' in filename:
        cfDNA_bams.append(os.path.join(DIR_bams, filename))
    elif filename.endswith('.bam') and 'gDNA' in filename:
	    wbc_bams.append(os.path.join(DIR_bams, filename))
#    
cfDNA_df=pd.DataFrame(cfDNA_bams, columns=['Path_bam_t'])
gDNA_df=pd.DataFrame(wbc_bams, columns=['Path_bam_n'])

cfDNA_df['Patient_ID'] = cfDNA_df['Path_bam_t'].apply(lambda x: os.path.basename(x).split('_cfDNA')[0])
gDNA_df['Patient_ID'] = gDNA_df['Path_bam_n'].apply(lambda x: os.path.basename(x).split('_gDNA')[0])

variants_df = pd.merge(df, cfDNA_df)
variants_df = pd.merge(variants_df, gDNA_df)

# Make sure an older version of the script doesnt exist 
try:
	os.remove(PATH_batch)
except:
	print("Error while deleting batch files")

make_igv_batch_script(variants_df, PATH_batch, DIR_snapshots, add_prefix, add_suffix, given_range)

# to run this script execute on the terminal
# /home/amunzur/IGV_Linux_2.11.3/igv.sh --batch /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/IGV_batch_scripts/tumor_wbc/tumor_wbc.txt
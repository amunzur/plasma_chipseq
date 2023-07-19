"""
This script will go through the combined SVs and generate individual FASTAs for each SV. 
"""
import pandas as pd
import os
import re
import numpy as np
import subprocess

PATH_sv = "/home/vpc/datasets/plasma_chipseq/pipeline/analysis/Structural_variants_panel.csv"
PATH_blat_script = "/home/vpc/datasets/plasma_chipseq/pipeline/rearrangements/blat_script.bash"
PATH_hg38 = "/home/vpc/homo_sapiens/hg38.fa"
DIR_fasta = "/home/vpc/datasets/plasma_chipseq/pipeline/rearrangements/sv_fasta"
DIR_blat_result = "/home/vpc/datasets/plasma_chipseq/pipeline/rearrangements/blat_results"

def QC_blat_output(path_blat):
    """
    Given the path to a BLAT result related to an SV, go over the BLAT result.
    """
    df = pd.read_csv(path_blat, skiprows = 5, sep = "\t", header = None)
    df.columns = ["match", "mismatch", "rep.match", "N's"," Q gap count", "Q gap bases", "T gap count", "T gap bases", 
                  "strand", "Q name", "Q size", "Q start", "Q end", "T name", "T size", "T start", "T end", "block count", "blockSizes", "qStarts", "tStarts"]

sv = pd.read_csv(PATH_sv)

# Make FASTAs for each SV containing the supporting reads 
for idx,row in sv.iterrows():
    sample = row["SAMPLE"]
    reads = row["SUPPORTING READS"].split(";") # a list where each element is a read
    with open(os.path.join(DIR_fasta, sample+"_SV"+str(idx)+".fa"), "w") as file:
        for idx,read in enumerate(reads):
            file.write(">" + str(idx) + "\n")
            file.write(read + "\n")
            
# Run BLAT - this part saves each blat command to a script
for fasta in os.listdir(DIR_fasta): 
    blat_result = os.path.join(DIR_blat_result, fasta.replace(".fa", ".txt"))
    bash_code = "/home/vpc/tools/blat " + PATH_hg38 + " " + os.path.join(DIR_fasta, fasta) + " " + blat_result
    with open(PATH_blat_script, "a") as file: # the bash script that will have the blat commands
        file.write(bash_code + " &" + "\n")

# this part is for QCing the blat outputs
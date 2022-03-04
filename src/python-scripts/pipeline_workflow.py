# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: 'Python 3.8.10 64-bit (''satay-dev'': conda)'
#     name: python3
# ---

# ## Notebook showing a general workflow for the postprocessing analysis for the satay pipeline 

# +

import os, sys
import warnings
import timeit
import numpy as np
import pandas as pd 


from transposonmapper.processing.clean_bedwigfiles import cleanfiles
from transposonmapper.processing.transposonread_profileplot_genome import profile_genome
from transposonmapper.processing.genomicfeatures_dataframe import dna_features
from transposonmapper.statistics import volcano




# +
wig_files=[]
bed_files=[]
pergene_files=[]
#data_dir= "../satay/data_files/data_unmerged/"
#data_dir="../transposonmapper/data_files/files4test/"
data_dir="../data/"
#data_dir="../transposonmapper/data_files/"
for root, dirs, files in os.walk(data_dir):
    for file in files:
        if file.endswith("sorted.bam.wig"):
            wig_files.append(os.path.join(root, file))
        elif file.endswith("sorted.bam.bed"):
             bed_files.append(os.path.join(root, file))
        elif file.endswith('sorted.bam_pergene_insertions.txt'):
            pergene_files.append(os.path.join(root, file))





# -

# ### Clean the wig and bed files generated by transposon mapper (it is ok if run in spyder)

## clean wig files for proper visualization in the genome Browser http://genome-euro.ucsc.edu/cgi-bin/hgGateway
custom_header = ""
split_chromosomes = False
for files in zip(wig_files,bed_files):
    cleanfiles(filepath=files[0], custom_header=custom_header, split_chromosomes=split_chromosomes)
    cleanfiles(filepath=files[1], custom_header=custom_header, split_chromosomes=split_chromosomes)


# +
data_dir="../data/"
cleanbed_files=[]
for root, dirs, files in os.walk(data_dir):
    for file in files:
        if file.endswith("clean.bed"):
            cleanbed_files.append(os.path.join(root, file))

cleanwig_files=[]
for root, dirs, files in os.walk(data_dir):
    for file in files:
        if file.endswith("clean.wig"):
            cleanwig_files.append(os.path.join(root, file))
# -

cleanbed_files




# +
#transposonread_profileplot_genome.py (to check the insertion and read distribution throughout the genome)
#example of file to analyse with profile_genome
# bed_file=cleanbed_files[0]
variable="reads" #"reads" "transposons"
bar_width=None
savefig=True

for bed_file in cleanbed_files:

    profile=profile_genome(bed_file=bed_file, variable=variable, bar_width=bar_width, savefig=savefig,showfig=True)


# -

cleanbed_files=["../data/WT_1-Benoit/ERR1533147_trimmed.sorted.bam_clean.bed",
"../data/WT_2-Benoit/ERR1533148_trimmed.sorted.bam_clean.bed"]
pergene_files=["../data/WT_1-Benoit/ERR1533147_trimmed.sorted.bam_pergene_insertions.txt",
"../data/WT_2-Benoit/ERR1533148_trimmed.sorted.bam_pergene_insertions.txt"]
cleanwig_files=["../data/WT_1-Benoit/ERR1533147_trimmed.sorted.bam_clean.wig",
"../data/WT_2-Benoit/ERR1533148_trimmed.sorted.bam_clean.wig"]


# +
# genomic features 
# i=0

# to make it only for certain files 
#cleanwig_files=[cleanwig_files[6],cleanwig_files[7],cleanwig_files[13]]
for i in range(len(cleanwig_files)):
    wig_file = cleanwig_files[i]
    pergene_insertions_file = pergene_files[i]
    plotting=False
    variable="reads" #"reads" or "insertions"
    savefigure=False
    verbose=True

    dna_df2=[]
    for chrom in ['I', 'II', 'III', 'IV','V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']:
    #     region=chrom
    
        region = chrom #e.g. 1, "I", ["I", 0, 10000"], gene name (e.g. "CDC42")
        dna_df2.append(dna_features(region=region,
                    wig_file=wig_file,
                    pergene_insertions_file=pergene_insertions_file,
                    variable=variable,
                    plotting=plotting,
                    savefigure=savefigure,
                    verbose=verbose))
    data_genome=pd.concat(dna_df2, axis=0,ignore_index=True)
    data_genome.to_excel('../postprocessed-data/'+ cleanwig_files[i].split("/")[2] + ".xlsx",engine='openpyxl')
# -



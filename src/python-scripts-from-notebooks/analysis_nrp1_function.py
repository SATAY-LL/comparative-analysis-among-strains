# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: 'Python 3.9.7 64-bit (''transposonmapper'': conda)'
#     language: python
#     name: python3
# ---

# +
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys
from collections import defaultdict
from ast import literal_eval

from from_excel_to_list import from_excel_to_list
from transposonmapper.statistics import volcano

# +
## Plot transposons along the whole genome


## Importing pergene files 

pergene_files=[]
#data_dir= "../satay/data_files/data_unmerged/"
#data_dir="../transposonmapper/data_files/files4test/"
data_dir="../postprocessed-data/"
#data_dir="../transposonmapper/data_files/"
for root, dirs, files in os.walk(data_dir):
    for file in files:
        if file.endswith("pergene_insertions.xlsx"):
            pergene_files.append(os.path.join(root, file))

list_data=[]
for i in pergene_files:
    list_data.append(pd.read_excel(i,engine='openpyxl',index_col="Unnamed: 0"))

# -

keys= ['wt_merged','dbem3_merged']

# ## Volcano plots

# +
path_a = r"../data/"
filelist_a = ["wt_a/WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_pergene.txt",
"wt_b/WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_pergene.txt"]
path_b = r"../data/"
filelist_b = ["dnrp1_a/dnrp1-1_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt",
"dnrp1_b/dnrp1-2_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt"]


variable = 'tn_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.01 #set threshold above which p-values are regarded significant
normalize=True

trackgene_list = ['nrp1','bem3','bem1','bem2'] # ["cdc42"]


figure_title = "WT vs dnrp1"

volcano_df_nrp1_wt = volcano(path_a=path_a, filelist_a=filelist_a,
            path_b=path_b, filelist_b=filelist_b,
            variable=variable,
            significance_threshold=significance_threshold,
            normalize=normalize,
            trackgene_list=trackgene_list,
            figure_title=figure_title)

# +
volcano_df_nrp1_wt.sort_values(by=['fold_change'], inplace=True)

volcano_df_nrp1_wt[volcano_df_nrp1_wt["significance"]==True][0:20]

# +
## Plot number of normalized insertions from WT/ total number of insertions of the library 
# vs the same for dnrp1

# import excel file with the normalized data

data_norm_pd=pd.read_excel("./postprocessed-data/data_norm_linear_transformation_all_backgrounds.xlsx",
engine='openpyxl',index_col="Unnamed: 0")

# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3.8.10 ('satay-dev')
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
from scipy.stats import norm

from from_excel_to_list import from_excel_to_list
from functions_analysis_frompergene2fitness import reads_per_insertion_along_gene_length
from functions_analysis_frompergene2fitness import genes_discarded4fitness_from_flanking_regions
from functions_analysis_frompergene2fitness import protein_domains_info,reads_and_insertions_per_domain
from functions_analysis_frompergene2fitness import fitness_models

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

keys=[]
for i in np.arange(0,len(pergene_files)):
    keys.append(pergene_files[i].split("/")[-1].split("_")[0]+"_"+pergene_files[i].split("/")[-1].split("_")[1])

list_data_pd=pd.concat(list_data,axis=0,keys=keys)

# -

keys,len(keys)

wig_files=["../data/dbem3_b/yLIC137_8_merged_cleaned_forward_reads_trimmed.sorted.bam_clean.wig",
"../data/WT_1-Benoit/ERR1533147_trimmed.sorted.bam_clean.wig",
"../data/dnrp1_2/dnrp1-2_merged-DpnII-NlaIII-a_trimmed.sorted.bam_clean.wig",
"../data/bem1-aid_a/yWT03a_16_trimmed_out_restriction_sites_yWT03a_16_merged_cleaned_forward_reads_trimmed.sorted.bam_clean.wig",
"../data/dbem1dbem3_b/yTW001_6_merged_cleaned_forward_reads_trimmed.sorted.bam_clean.wig",
"../data/wt_merged/WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_clean.wig",
"../data/dbem3_a_trimmed_restriction_sites/ylic137_7_trimmed_out_restriction_sites_all_cleaned_fw_reads_trimmed.sorted.bam_clean.wig",
"../data/dbem1dbem3_a/yTW001_4_merged_cleaned_forward_reads_trimmed.sorted.bam_clean.wig",
"../data/bem1-aid-dbem3_a/yWT04a_14_merged_cleaned_forward_reads_trimmed.sorted.bam_clean.wig",
"../data/bem1-aid_merged/merged_datasets_trimmed_out_restriction_sites_merged_yTW003_new_trimmed.sorted.bam_clean.wig",
"../data/bem1-aid-dbem3_b/yWT04a_23_merged_cleaned_forward_reads_trimmed.sorted.bam_clean.wig",
"../data/dnrp1_1/dnrp1-1_merged-DpnII-NlaIII-a_trimmed.sorted.bam_clean.wig",
"../data/wt_b/WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_clean.wig",
"../data/wt_a/WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_clean.wig",
"../data/dnrp1_merged/dnrp1_merged_dnrp1-1_dnrp1-2_trimmed.sorted.bam_clean.wig",
"../data/WT_2-Benoit/ERR1533148_trimmed.sorted.bam_clean.wig",
"../data/bem1-aid_b/yWT0321_a_trimmed_out_restriction_sites_yWT0321_a_merged_cleaned_forward_reads_trimmed.sorted.bam_clean.wig",
"../data/dbem3_merged/merged_ylic137_trimmed.sorted.bam_clean.wig",
"../data/dbem3_a/all_cleaned_fw_reads_trimmed.sorted.bam_clean.wig"]

i=0
fitness_all=[]
for background in keys:
    
    r,gene_coordinates,reads_location,insertion_locations=reads_per_insertion_along_gene_length(list_data_pd,background,number_of_parts=10)
    flanking_regions,discarded_genes=genes_discarded4fitness_from_flanking_regions(list_data_pd,background,wig_files[i],windows_size=20000)
    data_domains=protein_domains_info(list_data_pd,background,gene_coordinates,reads_location,
                                  insertion_locations,discarded_genes)
    data_domains_extended=reads_and_insertions_per_domain(list_data_pd,background,
                                                        data_domains,reads_location,insertion_locations)
    fitness_models_pd=fitness_models(list_data_pd,background,data_domains_extended,r,discarded_genes)
    
    i=i+1
    fitness_all.append(fitness_models_pd)

# +
import pickle

with open("../postprocessed-data/fitness_models_all_backgrounds", "wb") as fp:   #Pickling
    pickle.dump(fitness_all, fp)
# -

# getting all discarded genes per background
i=0
discarded_genes_all=[]
for background in keys:
    

    flanking_regions,discarded_genes=genes_discarded4fitness_from_flanking_regions(list_data_pd,background,wig_files[i],windows_size=20000)
    discarded_genes_all.append(discarded_genes)
    i=i+1

discarded_genes_all

# +
import pickle

with open("../postprocessed-data/discarded_genes_all_backgrounds", "wb") as fp:   #Pickling
    pickle.dump(discarded_genes_all, fp)
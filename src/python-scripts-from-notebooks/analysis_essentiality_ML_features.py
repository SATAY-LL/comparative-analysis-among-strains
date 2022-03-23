# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: 'Python 3.8.10 64-bit (''satay-dev'': conda)'
#     language: python
#     name: python3
# ---

# ## This notebook will implement the major influencing features for essentiality prediction using ML , which are :
#
# - Neighboorhood index : Number of transposon insertions within the ORF, normalized by the length of the ORF and the surrounding 10kbp. 
# - Freedom index :  Length of the largest insertion-free region in the ORF (Open Reading Frame), divided by the ORF’s length.
#
# From paper : Levitan A, et al. (2020) Comparing the utility of in vivo transposon mutagenesis approaches in yeast species to infer gene essentiality.

# +
## Importing the required python libraries 
import os, sys
import warnings
import timeit
import numpy as np
import pandas as pd 
import pkg_resources
import matplotlib.pyplot as plt
import re
import seaborn as sns
from collections import defaultdict

from from_excel_to_list import from_excel_to_list

# +
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

# +
# import essential genes used in transposonmapper

essentials_satay=pd.read_csv("../postprocessed-data/Cerevisiae_AllEssentialGenes_List.txt",header=0
,sep="\t")

essentials_satay.columns=["gene name"]

# import conversion file from systematic names to standard names 
conversion=pd.read_csv("../postprocessed-data/from_systematic_genenames2standard_genenames.csv",
header=0,sep=",")

conversion.columns=["systematic name","standard  name"]

# save the standard names of the essential genes in a systematic format
standard_essentials=[]
for names in essentials_satay.loc[:,"gene name"]:
    
    if names in conversion["systematic name"].values:
        standard_essentials.append(conversion.loc[conversion["systematic name"]==names]["standard  name"].values[0])


# -

polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")
polarity_genes.fillna(0,inplace=True)

data_norm_pd=pd.read_excel("../postprocessed-data/data_norm_linear_transformation_per_background.xlsx",
engine='openpyxl',index_col="background")
data_norm_pd.drop(columns=["Unnamed: 0","Unnamed: 1"],inplace=True)


# +
list_data_pd_wt=list_data_pd.loc["wt_merged"]

ni_essentials=[] # neighborhood index
ni=data_norm_pd.loc["wt_merged"]["tr_normalized_windows"].values

ho=np.where(list_data_pd_wt.loc[:,"Gene name"]=="HO")[0][0]

ni_ho=data_norm_pd.loc["wt_merged"]["tr_normalized_windows"].values[ho]
ni_norm=ni/ni_ho


genes=list_data_pd_wt.loc[:,"Gene name"]
genes.reset_index(drop=True,inplace=True)

for i in standard_essentials:
    x=np.where(genes==i)[0]
    
    ni_essentials.append(ni_norm[x])


ni_essentials=np.array(ni_essentials)

ni_essentials_bellow_ho=ni_essentials[(ni_essentials<np.array(0.5))]

d=len(ni_essentials_bellow_ho)/len(standard_essentials)

# +
fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

axes.hist(np.concatenate(ni_essentials),bins=200,alpha=0.5,color="gray",label="essential genes");
axes.vlines(0.5,0,200,color="red",linestyle="dashed",linewidth=2,label="HO/2");
axes.annotate(f"{d*100:.2f}" +"%",xy=(0.1,175),fontsize=16)
axes.set_title("Distribution of normalized insertions for annotated essential genes",fontsize=16)
axes.tick_params(axis="both",labelsize=16)
axes.set_xlabel("Normalized insertions by a 10kB windows compared to HO locus",fontsize=16)
axes.set_xlim(0,3)
axes.set_ylabel("Count",fontsize=16)
axes.legend(loc="best")
plt.tight_layout()

fig.savefig("../figures/figures_thesis_chapter_2/fig_distribution_of_normalized_insertions_for_annotated_essential_genes.png",dpi=400)

# -



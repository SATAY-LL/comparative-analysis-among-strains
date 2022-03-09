# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
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

keys= ['wt_merged','dnrp1_merged']

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

# import excel file with the normalized data

data_norm_pd=pd.read_excel("../postprocessed-data/data_norm_linear_transformation_per_background.xlsx",
engine='openpyxl',index_col="background")
data_norm_pd.drop(columns=["Unnamed: 0","Unnamed: 1"],inplace=True)
# -

polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")

polarity_genes.fillna(0,inplace=True)
polarity_genes

# +
## Plot number of normalized insertions from WT/ total number of insertions of the library 
# vs the same for dnrp1

keys= ['wt_merged','dnrp1_merged']
data=data_norm_pd.loc[keys]

variable=["tr_normalized_windows" ,"reads_normalized_windows"]
var_norm=["Insertions" ,"Reads"]
x_1=data.loc["wt_merged"][variable[0]]/data.loc["wt_merged"][var_norm[0]].sum()
y_1=data.loc["dnrp1_merged"][variable[0]]/data.loc["dnrp1_merged"][var_norm[0]].sum()

x_0=data.loc["wt_merged"][variable[1]]/data.loc["wt_merged"][var_norm[1]].sum()
y_0=data.loc["dnrp1_merged"][variable[1]]/data.loc["dnrp1_merged"][var_norm[1]].sum()

# x_1=data.loc["wt_merged"]["reads_over_windows"]
# y_1=data.loc["dnrp1_merged"]["reads_over_windows"]

# x_0=data.loc["wt_merged"]["insertions_over_windows"]
# y_0=data.loc["dnrp1_merged"]["insertions_over_windows"]

# x=data.loc["wt_merged"]["Reads"]
# y=data.loc["dnrp1_merged"]["Reads"]

fig , ax = plt.subplots(ncols=2,figsize=(6,2))

plt.subplots_adjust(wspace=0.5)

#sns.regplot(x,y,fit_reg=True,color="black",marker="o")

ax[0].scatter(x_0,y_0,color="black",marker="o",alpha=0.5,label="Insertions")

ax[1].scatter(x_1,y_1,color="black",marker="o",alpha=0.5,label="Reads")

for axes in ax:
    # axes.set_xscale("log")
    # axes.set_yscale("log")    

    axes.set_xlabel("WT")
    axes.set_ylabel("$\Delta$nrp1")
    axes.legend()

# -

test=["wt_merged","bem1-aid_b"]
data_test=data_norm_pd.loc[test]

# +
variable=["tr_normalized_windows" ,"reads_normalized_windows"]
var_norm=["Insertions" ,"Reads"]
data=data_test

x_1=data.loc["wt_merged"][variable[0]]/data.loc["wt_merged"][var_norm[0]].sum()
y_1=data.loc["bem1-aid_b"][variable[0]]/data.loc["bem1-aid_b"][var_norm[0]].sum()

x_0=data.loc["wt_merged"][variable[1]]/data.loc["wt_merged"][var_norm[1]].sum()
y_0=data.loc["bem1-aid_b"][variable[1]]/data.loc["bem1-aid_b"][var_norm[1]].sum()
# x_1=data_test.loc["wt_merged"]["reads_over_windows"]
# y_1=data_test.loc[test[1]]["reads_over_windows"]

# x_0=data_test.loc["wt_merged"]["insertions_over_windows"]
# y_0=data_test.loc[[test[1]]]["insertions_over_windows"]
# x=data.loc["wt_merged"]["Reads"]
# y=data.loc["dnrp1_merged"]["Reads"]

fig , ax = plt.subplots(ncols=2,figsize=(6,2))

plt.subplots_adjust(wspace=0.5)

# sns.regplot(x_0,y_0,fit_reg=True,color="black",marker="o",ax=ax[0])
# sns.regplot(x_1,y_1,fit_reg=True,color="black",marker="o",ax=ax[1])

ax[0].scatter(x_0,y_0,color="black",marker="o",alpha=0.5,label="Insertions")

ax[1].scatter(x_1,y_1,color="black",marker="o",alpha=0.5,label="Reads")


for axes in ax:
    # axes.set_xscale("log")
    # axes.set_yscale("log")    
    # axes.set_xlim(0,0.01)
    # axes.set_ylim(0,0.01)
    axes.set_xlabel("WT")
    axes.set_ylabel("$\Delta$bem1")
    axes.legend()


# +
## Analysis of essential genes in WT and dnrp1
# - import the genes that are duplicated 
# - import the genes that have low insertions in the flanking regions

# To extract potential essential genes take :

# - genes with low normalized transposons counts (given by a 10kb window) NI
# - genes with the longest region void  of transposons (FI) 
# - Highest domain likelihood scores (Benoit score) 

# Compare low NI with known essential genes. 


# +
## import essential genes used in transposonmapper

essentials_satay=pd.read_csv("../postprocessed-data/Cerevisiae_AllEssentialGenes_List.txt",header=0,sep="\t")

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

data_wt=data_norm_pd.loc["wt_merged"]
data_wt.reset_index(inplace=True)
data_wt[data_wt.loc[:,"Gene name"]==standard_essentials[0]].index[0]

# +
data_wt=data_norm_pd.loc["wt_merged"].copy()
data_wt.reset_index(inplace=True)
data_wt["true essential"]=np.zeros(len(data_wt))

for true_genes in standard_essentials:
    
    index_essential=data_wt[data_wt.loc[:,"Gene name"]==true_genes].index
    if len(index_essential)>0:
        if type(index_essential[0])==np.int64:
            #scores_wt.loc[index_essential[0],"true essential"]=1
            tmp=data_wt.index[index_essential]

            data_wt.loc[tmp,"true essential"]=1
    
    
data_wt.fillna(0,inplace=True)
# -

data_wt.groupby(["true essential"]).count()

# +


g=sns.pairplot(data_wt,hue="true essential",vars=["Reads","tr_normalized_windows",
"reads_normalized_windows","tr-density","tr_normalized_windows"],
diag_kind="kde",diag_kws={"shade":True},plot_kws = {'alpha': 0.6, 's': 80, 'edgecolor': 'k'},
size = 4)

g.axes[0,0].set_xlim(0,10000)
# g.axes[0,1].set_xlim(0,100)
g.axes[0,1].set_xlim((0,0.05))
g.axes[0,2].set_xlim((0,0.05))
g.axes[0,3].set_xlim((0,0.2))

plt.tight_layout()
# -



# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
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

from functions_normalizations import linear_transformation_per_background,linear_transformations_plots,linear_transformations_plots_per_chrom



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
## Importing postproccesed data from pipeline
from ast import literal_eval
data_post=[]

for i in keys:
    data_post.append(pd.read_excel(data_dir+i+".xlsx",index_col="Unnamed: 0",engine='openpyxl'))
    
list_data_post=pd.concat(data_post,axis=0,keys=keys)
list_data_post.drop(columns=["Feature_alias","Feature_name"],inplace=True)
list_data_post.fillna(1,inplace=True)

# -

# ## Linear normalization procedure
#
# - Normalize the transposon density data per chromosome over the density of the chromosome. 
# - Between datasets : normalize datasets such that their read-counts have the same mean(e.g. by dividing each transposon site by the total read-count).
# - Refinement of this method: 
#
# *A
# refinement of this approach that is specific to TnSeq is to scale the read counts to have the
# same mean over non-zero sites (which we call ‘Non-Zero Mean Normalization’ or
# NZMean), since different datasets can have widely varying levels of saturation, and
# distributing the same number of reads over fewer TA sites will naturally inflate the mean
# read count among them.*
#
# - Limitations: 
#
# *One significant limitation of methods that linearly transform datasets is that they are susceptible to large spikes in read-counts. Because these methods multiply read-counts by a
# constant scalar value, they cannot reduce large outliers without also affecting small read-
# counts which are more common. Even if the datasets share the same mean, for instance, any
# skew in distribution of read-counts itself would still be present.* From 1.DeJesus, M. A. & Ioerger, T. R. Normalization of transposon-mutant library sequencing datasets to improve identification of conditionally essential genes. J. Bioinform. Comput. Biol. 14, 1642004 (2016).
#
# - The distribution of read-counts in most TnSeq datasets resembles a **Geometric-like
# distribution, in that read-counts at most sites are small** (i.e. 1–50), with a (rapidly)
# decreasing probability of sites with large counts. Ideally, a normalization method would
# improve detection of conditionally essential genes between conditions by eliminating any
# skew and making the datasets more closely fit this Geometric-like distribution.

backgrounds= ['wt_merged','dnrp1_merged','dbem3_merged','bem1-aid_merged',
'dbem1dbem3_a','dbem1dbem3_b']

backgrounds=keys

chrom_length=pd.read_excel("../postprocessed-data/chromosome_length.xlsx",index_col="Chromosome")
chrom_length.drop(columns=["Unnamed: 0"],inplace=True)
chrom_length.loc["I"]

# +
data_norm=[]
reads_per_chrom_list=[]
insertions_per_chrom_list=[]
for key in backgrounds:
    data_transformed,reads_per_chrom,insertions_per_chrom=linear_transformation_per_background(list_data_pd,key,
    chrom_length,windows_size=10000)
    data_norm.append(data_transformed)
    reads_per_chrom_list.append(reads_per_chrom)
    insertions_per_chrom_list.append(insertions_per_chrom)

data_norm_pd=pd.concat(data_norm,axis=0,keys=backgrounds)
reads_per_chrom_pd=pd.concat(reads_per_chrom_list,axis=0,keys=backgrounds)
insertions_per_chrom_pd=pd.concat(insertions_per_chrom_list,axis=0,keys=backgrounds)

# +
## To properly export to excel, lets add a column with the name of the background

data2excel=data_norm_pd.copy()
data2excel["background"]=data2excel.index.get_level_values(0)
data2excel.to_excel("../postprocessed-data/data_norm_linear_transformation_per_background.xlsx")
# -

data_norm_pd.loc["wt_merged"][0:3]

data_norm_pd.columns


# comparison of libraries in terms of normalized transposon densityies 
mutant="bem1-aid_merged"
data_wt= data_norm_pd.loc["wt_merged","Reads"]/data_norm_pd.loc["wt_merged","Reads"].sum()
data_mutant=data_norm_pd.loc[mutant,'Reads']/data_norm_pd.loc[mutant,'Reads'].sum()
plt.scatter(data_wt,data_mutant,color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("wt_merged")
plt.title("Reads comparison across libraries")
plt.ylabel(mutant)
# plt.ylim(0.1,100)
# plt.xlim(0.1,100)

# +
## make a function that looks for the gene names that have more than 10 as transposon density
# in the mutant background and less than one in WT, and viceversa. 
# -

data_norm_pd.loc["wt_merged"].columns

linear_transformations_plots(data_norm_pd,type="plot-genome-insertions",background="bem1-aid_merged")
linear_transformations_plots(data_norm_pd,type="plot-genome-insertions",background="wt_merged")

types=["insertions","reads","transposon density","plot-genome-reads",
"plot-genome-insertions","plot-genome-density"]
for i in types:
    linear_transformations_plots(data_norm_pd,i,"wt_merged",saveFigure=False)

# +
#linear_transformations_plots_per_chrom(data_norm_pd,"wt_merged","plot-genome-insertions",chrom="XV",saveFigure=True)

# +
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(10,5))

plt.subplots_adjust(wspace=0.8,hspace=0.6)


data=data_norm_pd.loc["wt_merged"]
data=data[~data.isin([np.nan, np.inf, -np.inf]).any(1)]    

values=["Linear-norm-reads","Reads","Linear-norm-insertions","Insertions",
"Linear-norm-tr-density","tr-density","reads_normalized_windows","Reads",
"tr_normalized_windows","Insertions","Insertions","tr-density"]

labels=["Norm. over chrom.","Raw reads","Norm. over chrom.","Raw insertions",
"Norm. over chrom.","Insertion density","Norm. over 10kb","Raw reads",
"Norm. over 10kb","Raw insertions","Raw insertions","Insertion density"]

fig.suptitle("Relationships between raw and transformed data",fontsize=18)

j=0
for i in np.arange(0,6):
    #plt.subplot(2,3,i+1)
    x=data.loc[:,values[j]]
    y=data.loc[:,values[j+1]]
    sns.regplot(x=x,y=y,ax=ax[i//3,i%3],color="black",seed=1,
    ci=95,order=1,scatter_kws={'alpha':0.3})
    # plt.scatter(x,y,color="black",alpha=0.3)
    
    ax[i//3,i%3].set_xlabel(labels[j],fontsize=16)
    ax[i//3,i%3].set_ylabel(labels[j+1],fontsize=16)
    # # #plt.xscale("log")
    # # #plt.yscale("log")
    # # plt.ylim(0,np.max(y)/2)
    # # plt.xlim(0,np.max(x)/2)
    # plt.yticks(fontsize=14)
    # plt.xticks(fontsize=14)
    j=j+2

for axes in ax.flatten():
    axes.tick_params(axis='x', labelsize=14)
    axes.tick_params(axis='y', labelsize=14)
    

fig.savefig("../figures/figures_thesis_chapter_2/linear_relationships-data.png",bbox_inches="tight",dpi=400)
#ax[1,2].set_axis_off()

# -

# ## There is a problem with the pergene insertion file!!!!!
#
# ### Problem:
#  **The data does not go continously from chromosome 1 to 16.** After chromosome 4 it jumps to chromosome 9,  then goes to Mitochondrial chromosome and then continues to chromosome 5 until 16. 
#
# ### Consequence:
#
# **We can not plot the pergene insertion file as it is because the data is not continously from chromosome 1 to 16.**
#
# ## Possible solutions:
#
# - Plot the data, from the pergene insertions file by matching the data to every chromosome number , not as it is. 
#
#

# +
## Compare different libraries 



# +
## Make the following plots:

# 1. Plot the reads per chrom 
# 2. Plot the insertions per chrom 
# 3. Plot the reads over the 10kb windows 
# 4. Plot the insertions over the 10kb windows 

data=data_norm_pd.loc["wt_merged"]

fig,ax=plt.subplots(nrows=2,ncols=2,figsize=(10,5))
plt.subplots_adjust(hspace=0.6,wspace=0.3)

chrom=["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]
x_reads=[]
x_insertions=[]
for chr in chrom:
    x_reads.append(reads_per_chrom_pd.loc["wt_merged"][chr])
    x_insertions.append(insertions_per_chrom_pd.loc["wt_merged"][chr])


ax[0,0].bar(np.arange(0,len(chrom),1),x_reads,color="black",alpha=0.5)
# sns.barplot(x=chrom,y=reads_per_chrom_pd.loc["wt_merged"][chrom],
# data=reads_per_chrom_pd.loc["wt_merged"],ax=ax[0,0],color="black",ci=0.95);
ax[0,0].set_xticks(np.arange(0,len(chrom),1))
ax[0,0].set_xticklabels(chrom,rotation=45,fontsize=12);
ax[0,0].set_title("Reads per chrom",fontsize=16)
ax[0,0].set_xlabel(" ",fontsize=16)


ax[0,1].bar(np.arange(0,len(chrom),1),x_insertions,color="black",alpha=0.5)
# sns.barplot(x=chrom,y=insertions_per_chrom_pd.loc["wt_merged"][chrom],
# data=insertions_per_chrom_pd.loc["wt_merged"],ax=ax[0,1],color="black",ci=0.95);
ax[0,1].set_xticks(np.arange(0,len(chrom),1))
ax[0,1].set_xticklabels(chrom,rotation=45,fontsize=14);
ax[0,1].set_title("Insertions per chrom",fontsize=16)

ax[0,1].set_xlabel(" ",fontsize=16)


ax[1,0].bar(np.arange(0,len(data)),data["reads_over_windows"],color="black",alpha=0.5)
# sns.lineplot(data=data["reads_over_windows"],ax=ax[1,0],color="black",ci="sd");

ax[1,0].set_title("Reads over 10kb windows",fontsize=16)
ax[1,0].set_ylabel(" ")
ax[1,0].set_xlabel("Genomic locations")


# sns.lineplot(data=data["insertions_over_windows"],ax=ax[1,1],color="black");
ax[1,1].bar(np.arange(0,len(data)),data["insertions_over_windows"],color="black",alpha=0.5)
ax[1,1].set_title("Insertions over 10kb windows",fontsize=16)
ax[1,1].set_ylabel(" ")
ax[1,1].set_xlabel("Genomic locations",fontsize=14)


for i in range(ax.shape[0]):
    for j in range(ax.shape[1]):
        ax[i,j].tick_params(axis='x', labelsize=14)
        ax[i,j].tick_params(axis='y', labelsize=14)

fig.savefig("../figures/figures_thesis_chapter_2/reads_insertions_per_chrom.png",bbox_inches="tight",dpi=400)
# -

reads_per_chrom_pd.loc["wt_merged"]

# ## Make a quantile quantile plot of the reads counts over the genes to have an idea of how deviated from a normal distribution is this. 
#
# - Sources: https://towardsdatascience.com/understand-q-q-plot-using-simple-python-4f83d5b89f8f
# - https://seaborn-qqplot.readthedocs.io/en/latest/

# ## Non linear normalization procedure
#
# - Normalize the reads data per chromosome to fit a beta geometric distribution with parameters alpha and beta.
# - The beta geometric distribution (also called the Type I Geometric) is a type of geometric distribution, where the probability of success parameter, p, has a Beta distribution with shape parameters alpha(α) and beta(β); both shape parameters are positive (α > 0 and β > 0).
# - It is a type of compound distribution.
#
# - See an example here: https://stackoverflow.com/questions/59308441/fitting-for-discrete-data-negative-binomial-poisson-geometric-distribution
#
# - https://www.geeksforgeeks.org/python-discrete-geometric-distribution-in-statistics/

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
# -

polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")
polarity_genes.fillna(0,inplace=True)


def rates_non_coarsed_model(data,background):

    data_background=data.loc[background]
    data_total_reads=np.sum(data.loc[: , "Reads"])
    Time=90 # reseeding time

    data_background_reads=[]

    rates=defaultdict(dict)
    for i in np.arange(0,len(data_background)):
        data_background_reads=from_excel_to_list(data_background.loc[i,"Reads per insertion location"])

        data_background_reads=np.array(data_background_reads)
        if data_background_reads!=[]:
            # take the data between the 10% and the 90% of the gene
            start=int(len(data_background_reads)*0.1)+1
            end=int(len(data_background_reads)*0.9)-1
            
            data_background_reads=data_background_reads[start:end]
            rates[background][i]=1/Time*np.log(data_background_reads/(1-data_background_reads/np.sum(data_background_reads)))*1/data_total_reads
        else:
            rates[background][i]=0
    rates_pd=pd.DataFrame(rates)
    return rates_pd


# +
rates=[]

for i in keys:
    rates.append(rates_non_coarsed_model(list_data_pd,background=i))

rates_pd=pd.concat(rates,axis=1)
# -

x = np.array([])
x.size

# +
mean_values_wt=[]
std_values_wt=[]

for i in np.arange(0,len(rates_pd)):
    mean_values_wt.append(np.mean(rates_pd.loc[i , "wt_merged"]))
    std_values_wt.append(np.std(rates_pd.loc[i , "wt_merged"]))
# -

list_data_pd_wt=list_data_pd.loc["wt_merged"]
gene=np.where(list_data_pd_wt.loc[:,"Gene name"]=="BEM1")[0][0]
rates_pd.loc[gene,"wt_merged"]

# +
mean_values=np.zeros((len(rates_pd),len(keys)))
std_values=np.zeros((len(rates_pd),len(keys)))

for i in np.arange(0,len(keys)):
    for j in np.arange(0,len(rates_pd)):
        mean_values[j,i]=(np.mean(rates_pd.loc[j, keys[i]]))
        std_values[j,i]=(np.std(rates_pd.loc[j , keys[i]]))



# +
index_wt_merged=np.where(np.array(keys)=="wt_merged")[0][0]

norm_mean_values=[]
norm_std_values=[]
for i in np.arange(0,len(mean_values)):
    norm_mean_values.append(mean_values[i]/mean_values[i,index_wt_merged])
    norm_std_values.append(std_values[i]/std_values[i,index_wt_merged])

norm_mean_values=np.array(norm_mean_values)
norm_std_values=np.array(norm_std_values)

norm_mean_values[~np.isfinite(norm_mean_values)] = 0
norm_std_values[~np.isfinite(norm_std_values)] = 0

norm_mean_values_pd=pd.DataFrame(norm_mean_values,columns=keys)
norm_std_values_pd=pd.DataFrame(norm_std_values,columns=keys)

# +
list_data_pd_wt=list_data_pd.loc["wt_merged"]
ho=np.where(list_data_pd_wt.loc[:,"Gene name"]=="HO")[0][0]

fitness_wt_ho=mean_values_wt[ho]
std_wt_ho=std_values_wt[ho]

values=np.array(mean_values_wt)/fitness_wt_ho
values=values[np.where(values!=-np.inf)]
values=values[np.where(values!=np.nan)]
values=values[np.where(values!=np.inf)]
#values[~np.isfinite(values)] = 0

std_values_wt=np.array(std_values_wt)/std_wt_ho
std_values_wt=std_values_wt[np.where(std_values_wt!=-np.inf)]
std_values_wt=std_values_wt[np.where(std_values_wt!=np.nan)]
std_values_wt=std_values_wt[np.where(std_values_wt!=np.inf)]

# -

fig, axes = plt.subplots(2, 1,  gridspec_kw={"height_ratios":(.10, .30)}, figsize = (8, 5))
sns.violinplot(values, ax=axes[0],color="gray",orient="h",inner="quartile")
sns.violinplot(std_values_wt, ax=axes[0],color="pink",orient="h",inner="quartile")
sns.histplot(values,bins=200,color="gray",ax=axes[1],stat="percent",label="mean-WT",kde=True,element="step")
sns.histplot(std_values_wt,bins=200,color="pink",ax=axes[1],stat="percent",label="std-WT",kde=True,element="step")
axes[1].set_xlabel("Fitness values compared to HO locus",fontsize=16)
axes[1].set_ylabel("Percent",fontsize=16)
axes[1].tick_params(axis="both",labelsize=16)
axes[0].tick_params(axis="both",labelsize=16)
axes[1].legend(fontsize=16)

# +
## Histplot for individual gene fitness gaussian distribution 
backgrounds_heatmap=["wt_a","wt_b",'wt_merged', 'dnrp1_a','dnrp1_b','dnrp1_merged', 'dbem3_a',
'dbem3_b','dbem3_merged','bem1-aid_a', 'bem1-aid_b', "bem1-aid_merged",'dbem1dbem3_a', 'dbem1dbem3_b']

gene=np.where(list_data_pd_wt.loc[:,"Gene name"]=="CDC42")[0][0]

mean=norm_mean_values_pd.loc[gene,backgrounds_heatmap]
std=norm_std_values_pd.loc[gene,backgrounds_heatmap]

fig, axes = plt.subplots(len(backgrounds_heatmap), 1,  figsize = (6,35))

for i in np.arange(0,len(backgrounds_heatmap)):
    mean=norm_mean_values_pd.loc[gene,backgrounds_heatmap[i]]
    std=norm_std_values_pd.loc[gene,backgrounds_heatmap[i]]
    sns.histplot(np.random.normal(mean,std,len(norm_mean_values_pd)),bins=200,
    color="gray",stat="percent",ax=axes[i],label=backgrounds_heatmap[i],kde=True,element="step")

    axes[i].set_xlabel(backgrounds_heatmap[i],fontsize=16)
    axes[0].set_title(list_data_pd_wt.loc[gene,"Gene name"],fontsize=16)
    axes[i].tick_params(axis="both",labelsize=16)
    axes[i].legend(fontsize=16)



# +
## Heatmaps 

backgrounds_heatmap=["wt_a","wt_b",'wt_merged', 'dnrp1_a','dnrp1_b','dnrp1_merged', 'dbem3_a',
'dbem3_b','dbem3_merged','bem1-aid_a', 'bem1-aid_b', "bem1-aid_merged",'dbem1dbem3_a', 'dbem1dbem3_b']
norm_mean_values_pd=norm_mean_values_pd.loc[:,backgrounds_heatmap]

list_data_pd_wt=list_data_pd.loc["wt_merged"]

index_polarity_genes=[]
chosen_polarity_genes=[]
for i in np.arange(0,len(polarity_genes)):
    tmp=np.where(list_data_pd_wt.loc[:,"Gene name"]==polarity_genes.index[i])
    if tmp[0].size!=0:
        index_polarity_genes.append(tmp[0][0])
        chosen_polarity_genes.append(polarity_genes.index[i])
        


# +
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 15))

sns.heatmap(data=norm_mean_values_pd.loc[index_polarity_genes],cmap="seismic",vmin=0,vmax=1,xticklabels=backgrounds_heatmap,
yticklabels=chosen_polarity_genes ,cbar=True,annot=True,cbar_kws={'label': 'Fitness compared to WT'})

labels=[]
for i in np.arange(0,len(chosen_polarity_genes)):
    labels.append(chosen_polarity_genes[i]+"-"+str(polarity_genes.loc[chosen_polarity_genes[i],:].unique()))
ax.set_yticklabels(labels);
# -

g=sns.clustermap(data=norm_mean_values_pd.loc[index_polarity_genes],cmap="seismic",vmin=0,vmax=1,xticklabels=backgrounds_heatmap,
yticklabels=chosen_polarity_genes ,cbar=True,annot=True,cbar_kws={'label': 'Fitness compared to WT'})
g.fig.set_size_inches((15,15))



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


a=data_norm_pd.loc["wt_merged"]["tr-density"].tolist()
sorted_a=sorted(a,reverse=True)
sorted_a[0:10]

a=np.where(data_norm_pd.loc["wt_merged"]["tr-density"]==sorted_a[2900])
b=data_norm_pd.loc["wt_merged"]
a[0][0],sorted_a[2900]

b.loc[:,"Gene name"][a[0][0]]


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

# +
# export to excel 
# rates2excel=rates_pd.copy()
# rates2excel.index=list_data_pd.loc["wt_merged"]["Gene name"]

# rates2excel.to_excel("../postprocessed-data/rates_non_coarsed_fitness_model.xlsx")
# -

list_data_pd_wt=list_data_pd.loc["wt_merged"]
gene=np.where(list_data_pd_wt.loc[:,"Gene name"]=="BEM1")[0][0]
rates_pd.loc[gene,"wt_merged"]

# +
from scipy import stats

mean_values=np.zeros((len(rates_pd),len(keys)))
std_values=np.zeros((len(rates_pd),len(keys)))
sem_values=np.zeros((len(rates_pd),len(keys)))

for i in np.arange(0,len(keys)):
    for j in np.arange(0,len(rates_pd)):
        gene_values=rates_pd.loc[j , keys[i]]
        mean_values[j,i]=(np.mean(gene_values))
        std_values[j,i]=(np.std(gene_values))
        
        
        if type(gene_values)!=int: 
            sem_values[j,i]=std_values[j,i]/np.sqrt(len(gene_values))
        else:
            sem_values[j,i]=0



# +
# values2excel_pd=pd.DataFrame(mean_values,columns=keys)
# values2excel_pd=pd.DataFrame(std_values,columns=keys)
# values2excel_pd=pd.DataFrame(sem_values,columns=keys)

# values2excel_pd.to_excel("../postprocessed-data/sem_values_non_coarsed_fitness_model.xlsx")

# +
index_wt_merged=np.where(np.array(keys)=="wt_merged")[0][0]

norm_mean_values=[]
norm_std_values=[]
norm_sem_values=[]

for i in np.arange(0,len(mean_values)):
    norm_mean_values.append(mean_values[i]/mean_values[i,index_wt_merged])
    norm_std_values.append(std_values[i]/std_values[i,index_wt_merged]) # fix this , we should propagate the error
    norm_sem_values.append(sem_values[i]/sem_values[i,index_wt_merged]) # fix this , we should propagate the error
   

norm_mean_values=np.array(norm_mean_values)
norm_std_values=np.array(norm_std_values)
norm_sem_values=np.array(norm_sem_values)

norm_mean_values[~np.isfinite(norm_mean_values)] = 0
norm_std_values[~np.isfinite(norm_std_values)] = 0
norm_sem_values[~np.isfinite(norm_sem_values)] = 0

norm_mean_values_pd=pd.DataFrame(norm_mean_values,columns=keys)
norm_std_values_pd=pd.DataFrame(norm_std_values,columns=keys)
norm_sem_values_pd=pd.DataFrame(norm_sem_values,columns=keys)

# +
list_data_pd_wt=list_data_pd.loc["wt_merged"]
index_wt_merged=np.where(np.array(keys)=="wt_merged")[0][0]
ho=np.where(list_data_pd_wt.loc[:,"Gene name"]=="HO")[0][0]

fitness_wt_ho=mean_values[ho,index_wt_merged]
std_wt_ho=std_values[ho,index_wt_merged]
sem_wt_ho=sem_values[ho,index_wt_merged]

mean_values_wt=mean_values[:,index_wt_merged]
mean_values_wt[~np.isfinite(mean_values_wt)] = 0

std_values_wt=std_values[:,index_wt_merged]
std_values_wt[~np.isfinite(std_values_wt)] = 0

values=mean_values[:,index_wt_merged]/fitness_wt_ho
# values=values[np.where(values!=-np.inf)]
# values=values[np.where(values!=np.nan)]
# values=values[np.where(values!=np.inf)]
values[~np.isfinite(values)] = 0

std_values_wt_2HO=values*np.sqrt((np.array(std_values_wt)/np.array(mean_values_wt))**2+(np.array(std_wt_ho)/fitness_wt_ho)**2)
# std_values_wt=std_values_wt[np.where(std_values_wt!=-np.inf)]
# std_values_wt=std_values_wt[np.where(std_values_wt!=np.nan)]
# std_values_wt=std_values_wt[np.where(std_values_wt!=np.inf)]
std_values_wt[~np.isfinite(std_values_wt)] = 0

sem_values_wt=sem_values[:,index_wt_merged]
sem_values_wt[~np.isfinite(sem_values_wt)] = 0

# -

np.mean(std_values_wt)/np.mean(mean_values_wt)

# +
fig, axes = plt.subplots(2, 1,  gridspec_kw={"height_ratios":(.10, .30)}, figsize = (8, 5))
sns.violinplot(values, ax=axes[0],color="gray",orient="h",inner="quartile")
#sns.violinplot(sem_values_wt, ax=axes[0],color="pink",orient="h",inner="quartile")
sns.histplot(values,bins=200,color="gray",ax=axes[1],stat="percent",label="mean-WT",kde=True,element="step")
#sns.histplot(sem_values_wt,bins=200,color="pink",ax=axes[1],stat="percent",label="std-WT",kde=True,element="step")
axes[1].set_xlabel("Fitness values compared to HO locus",fontsize=16)
axes[1].set_ylabel("Percent",fontsize=16)
axes[1].tick_params(axis="both",labelsize=16)
axes[0].tick_params(axis="both",labelsize=16)
axes[1].legend(fontsize=16)

#fig.savefig("../figures/figures_thesis_chapter_2/fig_fitness_distribution_normalized_wt2HO_non_coarse_model_mean.png",dpi=400)

# +
fig=plt.figure(figsize=(8,5))



plt.hist2d(mean_values_wt,sem_values_wt,density=True,bins=(20,60), cmap=plt.cm.Blues);
cbar = plt.colorbar()
for t in cbar.ax.get_yticklabels():
     t.set_fontsize(16)
plt.xlabel("Mean of fitness values ",fontsize=16)
plt.ylabel("Sem  ",fontsize=16)
plt.xticks(fontsize=16);
plt.yticks(fontsize=16);
plt.xlim(0,3e-10)
plt.ylim(0,1e-10)

plt.title("2D histogram : SEM vs MEAN",fontsize=18)
plt.tight_layout()
#plt.savefig("../figures/figures_thesis_chapter_2/fig_sem_mean_fitness_distribution_normalized_wt2HO_non_coarse_model.png",dpi=400)

# +
fig=plt.figure(figsize=(8,5))

plt.hist2d(mean_values_wt,std_values_wt,density=True,bins=(20,60), cmap=plt.cm.Blues);
cbar = plt.colorbar()
for t in cbar.ax.get_yticklabels():
     t.set_fontsize(16)
plt.xlabel("Mean of fitness values",fontsize=16)
plt.ylabel("STD ",fontsize=16)
plt.xticks(fontsize=16);
plt.yticks(fontsize=16);
plt.xlim(0,3e-10)
plt.ylim(0,2e-10)
plt.title("2D histogram : STD vs MEAN",fontsize=18)
plt.tight_layout()
#plt.savefig("../figures/figures_thesis_chapter_2/fig_std_mean_fitness_distribution_normalized_wt2HO_non_coarse_model.png",dpi=400)

# +
high=[1.1,np.max(values)]
low=[high[0]*0.5,high[0]*0.7]
neutral=[low[1],high[0]]
essential=low[0]

high,low,neutral,essential

# +
## make a pie chart with the mean fitness values compared to HO

high=[1.1,np.max(values)]
low=[high[0]*0.5,high[0]*0.7]
neutral=[low[1],high[0]]
essential=low[0]



values_neutral=values[np.where((values>neutral[0]) & (values<neutral[1]))]
values_high=values[np.where((values>high[0]) & (values<high[1]))]
values_low=values[np.where((values<low[1]) & (values>low[0]))]
values_essential=values[np.where(values<essential)]

fig,axes=plt.subplots(1,1,figsize=(8,8))
#define data
data = [len(values_neutral),len(values_high),len(values_low),len(values_essential)]
labels = ['Neutral', 'Advantageous', 'Disadvantageous', 'Possible essentials']

#define Seaborn color palette to use
#colors = sns.color_palette('pastel')[0:4]
colors=["blue","green","red","pink"]

#create pie chart
plt.pie(data, labels = labels, colors = colors, autopct='%.0f%%',textprops={'fontsize': 18});
plt.tight_layout()
#fig.savefig("../figures/figures_thesis_chapter_2/fig_pie_fitness_distribution_normalized_wt2HO_non_coarse_model.png",dpi=400)


# +
## Mean fitness for essential genes

## Distribution of fitness values for annotated essential genes
genes=list_data_pd_wt.loc[:,"Gene name"]
genes.reset_index(drop=True,inplace=True)
values_essentials=[]


for i in standard_essentials:
    x=np.where(genes==i)[0]
    values_essentials.append(values[x])
   

values_essentials=np.array(values_essentials)

# values_essentials_below_low_fit=np.where(np.array(values_essentials)<np.array(low[0]))[0]
# values_essentials_above_low_fit=np.where(np.array(values_essentials)>np.array(low[0]) & np.array(values_essentials)< np.array(1))[0]

values_essentials_below_low_fit=values_essentials[(values_essentials<np.array(low[0]))]
values_essentials_above_low_fit=values_essentials[(values_essentials>np.array(low[0])) & (values_essentials< np.array(1))]
values_essentials_high=values_essentials[(values_essentials>np.array(1))]

a=len(values_essentials_below_low_fit)/len(standard_essentials)
b=len(values_essentials_above_low_fit)/len(standard_essentials)
c=len(values_essentials_high)/len(standard_essentials)



# +
fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

axes.hist(np.concatenate(values_essentials),bins=30,alpha=0.5,color="gray");
axes.vlines(1,0,150,color="red",linestyle="dashed",linewidth=2,label="HO");
axes.vlines(low[0],0,150,color="blue",linestyle="dashed",linewidth=2,label="Low fitness threshold")
axes.annotate(f"{a*100:.2f}" +"%",xy=(0.1,150),fontsize=16)
axes.annotate(f"{b*100:.2f}" +"%",xy=(0.6,150),fontsize=16)
axes.annotate(f"{c*100:.2f}" +"%",xy=(1.1,150),fontsize=16)
axes.set_title("Distribution of fitness values for annotated essential genes",fontsize=16)
axes.tick_params(axis="both",labelsize=16)
axes.set_xlabel("Fitness values compared to HO locus",fontsize=16)
axes.set_ylabel("Count",fontsize=16)
axes.legend(loc="best")
plt.tight_layout()

#fig.savefig("../figures/figures_thesis_chapter_2/fig_non_coarse_fitness_distribution_normalized_wt2HO_essentials.png",dpi=400)

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
## Comparison between technical replicates

replicate_a=["wt_a"]
replicate_b=["wt_b"]

data_replicate_a=norm_mean_values_pd.loc[:,replicate_a]
data_replicate_b=norm_mean_values_pd.loc[:,replicate_b]

values_diff=np.array(data_replicate_a[replicate_a[0]])-np.array(data_replicate_b[replicate_b[0]])

fig, axes = plt.subplots(3, 1,  gridspec_kw={"height_ratios":(.10, .40,.40)}, figsize = (8, 8))


sns.violinplot(np.abs(values_diff), ax=axes[0],color="gray",orient="h",inner="quartile")

g=sns.histplot(np.abs(values_diff),bins=200,color="gray",ax=axes[1],stat="percent",label="wt_a-wt_b",kde=True,element="step")

axes[1].set_xlabel("Fitness values of " + replicate_a[0] +"-" + replicate_b[0]+ " compared to wt",fontsize=16)
axes[1].set_ylabel("Percent",fontsize=16)


sns.regplot(data_replicate_a,data_replicate_b,ax=axes[2],scatter_kws={"s":30,"color":"gray"})

axes[2].set_xlabel("Fitness values of "+replicate_a[0],fontsize=16)
axes[2].set_ylabel("Fitness values of "+ replicate_a[0],fontsize=16)
axes[2].set_xlim(0,1)
axes[2].set_ylim(0,1)

for axes in axes:
    axes.tick_params(axis="both",labelsize=16)

plt.tight_layout()
# -

np.array(data_replicate_a["wt_a"])

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
g.savefig("../figures/fig_fitness_heatmap_normalized2WT_non_coarse_model.png",dpi=400)



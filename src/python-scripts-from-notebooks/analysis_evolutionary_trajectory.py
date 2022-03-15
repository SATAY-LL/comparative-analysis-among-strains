# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
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


# +
## Plot the number of genes with less than 5 transposons in all mutants 

pergene_files=[]
#data_dir= "../satay/data_files/data_unmerged/"
#data_dir="../transposonmapper/data_files/files4test/"
data_dir="../postprocessed-data/"
#data_dir="../transposonmapper/data_files/"
for root, dirs, files in os.walk(data_dir):
    for file in files:
        if file.endswith("pergene_insertions.xlsx"):
            pergene_files.append(os.path.join(root, file))
# -

list_data=[]
for i in pergene_files:
    list_data.append(pd.read_excel(i,engine='openpyxl',index_col="Unnamed: 0"))

keys=[]
for i in np.arange(0,len(pergene_files)):
    keys.append(pergene_files[i].split("/")[-1].split("_")[0]+"_"+pergene_files[i].split("/")[-1].split("_")[1])

keys

backgrounds=["wt_a","wt_b",'wt_merged', 'bem1-aid_a', 'bem1-aid_b', "bem1-aid_merged",
 'dbem1dbem3_a', 'dbem1dbem3_b','dnrp1_merged', 'dbem3_merged']

backgrounds=keys

list_data_pd=pd.concat(list_data,axis=0,keys=keys)

list_data_pd=list_data_pd.loc[backgrounds]

# +
# import excel file with the normalized data

data_norm_pd=pd.read_excel("../postprocessed-data/data_norm_linear_transformation_per_background.xlsx",
engine='openpyxl',index_col="background")
data_norm_pd.drop(columns=["Unnamed: 0","Unnamed: 1"],inplace=True)
# -

data_norm_pd.index.unique()

polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")
polarity_genes.fillna(0,inplace=True)

data_norm_pd_wt=data_norm_pd.loc["wt_merged"]
data_norm_pd_wt.columns
x_norm=data_norm_pd_wt.loc[:,"reads_normalized_windows"]/data_norm_pd_wt.loc[:,"tr_normalized_windows"]
x=data_norm_pd_wt.loc[:,"Reads"]/data_norm_pd_wt.loc[:,"Insertions"]
# plt.subplot(1,2,1)
# plt.hist(x_norm,bins=100);
# plt.subplot(1,2,2)
# plt.hist(x,bins=100);
x_norm.std()/x_norm.mean(),x.std()/x.mean()

# +
# Number of genes with less than 2 transposons and 2 reads in all mutants
Q=[]
for i in np.arange(0,len(backgrounds)):
    tmp=(list_data_pd.loc[backgrounds[i]])
    #L.append(len(tmp[(tmp.loc[:,"Insertions"]<5) & (tmp.loc[:,"Reads"]<2)]))
    Q.append(len(tmp[(tmp.loc[:,"Insertions"]<5)]))



# +
Q_dict=dict(zip(backgrounds,Q))

A={k: v for k, v in sorted(Q_dict.items(), key=lambda item: item[1])}

fig = plt.figure(figsize=(10, 7))
plt.bar(A.keys(),A.values())
plt.xticks(rotation=90);

plt.ylabel("# of genes with less than 5 insertions ",fontsize=9)
plt.tight_layout(pad=3)
#plt.savefig("../figures/fig_number_genes_with_less_than_5_reads.png",dpi=300)


# +
# Number of transposons in each library 
L=[]
R=[]
for i in np.arange(0,len(backgrounds)):
    tmp=(list_data_pd.loc[backgrounds[i]])
    L.append(np.sum(tmp.loc[:,"Insertions"]))
    R.append(np.sum(tmp.loc[:,"Reads"]))

    
L_dict=dict(zip(backgrounds,L))
R_dict=dict(zip(backgrounds,R))

A={k: v for k, v in sorted(L_dict.items(), key=lambda item: item[1])}
B={k: v for k, v in sorted(R_dict.items(), key=lambda item: item[1])}

fig,axes = plt.subplots(nrows=1,ncols=2,figsize=(15, 7))
axes[0].bar(A.keys(),A.values())
axes[0].set_xticklabels(labels=A.keys(),rotation=90);
axes[0].set_ylabel("# transposons per library",fontsize=9)

axes[1].bar(B.keys(),B.values())
axes[1].set_xticklabels(labels=B.keys(),rotation=90);
axes[1].set_ylabel("# reads per library",fontsize=9)

#plt.savefig("../figures/fig_number_insertions_and_reads_for_all_libraries.png",dpi=300)

plt.tight_layout(pad=3)

# +
## Plot length of genes vs number of transposons
from scipy.optimize import curve_fit

def func(x, a,b):
    return a*x +b

fig,axes = plt.subplots(nrows=1,ncols=1,figsize=(10, 5))

data=list_data_pd.loc["wt_merged"]
xaxis=(data.loc[:,"End location"]-data.loc[:,"Start location"])
sns.regplot(x=xaxis,y=data.loc[:,"Insertions"],ax=axes,color="black",scatter_kws={"s":80,"alpha":0.2},label="WT")

# plt.plot(xaxis,data.loc[:,"Insertions"],"o",markersize=10,alpha=0.5)

# x = np.array(xaxis)
# popt, pcov = curve_fit(func, x, np.array(data.loc[:,"Insertions"]))
# plt.plot(x, func(x, *popt), 'ko', label="Fitted Curve->" + str(popt[0]) + '*x+'+str(popt[1]),alpha=0.6)

plt.ylabel("gene length",fontsize=18)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("# transposons",fontsize=18)
plt.legend()
#plt.savefig("../figures/fig_length_vs_number_of_transposons.png",dpi=300)

# +
## Plot the number of genes with less than 5 transposons in all mutants over the total number of transposons per library
from scipy.optimize import curve_fit

fig,axes = plt.subplots(nrows=1,ncols=1,figsize=(15, 10))
l_array=[]
q_array=[]
for keys in L_dict.keys():
    axes.plot(L_dict[keys],Q_dict[keys],'*r',markersize=20,alpha=0.5)
    axes.annotate(keys,(L_dict[keys],Q_dict[keys]))
    l_array.append(L_dict[keys])
    q_array.append(Q_dict[keys])

def func(x, a,b):
    return b+a / x

x = np.array(l_array)

popt, pcov = curve_fit(func, x, np.array(q_array))
axes.plot(x, func(x, *popt), 'ko', label="Fitted Curve (b+a/x)",alpha=0.6,markersize=10)

axes.set_xlabel("# of transposons in that library",fontsize=16)
axes.set_title("# of genes with less than 5 transposons",fontdict={"size":16})
axes.legend()
axes.set_xscale("log")
axes.set_yscale("log")
plt.grid(axis="both",which="both",ls="--",alpha=0.5)
plt.tick_params(axis='both', which='both', labelsize=16)

#fig.savefig("../figures/number_of_tr_less_than_5_vs_total_number_of_tr.png",dpi=300)

# +
from module_intergenic_model import adding_features2dataframe,getting_r

list_data_extended=[]
for i in np.arange(0,len(backgrounds)):
    tmp=(list_data_pd.loc[backgrounds[i]])
    list_data_extended.append(adding_features2dataframe(tmp))


# -

list_data_extended_pd=pd.concat(list_data_extended,axis=0,keys=backgrounds)

# +
# Number of reads per transposons per library
L=[]
R=[]
for i in np.arange(0,len(backgrounds)):
    tmp=(list_data_extended_pd.loc[backgrounds[i]])
    L.append(tmp.loc[:,"reads-per-tr"].mean())
    R.append(tmp.loc[:,"reads-per-tr"].std())

L_dict=dict(zip(backgrounds,L))
R_dict=dict(zip(backgrounds,R))

A={k: v for k, v in sorted(L_dict.items(), key=lambda item: item[1])}
B={k: v for k, v in sorted(R_dict.items(), key=lambda item: item[1])}

fig,axes = plt.subplots(nrows=1,ncols=2,figsize=(15, 7))
axes[0].bar(A.keys(),A.values())
axes[0].set_xticklabels(labels=A.keys(),rotation=90);
axes[0].set_ylabel("Mean number of reads per transposons per library",fontsize=9)

axes[1].bar(B.keys(),B.values())
axes[1].set_xticklabels(labels=B.keys(),rotation=90);
axes[1].set_ylabel("Standard deviation of number of reads per transposons per library",fontsize=9)

plt.tight_layout(pad=3)
#plt.savefig("../figures/fig_mean_and_std_number_reads_per_transposons_per_library.png",dpi=300)
# -

# ## Compute fitness from a coarse grained model, taking the total sum of reads per transposon per gene 

# +

list_data_rates=[]
for i in np.arange(0,len(backgrounds)):
    tmp=(list_data_extended_pd.loc[backgrounds[i]])
    list_data_rates.append(getting_r(tmp))


# -

# ## Export fitness values to excel to be used in other notebooks 

# +
# fitness_all=[]
# for i in np.arange(0,len(backgrounds)):

#     tmp=list_data_rates[i][0].to_frame()
#     tmp.columns=["fitness"]
#     tmp.index=list_data_pd.loc[backgrounds[i],"Gene name"]
#     fitness_all.append(tmp)

# fitness_all_pd=pd.concat(fitness_all,axis=0,keys=backgrounds)
# fitness_all_pd.reset_index(inplace=True)

# fitness_all_pd.rename(columns={"level_0":"background"},inplace=True)
# fitness_all_pd.to_excel("../postprocessed-data/fitness_coarse_grained_all_pd.xlsx")


# +
rates_dict=dict(zip(backgrounds,list_data_rates))

rates_norm_dict=defaultdict(dict)

for i in np.arange(0,len(backgrounds)):
    tmp=(rates_dict[backgrounds[i]])
    if rates_dict["wt_merged"]!=0:
        rates_norm_dict[backgrounds[i]]["rates-intergenic"]=np.divide(tmp,rates_dict["wt_merged"])[0]
    else:
        rates_norm_dict[backgrounds[i]]["rates-intergenic"]=tmp[0]
# -

rates_norm_dict_pd=pd.DataFrame.from_dict(rates_norm_dict,orient="index")
rates_norm_dict_pd.head(2)

# +
## Fitness plots  normalized to the values of wt_merged

keys_fitness=["bem1-aid_a","bem1-aid_b"]
#keys_fitness=backgrounds

rates_norm_dict_pd=pd.DataFrame.from_dict(rates_norm_dict,orient="index")

plt.subplots(nrows=1,ncols=2,figsize=(8, 3))
j=1
for i in np.arange(0,len(keys_fitness),2):
    plt.subplot(1,2,j)
    plt.scatter(rates_norm_dict_pd.loc[keys_fitness[i]].tolist(),
    rates_norm_dict_pd.loc[keys_fitness[i+1]].tolist(),s=30,alpha=0.5)
    plt.plot(np.arange(0,2.1,0.1),np.arange(0,2.1,0.1),color="black")
   
    plt.ylabel("Fitness_" + keys_fitness[i+1],fontsize=9)
    plt.xlabel("Fitness_" + keys_fitness[i],fontsize=9)
    plt.tight_layout(pad=3)
    
    j=j+1
#plt.savefig("../figures/fig_fitness_scatter_normalized_wt_merged.png",dpi=300)
# -

backgrounds

# +
## Distribution of fitness values from WT

list_data_pd_wt=list_data_pd.loc["wt_merged"]
ho=np.where(list_data_pd_wt.loc[:,"Gene name"]=="HO")[0][0]

index_wt_merged=np.where(np.array(backgrounds)=="wt_merged")[0][0]

fitness_wt=list_data_rates[index_wt_merged]
fitness_wt_ho=fitness_wt[0][ho]

values=np.array(fitness_wt[0])/fitness_wt_ho

genes=list_data_pd_wt.loc[np.where(values!=-np.inf)[0],"Gene name"]
values=values[np.where(values!=-np.inf)]

# -

low_fit=np.where(values<0.5)[0]
len(low_fit)
genes.reset_index(drop=True,inplace=True)
x=np.where(genes=="CDC24")[0]
values[x]

# +
from statistics import NormalDist

def confidence_interval(data, confidence=0.95):
  dist = NormalDist.from_samples(data)
  z = NormalDist().inv_cdf((1 + confidence) / 2.)
  h = dist.stdev * z / ((len(data) - 1) ** .5)
  return dist.mean - h, dist.mean + h

values=values[np.where(values!=-np.inf)]

a,b=confidence_interval(values)
# -

a,b,np.mean(values)

# +



fig, axes = plt.subplots(2, 1,  gridspec_kw={"height_ratios":(.10, .30)}, figsize = (8, 5))


sns.violinplot(values, ax=axes[0],color="gray",orient="h",inner="quartile")

g=sns.histplot(values,bins=200,color="gray",ax=axes[1],stat="percent",label="WT",kde=True,element="step")

axes[1].set_xlabel("Fitness values compared to HO locus",fontsize=16)
axes[1].set_ylabel("Percent",fontsize=16)
axes[1].tick_params(axis="both",labelsize=16)
axes[0].tick_params(axis="both",labelsize=16)

#axes[1].vlines(x=1,ymin=0,ymax=2.5,color="red",linestyle="--",linewidth=2,label="HO")
#axes[1].vlines(x=np.mean(values),ymin=0,ymax=2.5,color="black",linestyle="--",linewidth=1)

axes[1].legend()
#fig.savefig("../figures/figures_thesis_chapter_2/fig_fitness_distribution_normalized_wt2HO.png",dpi=400)

# +
half_maximun=np.max(g.get_lines()[0].get_data()[1])/2

index_half_maximun=[0,41]

x0,x1=g.get_lines()[0].get_data()[0][0],g.get_lines()[0].get_data()[0][41]
xf=g.get_lines()[0].get_data()[0][-1]

print("The full width at half maximum is:",x1-x0, 
"which represents the", (x1-x0)/xf*100,"% of the total width")

# +
## heatmap of fitnes values across backgrounds
list_data_pd_wt=list_data_pd.loc["wt_merged"]

index_polarity_genes=[]
chosen_polarity_genes=[]
for i in np.arange(0,len(polarity_genes)):
    tmp=np.where(list_data_pd_wt.loc[:,"Gene name"]==polarity_genes.index[i])
    if tmp[0].size!=0:
        index_polarity_genes.append(tmp[0][0])
        chosen_polarity_genes.append(polarity_genes.index[i])
# -

backgrounds_heatmap=["wt_a","wt_b",'wt_merged', 'dnrp1_a','dnrp1_b','dnrp1_merged', 'dbem3_a',
'dbem3_b','dbem3_merged','bem1-aid_a', 'bem1-aid_b', "bem1-aid_merged",'dbem1dbem3_a', 'dbem1dbem3_b']

# +
array2heatmap=np.zeros((len(chosen_polarity_genes),len(backgrounds_heatmap)))

for i in np.arange(0,len(backgrounds_heatmap)):
    for j in np.arange(0,len(chosen_polarity_genes)):
        array2heatmap[j,i]=rates_norm_dict_pd.loc[backgrounds_heatmap[i]].tolist()[0][index_polarity_genes[j]]
  

# +
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 15))
array2heatmap[~np.isfinite(array2heatmap)] = 0
sns.heatmap(array2heatmap,cmap="seismic",vmin=0,vmax=1,xticklabels=backgrounds_heatmap,
yticklabels=chosen_polarity_genes ,cbar=True,annot=True,cbar_kws={'label': 'Fitness compared to WT'})

labels=[]
for i in np.arange(0,len(chosen_polarity_genes)):
    labels.append(chosen_polarity_genes[i]+"-"+str(polarity_genes.loc[chosen_polarity_genes[i],:].unique()))
ax.set_yticklabels(labels);
fig.savefig("../figures/fig_heatmap_fitness_normalized_wt_with_replicates-and-function.png",dpi=300)
# -

g=sns.clustermap(array2heatmap,cmap="seismic",vmin=0,vmax=1,xticklabels=backgrounds_heatmap,
yticklabels=chosen_polarity_genes ,cbar=True,annot=True,cbar_kws={'label': 'Fitness compared to WT'})
g.fig.set_size_inches((15,15))
g.savefig("../figures/fig_heatmap_fitness_standarized_clustermap.png",dpi=300)

# +
import math

def sigmoid(x,k,n):
    a = []
    for item in x:
        a.append(n/(n+math.exp(-item/k)))
    return a

import matplotlib.pyplot as plt
import numpy as np

x = np.arange(-8, 8, 0.2)
for k in np.arange(1,2,0.2):
    plt.plot(x, sigmoid(x,k,1), color='black', linewidth=1.0, linestyle='-')
    plt.xticks([])
    plt.yticks([])

plt.savefig("../figures/figures_thesis_chapter_2/fig_sigmoid_function.png",dpi=400)


# +
# from from_excel_to_list import from_excel_to_list
# list_data_extended_pd=pd.concat(list_data_extended,axis=0,keys=keys)
# list_data_rates_pd=list_data_extended_pd.copy()
# for i in np.arange(0,len(keys)):
#     list_data_rates_pd.loc[keys[i],"rates-intergenic"]=rates_norm_dict[keys[i]]["rates-intergenic"]

# # Convert reads per location to numeric arrays
# list_data_rates_pd.loc[keys[0],"Reads per insertion location"][0]=from_excel_to_list(list_data_rates_pd.loc[keys[0],"Reads per insertion location"][0])

# TO DO:
# - Replace the column "Reads per insertion location" with a numeric array
# - Replace the colum "Insertion locations" with a numeric array
# - Implement the uncertainty of the rates

# +
#%% uncertainty related to the rates 
## expand this expression : np.log(N*K/(K-N))/T taking into account that N=p/q
## uncertainty of p= std(p) per gene and uncertainty of q= 1-density(q) 


# data=data_wt_extended
# T=90
# #cte=T*np.sum(data["reads-per-tr"])
# cte=np.sum(data["reads-per-tr"])
# P=data["reads-per-tr"]/cte
# uncertainty_r=[]
# for i in data.index:
#     p=data["Reads per insertion location"][i]
#     q=data["tr-density"][i]
#     uncertainty_r.append(P[i]*np.sqrt((np.std(p)/data["reads-per-tr"][i])**2+((1-q)/q)**2))

# data["uncertainty_r"]=uncertainty_r

# +
# Implement the std of the fitness values



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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys
from collections import defaultdict
from ast import literal_eval
from scipy.stats import norm

from from_excel_to_list import from_excel_to_list

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

keys= ['wt_merged','wt_a','wt_b','dnrp1_1','dnrp1_2']

data=list_data_pd.loc[keys] # Take only data from targeted genotypes


data.columns

# +
from sklearn.linear_model import LinearRegression

#initiate linear regression model
model = LinearRegression()

#define predictor and response variables



figure,ax=plt.subplots(1,2,figsize=(10,5))
plt.subplots_adjust(wspace=0.4,hspace=0.5)

magnitudes=["Reads","Insertions"]
strains=["wt_a","wt_b","dnrp1_1","dnrp1_2"]

for i in np.arange(0,len(magnitudes)):
    tmp_0=data.loc[strains[0],magnitudes[i]]/data.loc[strains[0],magnitudes[i]].sum()
    tmp_1=data.loc[strains[1],magnitudes[i]]/data.loc[strains[1],magnitudes[i]].sum()

    X, y = tmp_0.values.reshape(-1,1), tmp_1.values.reshape(-1,1)

    #fit regression model
    model.fit(X, y)

    #calculate R-squared of regression model
    r_squared = model.score(X, y)


    ax[i].scatter(tmp_0,tmp_1,s=20,alpha=0.5)
    ax[i].set_title("Normalized "+magnitudes[i])
    ax[i].set_xlabel("tech replicate 1")
    ax[i].set_ylabel("tech replicate 2")
    ax[i].set_xlim(0,0.001)
    ax[i].set_ylim(0,0.001)
    ax[i].plot([0,0.001],[0,0.001],color="black",linestyle="--")
    ax[i].text(0, 0.0005, '$R^2=%.3f$' % (r_squared),fontsize=12)

    # ax[i,1].set_title( magnitudes[i] + " difference")
    # ax[i,1].hist(tmp_1-tmp_0,bins=100)
    # data2fit = tmp_1-tmp_0
    # mu, std = norm.fit(data2fit)
    
    # xmin, xmax = ax[i,1].get_xlim()
    # x = np.linspace(xmin, xmax, 100)
    # p = norm.pdf(x, mu, std)

    # ax[i,1].plot(x, p, 'k', linewidth=2)

   
figure,ax=plt.subplots(1,2,figsize=(10,5))
plt.subplots_adjust(wspace=0.4,hspace=0.5)

for i in np.arange(0,len(magnitudes)):
    tmp_0=data.loc[strains[2],magnitudes[i]]/data.loc[strains[2],magnitudes[i]].sum()
    tmp_1=data.loc[strains[3],magnitudes[i]]/data.loc[strains[3],magnitudes[i]].sum()

    X, y = tmp_0.values.reshape(-1,1), tmp_1.values.reshape(-1,1)

    #fit regression model
    model.fit(X, y)

    #calculate R-squared of regression model
    r_squared = model.score(X, y)


    ax[i].scatter(tmp_0,tmp_1,s=20,alpha=0.5)
    ax[i].set_title("Normalized "+magnitudes[i])
    ax[i].set_xlabel("biolog replicate 1")
    ax[i].set_ylabel("biolog replicate 2")
    ax[i].set_xlim(0,0.001)
    ax[i].set_ylim(0,0.001)
    ax[i].plot([0,0.001],[0,0.001],color="black",linestyle="--")
    ax[i].text(0, 0.0005, '$R^2=%.3f$' % (r_squared),fontsize=12)


#figure.savefig("../figures/fig_differences_WT_technical_replicates.png",dpi=300)



# +
data_post=[]

for i in keys:
    data_post.append(pd.read_excel(data_dir+i+".xlsx",index_col="Unnamed: 0",engine='openpyxl'))
    

list_data_post=pd.concat(data_post,axis=0,keys=keys)
list_data_post.drop(columns=["Feature_alias","Feature_name"],inplace=True)
list_data_post.fillna(1,inplace=True)


# -

def romanize(number):
   n2rMap = {1000:'M', 900:'CM', 500:'D', 400:'CD', 100:'C', 90:'XC', 50:'L', 40:'XL', 10:'X', 9:'IX', 5:'V', 4:'IV', 1:'I'}
   roman = ""
   for key in n2rMap.keys():
      count = int(number / key)
      roman += n2rMap[key] * count
      number -= key * count
   return roman


chrom=[]
for i in np.arange(1,17):
    chrom.append(romanize(i))
    

data=list_data_post.loc["wt_merged"]
data.columns

# +
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
data=list_data_post.loc["wt_merged"]
plt.plot(data.index,data["Ninsertions"],color="purple",label="WT",alpha=0.5)

coordinates_chrom=[]

for i in chrom:
    coordinates_chrom.append((np.where(data["Chromosome"]==i)[0][0]))

plt.xticks(coordinates_chrom,chrom,fontsize=18,rotation=90);
plt.ylim(0,1000)

centromere_pos=np.where(list_data_post.loc["wt_merged"]["Feature_type"]=="Centromere")

for centromere in centromere_pos:
    plt.vlines(centromere,0,1000,color="black",linestyles="dashed",alpha=0.7,label="centromeres")

plt.ylabel("Transposon counts",fontsize=18)


# noncoding_pos=np.where(list_data_post.loc["wt_merged"]["Feature_type"]=="Noncoding_exon")

# plt.plot(noncoding_pos[0],data.loc[noncoding_pos,"Ninsertions"],alpha=0.5,
# color="purple",label="noncoding exons")

gene_pos=np.where(list_data_post.loc["wt_merged"]["Feature_type"]=="Gene; Verified")

plt.plot(gene_pos[0],data.loc[gene_pos,"Ninsertions"],alpha=0.5,color="green",label="CDS")
plt.legend(loc="upper right",fontsize=16)
plt.xlabel("Chromosomes",fontsize=18)
plt.yticks(fontsize=18);
plt.tight_layout()
# -

fig.savefig("../figures/figures_thesis_chapter_2/transposon_counts_pergene_wt_merged_with features.png",
dpi=400,transparent=True)


# +
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
data=list_data_post.loc["wt_merged"]
plt.plot(data.index,data["Nreads"],color="purple",label="WT",alpha=0.5)

coordinates_chrom=[]

for i in chrom:
    coordinates_chrom.append((np.where(data["Chromosome"]==i)[0][0]))

plt.xticks(coordinates_chrom,chrom,fontsize=18,rotation=90);


centromere_pos=np.where(list_data_post.loc["wt_merged"]["Feature_type"]=="Centromere")

for centromere in centromere_pos:
    plt.vlines(centromere,0,200000,color="black",linestyles="dashed",alpha=0.7,label="centromeres")




# noncoding_pos=np.where(list_data_post.loc["wt_merged"]["Feature_type"]=="Noncoding_exon")

# plt.plot(noncoding_pos[0],data.loc[noncoding_pos,"Nreads"],alpha=0.5,
# color="purple",label="noncoding exons")

gene_pos=np.where(list_data_post.loc["wt_merged"]["Feature_type"]=="Gene; Verified")

plt.plot(gene_pos[0],data.loc[gene_pos,"Nreads"],alpha=0.5,color="green",label="CDS")
plt.ylabel("Read counts",fontsize=18)
plt.xlabel("Chromosomes",fontsize=18)
plt.yticks(fontsize=18)
plt.legend(loc="upper right",fontsize=16)
plt.ylim(0,200000)
plt.tight_layout()
# -

fig.savefig("../figures/figures_thesis_chapter_2/read_counts_pergene_wt_merged_with features.png",
dpi=400,transparent=True)

## Reads distribution per chromosome:
fig,ax=plt.subplots(nrows=4,ncols=4,figsize=(20,16))
plt.subplots_adjust(hspace=0.5,wspace=0.5)
for i in np.arange(1,17):
    plt.subplot(4,4,i)
    #data[data.loc[:,"Chromosome"]==chrom[i-1]]["Nreads"].hist(bins=100,alpha=0.5,color="black",label=chrom[i-1])
    x=data[data.loc[:,"Chromosome"]==chrom[i-1]]["Nreads"]
    sns.histplot(x,label=chrom[i-1],
    stat="density",binwidth=200,kde=True,color="black",alpha=0.5)
    plt.xlabel("Reads",fontsize=18)
    plt.ylabel("Density",fontsize=18)
    plt.xticks(ticks= [min(x),max(x)],fontsize=16)
    plt.ylim(0,0.0008)
    plt.yticks(ticks= [0,0.0008],fontsize=16)
    plt.legend()
    # plt.tick_params(
    #     axis='x',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     bottom=False,      # ticks along the bottom edge are off
    #     top=False,         # ticks along the top edge are off
    #     labelbottom=False) # labels along the bottom edge are off

fig.savefig("../figures/figures_thesis_chapter_2/read_distribution_per_chromosome.png",
dpi=400,transparent=True)

# +
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(5,5))

sns.histplot(x="Nreads",data=data,stat="percent",binwidth=500,kde=False,color="black",
alpha=0.5)
plt.xlabel("Reads",fontsize=18)
plt.ylabel("Percent",fontsize=18)
plt.xlim(0,10000)
plt.ylim(0,50)
plt.yticks(ticks= [0,10,20,30,40],fontsize=16)# 
plt.xticks(ticks= [0,5000,10000],fontsize=16)

plt.title("Reads distribution WT genome",fontsize=18)
plt.tight_layout()

# +
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(5,5))

sns.histplot(x="Ninsertions",data=data,stat="percent",binwidth=10,kde=False,color="black",
alpha=0.5)
plt.xlabel("Insertions",fontsize=18)
plt.ylabel("Percent",fontsize=18)
plt.xlim(0,200)
#plt.ylim(0,50)
plt.yticks(ticks= [0,10,20,30],fontsize=16)# 
plt.xticks(ticks= [0,100,200],fontsize=16)

plt.title("Insertions distribution WT genome",fontsize=18)
plt.tight_layout()
# -

fig.savefig("../figures/figures_thesis_chapter_2/read_distribution_genome.png",
dpi=400,transparent=True)



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

backgrounds= ['wt_merged','bem1-aid_a','bem1-aid_b','dbem1dbem3_a','dbem1dbem3_b',
'dnrp1_merged','dbem3_merged']

chrom_length=pd.read_excel("../postprocessed-data/chromosome_length.xlsx",index_col="Chromosome")
chrom_length.drop(columns=["Unnamed: 0"],inplace=True)
chrom_length.loc["I"]


def linear_transformation_per_background(pergene_insertions_all_data,background,chrom_length,
windows_size=10000):
    """Executes the linear transformation procedure for a given background.
    It will divide the insertions and reads per gene over the total amount in each chromosome.

    Parameters
    ----------
    pergene_insertions_all_data : pandas.dataframe
        _description_
    background : list
        _description_
    chrom_length : pandas.dataframe
        _description_

    Returns
    -------
    _type_
        _description_
    """

    chrom=["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV",
    "XV","XVI"]

    data=pergene_insertions_all_data.copy()

    data_background=data.loc[background]
    

    reads_per_chrom=data_background.groupby(by=["Chromosome"]).sum()["Reads"]
    insertions_per_chrom=data_background.groupby(by=["Chromosome"]).sum()["Insertions"]

    data_background.index=data_background["Chromosome"]

    
    for i in chrom:

        data_background_chrom=data_background.loc[data_background["Chromosome"]==i]
        lengths_genes=data_background_chrom.loc[:,"End location"]-data_background_chrom.loc[:,"Start location"]

        data_background.loc[data_background["Chromosome"]==i,"Linear-norm-reads"]=data_background_chrom.loc[:,"Reads"]/reads_per_chrom[i]
        data_background.loc[data_background["Chromosome"]==i,"Linear-norm-insertions"]=data_background_chrom.loc[:,"Insertions"]/insertions_per_chrom[i]
        
        data_background.loc[data_background["Chromosome"]==i,"Linear-norm-tr-density"]=data_background.loc[data_background["Chromosome"]==i,"Linear-norm-insertions"]*np.divide(np.array(chrom_length.loc[i]),np.array(lengths_genes))

        data_background.loc[data_background["Chromosome"]==i,"tr-density"]=data_background_chrom.loc[:,"Insertions"]/lengths_genes
        data_background.loc[data_background["Chromosome"]==i,"length gene"]=lengths_genes
        
    data_without_index=data_background.copy()
    data_without_index.drop(columns=["Chromosome"],inplace=True)
    data_without_index.reset_index(inplace=True)

 #### Normalizing taking into account a windows size around the gene instead of the whole chromosome
    for gene in data_without_index["Gene name"].unique():
    
        location_gene=[data_without_index[data_without_index["Gene name"]==gene]["End location"].values[0],
                    data_without_index[data_without_index["Gene name"]==gene]["Start location"].values[0]]

       

        windows_location=[location_gene[0]-windows_size,location_gene[1]+windows_size]


        locations_ups=np.where((windows_location[0]<data_without_index["Start location"]) & (data_without_index["Start location"]<location_gene[0]))[0]
        locations_down=np.where((location_gene[1]<data_without_index["End location"]) & (data_without_index["End location"] < windows_location[1]))[0]

        locations_total=np.unique(np.concatenate((locations_ups,locations_down)))

        total_insertions=data_without_index.loc[locations_total,"Insertions"].sum()
        total_reads=data_without_index.loc[locations_total,"Reads"].sum()

        index_gene=np.where(data_without_index["Gene name"]==gene)[0][0]

        data_without_index.loc[index_gene,"tr_normalized_windows"]=data_without_index.loc[index_gene,"Insertions"]/total_insertions
        data_without_index.loc[index_gene,"reads_normalized_windows"]=data_without_index.loc[index_gene,"Reads"]/total_reads
        data_without_index.loc[index_gene,"insertions_over_windows"]=total_insertions
        data_without_index.loc[index_gene,"reads_over_windows"]=total_reads

    return data_without_index,reads_per_chrom,insertions_per_chrom
        

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
# data_norm_pd.to_excel("../postprocessed-data/data_norm_linear_transformation_all_backgrounds.xlsx")
# -

data_norm_pd.loc["wt_merged"][0:3]

reads_per_chrom_pd.loc["wt_merged"]["I"]

# comparison of libraries in terms of normalized transposon densityies 
mutant="dbem1dbem3_a"
plt.scatter(data_norm_pd.loc["wt_merged","Linear-norm-tr-density"],data_norm_pd.loc[mutant,"Linear-norm-tr-density"],color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("wt_merged")
plt.ylabel(mutant)
plt.ylim(0.1,100)
plt.xlim(0.1,100)

# +
## make a function that looks for the gene names that have more than 10 as transposon density
# in the mutant background and less than one in WT, and viceversa. 
# -

data_norm_pd.loc["wt_merged"].columns


# +
def linear_transformations_plots(normalized_data,type,background,saveFigure=False):
    """_summary_

    Parameters
    ----------
    normalized_data : _type_
        _description_
    type : _type_
        _description_
    background : _type_
        _description_
    saveFigure : bool, optional
        _description_, by default False
    """        


    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(8,2))
    plt.subplots_adjust(wspace=0.8)

    data=normalized_data.loc[background]

    data=data[~data.isin([np.nan, np.inf, -np.inf]).any(1)]      

    if type=="reads":
        
        
        cols=["Reads","Linear-norm-reads","reads_normalized_windows"]
        labels=["Raw reads","Norm. over chrom.","Norm. over 10kb windows"]

        fig.suptitle("Linear transformation of the reads",y=1.1,fontsize=18)
        
        for axes in ax:
            axes.set_ylim(0,3000)

        ax[0].set_xlim(0,np.max(data[cols[0]])/40)
        ax[1].set_xlim(0,np.max(data[cols[1]])/40)
        ax[2].set_xlim(0,np.max(data[cols[2]])/40)

        ax[0].set_ylabel("Frequency",fontsize=12)
        ax[1].set_ylabel(" ")   
        ax[2].set_ylabel(" ")

        sns.histplot(data[cols[0]],ax=ax[0],color="gray",binwidth=np.max(data[cols[0]])/1000)

        ax[0].set_xlabel(labels[0],fontsize=12)
        
        sns.histplot(data[cols[1]],ax=ax[1],color="blue",binwidth=np.max(data[cols[1]])/1000)
        
        ax[1].set_xlabel(labels[1],fontsize=12)

        sns.histplot(data[cols[2]],ax=ax[2],color="red", binwidth=np.max(data[cols[2]])/1000)
            
        ax[2].set_xlabel(labels[2],fontsize=12)
   

    elif type=="insertions":
        

        cols=["Insertions","Linear-norm-insertions","tr_normalized_windows"]
        labels=["Raw insertions","Norm. over chrom.","Norm. over 10kb windows"]
        
        fig.suptitle("Linear transformation of the insertions",y=1.1,fontsize=18)

        for axes in ax:
             axes.set_ylim(0,1000)

        ax[0].set_xlim(0,np.max(data[cols[0]])/15)
        ax[1].set_xlim(0,np.max(data[cols[1]])/15)
        ax[2].set_xlim(0,np.max(data[cols[2]])/15)

        ax[0].set_ylabel("Frequency",fontsize=12)
        ax[1].set_ylabel(" ")   
        ax[2].set_ylabel(" ")
        


        sns.histplot(data[cols[0]],ax=ax[0],color="gray",binwidth=np.max(data[cols[0]])/1000)

        ax[0].set_xlabel(labels[0],fontsize=12)
        
        sns.histplot(data[cols[1]],ax=ax[1],color="blue",binwidth=np.max(data[cols[1]])/1000)
        
        ax[1].set_xlabel(labels[1],fontsize=12)

        sns.histplot(data[cols[2]],ax=ax[2],color="red", binwidth=np.max(data[cols[2]])/1000)
            
        ax[2].set_xlabel(labels[2],fontsize=12)

    elif type=="transposon density":

        cols=["tr-density","Linear-norm-tr-density","length gene"]
        labels=["Insertion density (1/bp)","Norm. over chrom.(1/bp)","Gene length (bp)"]

        
        fig.suptitle("Transposon density normalization",y=1.1,fontsize=18)

        ax[0].set_xlim(0,np.max(data[cols[0]])/15)
        ax[1].set_xlim(0,np.max(data[cols[1]])/15)
        ax[2].set_xlim(0,np.max(data[cols[2]]/4))

        for axes in ax:
            axes.set_ylim(0,400)

        ax[0].set_ylabel("Frequency",fontsize=12)
        ax[1].set_ylabel(" ")   
        ax[2].set_ylabel(" ")
        


        sns.histplot(data[cols[0]],ax=ax[0],color="gray",binwidth=np.max(data[cols[0]])/1000)

        ax[0].set_xlabel(labels[0],fontsize=12)
        
        sns.histplot(data[cols[1]],ax=ax[1],color="blue",binwidth=np.max(data[cols[1]])/1000)
        
        ax[1].set_xlabel(labels[1],fontsize=12)

        sns.histplot(data[cols[2]],ax=ax[2],color="red", binwidth=np.max(data[cols[2]])/1000)
            
        ax[2].set_xlabel(labels[2],fontsize=12)

    elif type=="plot-genome-insertions": 
        
        coordinates_chrom=[]
        chrom=["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]

        y=[]
        y_linear=[]
        y_windows=[]

        cols=["Insertions","Linear-norm-insertions","tr_normalized_windows"]
        labels=["Raw insertions","Norm. over chrom.","Norm. over 10kb"]
        
        for i in chrom:
        
            coordinates_chrom.append((np.where(data["Chromosome"]==i)[0][0]))
            
            

            y.append(data[data["Chromosome"]==i][cols[0]])
            y_linear.append(data[data["Chromosome"]==i][cols[1]])
            y_windows.append(data[data["Chromosome"]==i][cols[2]])
        

        fig.suptitle("Genome wise data normalization",y=1.1,fontsize=18)

        sns.lineplot(data=np.concatenate(y),ax=ax[0],color="black")
        sns.lineplot(data=np.concatenate(y_linear),ax=ax[1],color="blue")
        sns.lineplot(data=np.concatenate(y_windows),ax=ax[2],color="red")

        ax[0].set_ylabel(labels[0],fontsize=12)
        ax[1].set_ylabel(labels[1],fontsize=12)
        ax[2].set_ylabel(labels[2],fontsize=12)

        for axes in ax:
            axes.set_xlabel("Gene positions",fontsize=12)
            

    elif type=="plot-genome-reads": 
    
        coordinates_chrom=[]
        chrom=["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]

        y=[]
        y_linear=[]
        y_windows=[]

        cols=["Reads","Linear-norm-reads","reads_normalized_windows"]
        labels=["Raw reads","Norm. over chrom.","Norm. over 10kb"]

        for i in chrom:
        
            coordinates_chrom.append((np.where(data["Chromosome"]==i)[0][0]))

            y.append(data[data["Chromosome"]==i][cols[0]])
            y_linear.append(data[data["Chromosome"]==i][cols[1]])
            y_windows.append(data[data["Chromosome"]==i][cols[2]])

        fig.suptitle("Genome wise data normalization",y=1.1,fontsize=18)

        sns.lineplot(data=np.concatenate(y),ax=ax[0],color="black")
        sns.lineplot(data=np.concatenate(y_linear),ax=ax[1],color="blue")
        sns.lineplot(data=np.concatenate(y_windows),ax=ax[2],color="red")

        ax[0].set_ylabel(labels[0],fontsize=12)
        ax[1].set_ylabel(labels[1],fontsize=12)
        ax[2].set_ylabel(labels[2],fontsize=12)

        for axes in ax:
            axes.set_xlabel("Gene positions",fontsize=12)


    elif type=="plot-genome-density": 
    
        coordinates_chrom=[]
        chrom=["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]

        y=[]
        y_linear=[]
        y_windows=[]

        cols=["tr-density","Linear-norm-tr-density","length gene"]
        labels=["Insertion density (1/bp)","Over chrom.(1/bp)","Gene length (bp)"]

        for i in chrom:
        
            coordinates_chrom.append((np.where(data["Chromosome"]==i)[0][0]))

            y.append(data[data["Chromosome"]==i][cols[0]])
            y_linear.append(data[data["Chromosome"]==i][cols[1]])
            y_windows.append(data[data["Chromosome"]==i][cols[2]])

        fig.suptitle("Genome wise data normalization",y=1.2,fontsize=18)

        sns.lineplot(data=np.concatenate(y),ax=ax[0],color="black")
        sns.lineplot(data=np.concatenate(y_linear),ax=ax[1],color="blue")
        sns.lineplot(data=np.concatenate(y_windows),ax=ax[2],color="red")

        ax[0].set_title(labels[0],fontsize=12)
        ax[1].set_title(labels[1],fontsize=12)
        ax[2].set_title(labels[2],fontsize=12)

        for axes in ax:
            axes.set_xlabel("Gene positions",fontsize=12)
            

            
    for axes in ax:
        axes.tick_params(axis='x', labelsize=16)
        axes.tick_params(axis='y', labelsize=16)

    

    

    if saveFigure:
        fig.savefig("../figures/linear_transformation_"+background+"_"+type+".png",
        bbox_inches="tight",dpi=400)

    
# -

linear_transformations_plots(data_norm_pd,type="plot-genome-insertions",background="wt_merged")

types=["insertions","reads","transposon density","plot-genome-reads",
"plot-genome-insertions","plot-genome-density"]
for i in types:
    linear_transformations_plots(data_norm_pd,i,"wt_merged",saveFigure=True)


# +
def linear_transformations_plots_per_chrom(normalized_data,background,type,chrom="I",
saveFigure=False):
    """_summary_

    Parameters
    ----------
    normalized_data : _type_
        _description_
    background : _type_
        _description_
    type : _type_
        _description_
    chrom : str, optional
        _description_, by default "I"
    saveFigure : bool, optional
        _description_, by default False
    """


    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,2))
    
    data=normalized_data.loc[background]
    data=data[data["Chromosome"]==chrom]
    data=data[~data.isin([np.nan, np.inf, -np.inf]).any(1)] 
    
    
    if type=="plot-genome-density": 
    
        cols=["tr-density","Linear-norm-tr-density"]
        labels=["Insertion density (1/bp)","Over chrom.(1/bp)"]

        

        fig.suptitle("Normalization over chrom:"+chrom,y=1.2,fontsize=18)

        sns.lineplot(data=data.loc[:,cols[0]],ax=ax[0],color="black")
        sns.lineplot(data=data.loc[:,cols[1]],ax=ax[1],color="red")

        ax[0].set_title(labels[0],fontsize=12)
        ax[1].set_title(labels[1],fontsize=12)

    elif type=="plot-genome-insertions":
            
            cols=["Insertions","tr_normalized_windows"]
            labels=["Raw insertions","Over 10kb on chrom."]
    
            fig.suptitle("Normalization over 10kb on chrom:"+chrom,y=1.2,fontsize=18)
    
            sns.lineplot(data=data.loc[:,cols[0]],ax=ax[0],color="black")
            sns.lineplot(data=data.loc[:,cols[1]],ax=ax[1],color="red")
    
            ax[0].set_title(labels[0],fontsize=12)
            ax[1].set_title(labels[1],fontsize=12)

    elif type=="plot-genome-reads":
        
        cols=["Reads","reads_normalized_windows"]
        labels=["Raw reads","Over 10kb on chrom."]

        fig.suptitle("Normalization over 10kb on chrom:"+chrom,y=1.2,fontsize=18)

        sns.lineplot(data=data.loc[:,cols[0]],ax=ax[0],color="black")
        sns.lineplot(data=data.loc[:,cols[1]],ax=ax[1],color="red")

        ax[0].set_title(labels[0],fontsize=12)
        ax[1].set_title(labels[1],fontsize=12)
        

    for axes in ax:
        axes.tick_params(axis='x', labelsize=16)
        axes.tick_params(axis='y', labelsize=16)
        axes.set_xlabel("Gene positions",fontsize=12)
        axes.set_ylabel(" ",fontsize=12)
    
    
            

    if saveFigure:
        fig.savefig("../figures/linear_transformation_"+background+"_"+type+"_"+chrom+".png",
        bbox_inches="tight",dpi=400)
            
    
        
# -

linear_transformations_plots_per_chrom(data_norm_pd,"wt_merged","plot-genome-insertions",chrom="XV",saveFigure=True)

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
    sns.regplot(x=x,y=y,ax=ax[i//3,i%3],color="black",seed=1,ci=95,order=1)
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


sns.lineplot(data=x_reads,ax=ax[0,0],color="black",ci="sd",err_style="band");
ax[0,0].set_xticks(np.arange(0,len(chrom),1))
ax[0,0].set_xticklabels(chrom,rotation=45,fontsize=12);
ax[0,0].set_title("Reads per chrom",fontsize=16)

sns.lineplot(data=x_insertions,ax=ax[0,1],color="black",ci="sd",err_style="band");
ax[0,1].set_xticks(np.arange(0,len(chrom),1))
ax[0,1].set_xticklabels(chrom,rotation=45,fontsize=12);
ax[0,1].set_title("Insertions per chrom",fontsize=16)


sns.lineplot(data=data["reads_over_windows"],ax=ax[1,0],color="black",ci="sd",err_style="band");

ax[1,0].set_title("Reads over 10kb windows",fontsize=16)
ax[1,0].set_ylabel(" ")
# ax[1,0].tick_params(axis='x', labelsize=12)
ax[1,0].set_xlabel("Genomic locations",fontsize=14)

sns.lineplot(data=data["insertions_over_windows"],ax=ax[1,1],color="black",ci="sd",err_style="band");

ax[1,1].set_title("Insertions over 10kb windows",fontsize=16)
ax[1,1].set_ylabel(" ")
ax[1,1].set_xlabel("Genomic locations",fontsize=14)

for i in range(ax.shape[0]):
    for j in range(ax.shape[1]):
        ax[i,j].tick_params(axis='x', labelsize=14)
        ax[i,j].tick_params(axis='y', labelsize=14)

fig.savefig("../figures/figures_thesis_chapter_2/reads_insertions_per_chrom.png",bbox_inches="tight",dpi=400)
# -

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

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

# ## Transposon insertion data
# ### Describe our first dataset on WT. 
# - Show the differences in read count and transposon counts in technical replicates( they are replicates splited before the PCR). Cite Enzo thesis to say that PCR amplification is one of the major noise sources in the read counts of satay libraries, to explain the differences across replicates.
# - Show the known biases of this type of transposon system to genes located near to centromeres. 
# - Plot the cumulative plot of transposons against the distance of every gene to the centromere of its chromosome. 
# - Show the tendency of transposons to insert in non coding regions and that this does not change  across different genetic backgrounds (wt,dnrp1,dbem3,dbem1dbem3) 
#

# +
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys
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

keys

data=list_data_pd.loc["wt_merged"]
data.head()

# ## Using read counts for fitness 
# - Simplified fitness model (simple malthusian model where the differences in mutant abundance is proportional to its growth rate, and there is no limiting factor)
#
# For the fitness calculation:
#
# - Take only the 80% central part of the gene (show the overall enrichment towards the first and last part in average across all genes.)
# - Annotate which genes "we can not analyse" because they dont have enough data. So, the challenge here is to define what is a region where we dont have sufficient data to say something about it. 
#     - use the density of the flanking regions as the reference for the gene transposon density.
#
# - For the "analysable genes" compute the fitness using the maltusian model of the:
#     - whole gene by averaging the read counts over the number of transposons
#     - subdivide the gene into known domains and compute the fitness of individual domains
#     - assign a new corrected fitness to the every gene by taking the fitness of the domain with strongest effect with respect to the average , so the max|f_mean-f_domain_i|
#  
#  Then :
#  - Plot the distributions of every fitness calculation referenced to the median fitness value of the population.
#  - Plot the average fitness vs the corrected fitness 
#
#
#

# +
reads_locations_center=[]
reads_locations=[]
insertion_locations=[]
insertion_locations_center=[]

gene_coordinates=[]
data=list_data_pd.loc["wt_merged"]
gene_len_data=[]

read_insert_pair=defaultdict(dict)
for j in data.index:
        coi=from_excel_to_list(data.loc[j]["Reads per insertion location"])
        coi_insertions=from_excel_to_list(data.loc[j]["Insertion locations"])
        gene_coordinates.append([data.loc[j,"Start location"],data.loc[j,"End location"]])
        reads_locations.append(coi)
        insertion_locations.append(coi_insertions)
        read_insert_pair[j]=(coi,coi_insertions)

        gene_len=data.loc[j,"End location"]-data.loc[j,"Start location"]
        gene_len_data.append(gene_len)
        if type(coi)!=int: # for genes that have more than one read
            start=gene_len*0.1
            end=gene_len*0.9
            data_center=coi[int(start):int(end)]
            data_center_insertions=coi_insertions[int(start):int(end)]
            reads_locations_center.append(data_center)
            insertion_locations_center.append(data_center_insertions)
        else:
            data_center=coi
            reads_locations_center.append(coi)
            insertion_locations_center.append(coi_insertions)

# +
# compute the reads per insertions for each gene along the gene length

gene_parts=np.linspace(0,1,21)
r=np.zeros(shape=(len(insertion_locations),len(gene_parts))) # reads_per_insertion_parts array , every gene in the rows 

for i in np.arange(0,len(insertion_locations)):
    if (insertion_locations[i])!=0:
        g=np.array(insertion_locations[i])
        f=np.array(gene_coordinates[i][0]+gene_coordinates[i][1]*gene_parts)
        binedges = g.searchsorted(f)

        rngs = [list(range(binedges[k], binedges[k+1])) for k in range(len(binedges)-1)]

        for k in np.arange(0,len(rngs)):
            readsperinsert=[]
            for j in np.arange(0,len(rngs[k])):
                readsperinsert.append(reads_locations[i][rngs[k][j]])
            r[i,k]=np.sum(readsperinsert)/len(readsperinsert)

    



# +
# replace nan values with zeros

r[np.isnan(r)]=0


# Sum all values of r along the rows

r_sum=np.sum(r,axis=0)

# Plot r values along gene parts

plt.figure(figsize=(5,5))

for i in np.arange(1,len(insertion_locations)):
    plt.plot(gene_parts,r[i,:],color="black",alpha=0.3,marker="o",markersize=5)

plt.plot(gene_parts,r[0,:],color="black",alpha=0.3,marker="o",markersize=5,
    label="Data per gene")
plt.yscale("log")
plt.xlabel("Gene length")
plt.ylabel("Reads per insertion per gene")

plt.bar(gene_parts,r_sum,color="black",alpha=0.1,width=0.05,label="Sum of all genes")
plt.xticks(np.arange(0,1.1,0.1),labels=["0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"]);

plt.legend()
plt.tight_layout()

plt.savefig("../figures/fig_reads_per_insertion_all_genes_along_gene_length.png",dpi=300,format="png")
## from the figure it follows that I should remove only the first 10% of the reads to compute the fitness of the whole gene


# +

keys= ['bem1-aid_a','dbem1dbem3_b','wt_merged','dbem1dbem3_a', 
'dnrp1_merged','bem1-aid_b','dbem3_merged']



# +
genes_out_by_neighborhood_pd=pd.read_excel("../postprocessed-data/genes_out_by_neighborhood.xlsx",index_col="Unnamed: 0")
genes_out_float=defaultdict(dict)
for key in keys:
   genes_out_float["sum upstream insertions"][key] =from_excel_to_list(genes_out_by_neighborhood_pd.loc[key,"sum upstream insertions"])
   genes_out_float["sum downstream insertions"][key] =from_excel_to_list(genes_out_by_neighborhood_pd.loc[key,"sum downstream insertions"])
   x=genes_out_by_neighborhood_pd.loc[key,"discarded_genes_neighborhood"]
   x=x.replace('[', '')
   x=x.replace(']', '')
   x=x.replace('list(', '')
   x=x.replace(')', '')
   x=x.replace("'", "")
   x=x.replace(" ", "")
   x=x.split(',')
   genes_out_float["discarded_genes_neighborhood"][key]=x
   
genes_out_float_pd=pd.DataFrame.from_dict(genes_out_float)

genes_out_float_pd.loc[:,"threshold coverage"]=genes_out_by_neighborhood_pd.loc[:,"threshold coverage"]

# +
no_data_genes=genes_out_float_pd.loc["wt_merged"]["discarded_genes_neighborhood"]

no_data_genes=np.array(no_data_genes)
index_gene=np.where(no_data_genes=='MF(ALPHA1')
no_data_genes[index_gene[0][0]]="MF(ALPHA)1"

index_gene=np.where(no_data_genes=='MF(ALPHA2')
no_data_genes[index_gene[0]]="MF(ALPHA)2"

index_gene=np.where(no_data_genes=='YTH1\nAAD15')
no_data_genes[index_gene[0]]="YTH1"
no_data_genes[index_gene[0]+1]="AAD15"
# -

fitness_data=pd.DataFrame(columns=["fitness","classification"],index=list_data_pd.loc[:,"Gene name"])
for i in no_data_genes:
    fitness_data.loc[i,"classification"]="insufficient_data"

fitness_data.head()

# +
no_data_genes



# +
## why is it insufficient data
up=genes_out_float_pd.loc["wt_merged"]["sum upstream insertions"] # sum of upstream insertions 10kb upstream of the gene
down=genes_out_float_pd.loc["wt_merged"]["sum downstream insertions"] # sum of downstream insertions 10kb downstream of the gene
threshold=genes_out_float_pd.loc["wt_merged"]["threshold coverage"] # threshold coverage based on the average density of the library for the gene to be considered as having sufficient data

goi=["AOS1"]
index_goi=np.where(no_data_genes==goi)
data_copy=list_data_pd.copy()
data_copy.index=data_copy.loc[:,"Gene name"]
data_goi=data_copy.loc[goi[0],"Insertions"]

plt.figure(figsize=(5,5))

plt.bar(["10kb up","gene","10kb down"],[up[index_goi[0][0]],data_goi.mean(),
down[index_goi[0][0]]],color=["green","black","purple"],alpha=0.5)
plt.hlines(threshold,0,2,linestyles="dashed",color="gray",label="threshold")
plt.text(0,up[index_goi[0][0]]+0.1,str(np.round(up[index_goi[0][0]],2)),ha="center")
plt.text(1,np.round(data_goi.mean(),2)+0.1,str(np.round(data_goi.mean(),2)),ha="center")
plt.text(2,down[index_goi[0][0]]+0.1,str(np.round(down[index_goi[0][0]],2)),ha="center")
plt.text(1,threshold+0.1,str(np.round(threshold,2)),ha="center",color="gray",fontsize=14)

#plt.ylim(0,np.amax([up[index_goi[0][0]],data_goi.mean(),down[index_goi[0][0]]])+1)
plt.yscale("log")
plt.ylabel("Number of insertions")
plt.title(goi[0])
plt.legend()

# +
## why is it insufficient data
up=genes_out_float_pd.loc["wt_merged"]["sum upstream insertions"] # sum of upstream insertions 10kb upstream of the gene
down=genes_out_float_pd.loc["wt_merged"]["sum downstream insertions"] # sum of downstream insertions 10kb downstream of the gene
threshold=genes_out_float_pd.loc["wt_merged"]["threshold coverage"] # threshold coverage based on the average density of the library for the gene to be considered as having sufficient data
data_copy=list_data_pd.copy()
data_copy.index=data_copy.loc[:,"Gene name"]

plt.subplots(4,4,figsize=(15,15))
plt.subplots_adjust(hspace=0.5,wspace=0.5)
j=1
for i in no_data_genes[0:16]:
    goi=i
    index_goi=np.where(no_data_genes==goi)
    data_goi=data_copy.loc[goi,"Insertions"]

    plt.subplot(4,4,j)

    plt.bar(["10kb up","gene","10kb down"],[up[index_goi[0][0]],data_goi.mean(),
    down[index_goi[0][0]]],color=["green","black","purple"],alpha=0.5)
    plt.hlines(threshold,0,2,linestyles="dashed",color="gray",label="threshold")
    plt.text(0,up[index_goi[0][0]]+0.1,str(np.round(up[index_goi[0][0]],2)),ha="center")
    plt.text(1,np.round(data_goi.mean(),2)+0.1,str(np.round(data_goi.mean(),2)),ha="center")
    plt.text(2,down[index_goi[0][0]]+0.1,str(np.round(down[index_goi[0][0]],2)),ha="center")

    plt.text(1,threshold+0.1,str(np.round(threshold,2)),ha="center",color="gray",fontsize=14)

    plt.yscale("log")
    plt.ylabel("Number of insertions")
    plt.title(goi)
    plt.legend()
    j=j+1

#plt.tight_layout()

plt.savefig("../figures/fig_examples_insufficient_data_genes.png",dpi=300)
# -



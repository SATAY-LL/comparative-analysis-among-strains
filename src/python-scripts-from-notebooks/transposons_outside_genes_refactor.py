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

from functions_transposons_outside_genes import get_genes_names
from functions_transposons_outside_genes import defining_threshold_given_tr_density
from functions_transposons_outside_genes import get_discarded_genes_by_duplication
from functions_transposons_outside_genes import get_amenable_genes_coverage_neighborhood
from functions_transposons_outside_genes import local_discrimination_genes_by_neighbors_coverage
# -

data_all=pd.read_excel('../postprocessed-data/postprocessed_data_all_backgrounds.xlsx',engine="openpyxl")

# +

keys=['dbem3_b',
 'dnrp1_b',
 'bem1-aid_a',
 'dnrp1_a',
 'dbem1dbem3_b',
 'wt_merged',
 'dbem1dbem3_a',
 'bem1-aid-dbem3_a',
 'bem1-aid-dbem3_b',
 'wt_b',
 'wt_a',
 'dnrp1_merged',
 'bem1-aid_b',
 'dbem3_merged',
 'dbem3_a']

# +

positions_float_pd=pd.read_csv("../postprocessed-data/genetic_positions_float_pd_all_backgrounds.csv",converters={'Positions_float': literal_eval,'Ninsertions': literal_eval})
positions_float_pd.rename(columns={'Unnamed: 0':'Gene name', "Unnamed: 1": "background"},inplace=True)
# -

data_all_modified,genes_names=get_genes_names(data_raw=data_all,keys=keys)

# +

# ## getting genes that are duplicated in the genome 


# discarded_genes_by_duplication=get_discarded_genes_by_duplication(positions_float_pd,genes_names=genes_names)

# np.savetxt("../postprocessed-data/discarded_genes_by_duplication.txt",discarded_genes_by_duplication,fmt="%s")

# -

discarded_genes_by_duplication=np.loadtxt("../postprocessed-data/discarded_genes_by_duplication.txt",dtype=str)

discarded_genes_by_duplication

# +
## Getting amenable genes for analysis of local neighborhood

# Independent variables 

windows_size=10000
targets,genes_not_discarded_by_location=get_amenable_genes_coverage_neighborhood(positions_float_pd,genes_names=genes_names,discarded_genes_by_duplication=discarded_genes_by_duplication,windows_size=windows_size)

# -

genes_not_discarded_by_location

# +

## example playground one gene
background="wt_merged"
## defining the Threshold for each background
threshold,density=defining_threshold_given_tr_density(data_all_modified,windows_size=windows_size,background=background)
## Discriminating genes by local neighborhood coverage
n=0
k=10
a,b,c,d=local_discrimination_genes_by_neighbors_coverage(positions_float_pd,background=background,gene_of_interest="BEM1",windows_size=windows_size,threshold=k*threshold)

# -

a,b,c,d,k*threshold,density

# +
## Big loop to get the local neighborhood coverage for each background 
background="wt_merged"
k=1 # Amplified factor for the threshold
windows_size=3000 # bp
genes_out=[]
tmp_c=[]
tmp_d=[]
genes_out_by_neighborhood=defaultdict(dict)

threshold,density=defining_threshold_given_tr_density(data_all_modified,windows_size=windows_size,background=background)

targets,genes_not_discarded_by_location=get_amenable_genes_coverage_neighborhood(positions_float_pd,genes_names=genes_names,discarded_genes_by_duplication=discarded_genes_by_duplication,windows_size=windows_size)

for amenable_genes in targets:
    
    a,b,c,d=local_discrimination_genes_by_neighbors_coverage(positions_float_pd,background=background,gene_of_interest=amenable_genes,windows_size=windows_size,threshold=k*threshold)
    if any((a!=[],b!=[])):
        
        tmp_c.append(c) # append sum upstream insertions
        tmp_d.append(d) # append sum downstream insertions
    if a!=[]:
        genes_out.append(a)
    elif b!=[]:
        genes_out.append(b)


genes_out_by_neighborhood["discarded_genes_neighborhood"][background]=np.unique(genes_out)
genes_out_by_neighborhood["sum upstream insertions"][background]=tmp_c
genes_out_by_neighborhood["sum downstream insertions"][background]=tmp_d
genes_out_by_neighborhood["threshold coverage"][background]=k*threshold    
# -

keys= ['bem1-aid_a','dbem1dbem3_b','wt_merged','dbem1dbem3_a', 
'dnrp1_merged','bem1-aid_b','dbem3_merged']

# +
#keys=["wt_a","wt_b"]

# +
### Run for a night #######DONT RUN THIS AGAIN IT TAKES 5 HOURS!!
k=1 # Amplified factor for the threshold
windows_size=10000 #bp
genes_out_by_neighborhood=defaultdict(dict)

for background in keys:
    tmp_a=[]
    tmp_b=[]
    tmp_c=[]
    tmp_d=[]

    threshold,density=defining_threshold_given_tr_density(data_all_modified,windows_size=windows_size,background=background)
    targets,genes_not_discarded_by_location=get_amenable_genes_coverage_neighborhood(positions_float_pd,genes_names=genes_names,discarded_genes_by_duplication=discarded_genes_by_duplication,windows_size=windows_size)

    for amenable_genes in targets:
       
        a,b,c,d=local_discrimination_genes_by_neighbors_coverage(positions_float_pd,background=background,gene_of_interest=amenable_genes,windows_size=windows_size,threshold=k*threshold)
        if a!=[]:
            tmp_a.append(a)

            
        if b!=[]:
            tmp_b.append(b)

        if any((a!=[],b!=[])):

            tmp_c.append(c)
            tmp_d.append(d)

    genes_out_by_neighborhood["discarded_genes_neighborhood"][background]=np.unique([tmp_a,tmp_b])
    genes_out_by_neighborhood["sum upstream insertions"][background]=tmp_c
    genes_out_by_neighborhood["sum downstream insertions"][background]=tmp_d
    genes_out_by_neighborhood["threshold coverage"][background]=k*threshold


# +
genes_out_by_neighborhood_pd=pd.DataFrame.from_dict(genes_out_by_neighborhood)

genes_out_by_neighborhood_pd

# +
#genes_out_by_neighborhood_pd.to_excel("../postprocessed-data/genes_out_by_neighborhood.xlsx",index=True)

# +

genes_out_by_neighborhood_pd=pd.read_excel("../postprocessed-data/genes_out_by_neighborhood.xlsx",index_col="Unnamed: 0")
# -

from from_excel_to_list import from_excel_to_list
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
   x=x.split(',')
   genes_out_float["discarded_genes_neighborhood"][key]=x
   

# +
genes_out_float_pd=pd.DataFrame.from_dict(genes_out_float)

genes_out_float_pd.loc[:,"threshold coverage"]=genes_out_by_neighborhood_pd.loc[:,"threshold coverage"]
# -

genes_out_float_pd

# +

for i in keys:
    print("There are",len(genes_out_float_pd.loc[i,"discarded_genes_neighborhood"]),
        "genes discarded in" ,i, "for essentiality scores due to its low coverage in its neighborhood.")
# -

for j in np.arange(0,len(keys)):
    for i in np.arange(0,len(keys)):
        if i!=j:
            tmp_1=np.setdiff1d(genes_out_float_pd.loc[keys[j],"discarded_genes_neighborhood"],genes_out_float_pd.loc[keys[i],"discarded_genes_neighborhood"])
            tmp_2=np.setdiff1d(genes_out_float_pd.loc[keys[i],"discarded_genes_neighborhood"],genes_out_float_pd.loc[keys[j],"discarded_genes_neighborhood"])

            print("these genes:",tmp_1,"are genes that are discarded in",keys[j],"that are not in",keys[i])
            print("these genes:",tmp_2,"are genes that are discarded in",keys[i],"that are not in",keys[j])

fig, ax = plt.subplots(nrows=len(keys), ncols=1, figsize=(5, 15))
for i in np.arange(0,len(keys)):
    ax[i].hist(genes_out_float_pd.loc[keys[i],"sum downstream insertions"], bins=400, alpha=0.5, label='downstream genes_'+keys[i],color='purple')
    # plt.legend(loc='upper right')
    ax[i].set_xlim([0,1000])
    ax[i].set_xlabel('Sum of downstream insertions per gene discarded')

    ax[i].vlines(x=genes_out_float_pd.loc[keys[i],"threshold coverage"] ,ymin=0, ymax=140, linestyles='dashed', colors='r', label='Threshold_'+ keys[i])

    ax[i].legend(loc='upper right')
#fig.savefig("../figures/fig_histogram_downstream_insertions_per_gene_discarded_by_neighborhood_all_keys.png",dpi=300)

fig, ax = plt.subplots(nrows=len(keys), ncols=1, figsize=(5, 15))
for i in np.arange(0,len(keys)):
    ax[i].hist(genes_out_float_pd.loc[keys[i],"sum upstream insertions"], bins=100, alpha=0.5, label='upstream genes_'+keys[i],color='green')
    # plt.legend(loc='upper right')
    ax[i].set_xlim([0,600000])
    ax[i].set_ylim([0,20])
    ax[i].set_xlabel('Sum of upstream insertions per gene discarded')

    ax[i].vlines(x=genes_out_float_pd.loc[keys[i],"threshold coverage"] ,ymin=0, ymax=140, linestyles='dashed', colors='r', label='Threshold_'+ keys[i])

    ax[i].legend(loc='upper right')
#fig.savefig("../figures/fig_histogram_upstream_insertions_per_gene_discarded_by_neighborhood_all_keys.png",dpi=300)



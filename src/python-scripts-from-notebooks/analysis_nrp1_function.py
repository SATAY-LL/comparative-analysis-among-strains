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
#     display_name: Python 3.9.7 ('transposonmapper')
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

from scipy import stats

plt.rc('font', family='serif',size=14)
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

<<<<<<< HEAD
from functions_interaction_computations import filter_fitness
from functions_interaction_computations import digenic_GI
from functions_interaction_computations import classify_GI

=======
>>>>>>> 1c5527a6c3a1e3cda5126ec83fb46e9a8a798f6f
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
    tmp=pd.read_excel(i,engine='openpyxl',index_col="Unnamed: 0")
    ## remove ADE2 genes
    tmp=tmp[tmp.loc[:,"Gene name"]!="ADE2"]
    tmp.index=np.arange(0,len(tmp))
    list_data.append(tmp)

keys=[]
for i in np.arange(0,len(pergene_files)):
    keys.append(pergene_files[i].split("/")[-1].split("_")[0]+"_"+pergene_files[i].split("/")[-1].split("_")[1])

list_data_pd=pd.concat(list_data,axis=0,keys=keys)

# -

keys,len(keys)

# +
import pickle
with open("../postprocessed-data/fitness_models_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

fitness_all_pd=pd.concat(b,axis=0,keys=keys)

with open("../postprocessed-data/over-regulated-genes-dnrp1", "rb") as fp:   # Unpickling
    over_nrp1 = pickle.load(fp)

with open("../postprocessed-data/down-regulated-genes-dnrp1", "rb") as fp:   # Unpickling
    under_nrp1 = pickle.load(fp)

# -

standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str)

# +

polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")
polarity_genes.fillna(0,inplace=True)
<<<<<<< HEAD
# -

data_fitness=filter_fitness(fitness_all_pd,backgrounds=keys,goi=["BEM1","BEM3","NRP1"],discard=["Not enough flanking regions"],set2zero=["Not enough reads",
    "Not enough insertions"],cols=["fitness_gene","fitness_domains_corrected"])

data_fitness.loc["wt_a"].loc["BEM1"]
=======

# +
data=[]
for backg in keys:
    f=fitness_all_pd.loc[backg]
    f=f[f.loc[:,"fitness_gene"]!="Not enough flanking regions"]
    for i in f.index:
        if f.loc[i,"fitness_gene"]=="Not enough reads":
            f.loc[i,"fitness_gene"]=0
        elif f.loc[i,"fitness_gene"]=="Not enough insertions":
            f.loc[i,"fitness_gene"]=0
        elif f.loc[i,"fitness_domains_corrected"]=="Not enough insertions":
            f.loc[i,"fitness_domains_corrected"]=0

    # f[f.loc[:,"fitness_gene"]=="Not enough reads"]=0
    # f[f.loc[:,"fitness_gene"]=="Not enough insertions"]=0
    # f[f.loc[:,"fitness_domains_corrected"]=="Not enough insertions"]=0
    
    # f=f[f.loc[:,"fitness_gene"]!="Not enough reads"]
    # f=f[f.loc[:,"fitness_gene"]!="Not enough insertions"]
    # f=f[f.loc[:,"fitness_domains_corrected"]!="Not enough insertions"]
    data.append(f)

data_fitness=pd.concat(data,axis=0,keys=keys)
>>>>>>> 1c5527a6c3a1e3cda5126ec83fb46e9a8a798f6f

# +
x=data_fitness.loc["wt_merged"].loc[:,"fitness_gene"]
y=data_fitness.loc["dnrp1_merged"].loc[:,"fitness_gene"]
z=data_fitness.loc["bem1-aid_merged"].loc[:,"fitness_gene"]

<<<<<<< HEAD
xa=data_fitness.loc["wt_a"].loc[:,"fitness_gene"]
xb=data_fitness.loc["wt_b"].loc[:,"fitness_gene"]

# +


# select the common indexe.setdefault()
common_index = set(x.index).intersection(set(y.index))
common_index=set(common_index).intersection(set(z.index))
common_index_ab=set(xa.index).intersection(set(xb.index))


plt.hist(xa.loc[common_index_ab]-xb.loc[common_index_ab],color="red",alpha=0.4,label="wt_a-wt_b");
plt.hist(y.loc[common_index]-x.loc[common_index],alpha=0.7,label="dnrp1-wt",color="green");
plt.hist(z.loc[common_index]-x.loc[common_index],alpha=0.3,label="bem1-aid-wt",color="blue");
plt.legend()

# +
bem1_aid=list_data_pd.loc["bem1-aid_merged"]
wt_merged=list_data_pd.loc["wt_merged"]
bem1_aid.index=bem1_aid.loc[:,"Gene name"]

bem1_aid.Reads.sum(),wt_merged.Reads.sum(),bem1_aid.Insertions.sum(),wt_merged.Insertions.sum()

# +
## Plotting that BEM1 is highly favored in bem1-aid background


bem1_aid=list_data_pd.loc["bem1-aid_merged"]
bem1_aid.index=bem1_aid.loc[:,"Gene name"]
x=bem1_aid.Insertions.sort_values(ascending=False).head(15)
y=bem1_aid.Reads.sort_values(ascending=False).head(15)
ax,fig=plt.subplots(1,2,figsize=(10,5))
plt.subplots_adjust(wspace=0.3)
plt.subplot(1,2,1)
plt.bar(x.index,x,color="gray")
plt.ylabel("Transposon insertions")
plt.xticks(rotation=90);
plt.subplot(1,2,2)
plt.bar(y.index,y,color="gray")
plt.ylabel("Reads")
plt.xticks(rotation=90);
plt.tight_layout()
plt.savefig("../figures/bem1_aid_bem1_overrepresented.png",dpi=300)

# +
## plot the represssion to bem1d from L Laan eLife paper from evolutionary experiments   

growth_rate_steps=[0.0082,0.001,0.006,0.007]

labels=["WT","dbem1","dbem1dbem3","dbem1dbem3dnrp1"]

plt.figure(figsize=(5,5))

#plt.plot(labels[1:],growth_rate_steps[1:],"o",color="black",markersize=10)
plt.bar(x=labels,height=growth_rate_steps,color="gray",alpha=0.6)
#plt.plot(growth_rate_steps[1:],growth_rate_steps[1:],"--",color="black",markersize=10,alpha=0.3)
#plt.hlines(growth_rate_steps[0],0,2,linestyles="dashed",label="Wild type ",color="purple")
plt.xticks(rotation=90);
plt.ylabel("Growth rate (min$^{-1}$)")
plt.grid(linewidth=0.2)

plt.ylim(0,0.012)
plt.tight_layout()
plt.savefig("../figures/fitness_laan_evolution.png",dpi=300)

# +
##Analysis evolutionary trajectory of dbem1 cells


## plot the represssion to bem1d from L Laan eLife paper from evolutionary experiments   

growth_rate_steps=[0.0082,0.001,0.006,0.007]

labels=["WT","dbem1","dbem1dbem3","dbem1dbem3dnrp1"]

plt.figure(figsize=(5,5))

#plt.plot(labels[1:],growth_rate_steps[1:],"o",color="black",markersize=10)
#plt.bar(x=labels[1:],height=growth_rate_steps[1:],color="gray",alpha=0.6)
plt.plot(labels,growth_rate_steps,"o",color="black",markersize=15,alpha=0.3)
#plt.hlines(growth_rate_steps[0],0,2,linestyles="dashed",label="Wild type ",color="purple")
plt.xticks(rotation=90);
plt.ylabel("Growth rate (min$^{-1}$)")
plt.grid(linewidth=0.2)

plt.ylim(0,0.012)
plt.tight_layout()
plt.savefig("../figures/evolutionary_trajectory_growth_rate_dots.png",dpi=300)

# +
backgrounds=["wt_a","wt_b","dnrp1_1","dnrp1_2","bem1-aid_a","bem1-aid_b"]

data=list_data_pd.loc[backgrounds]

## Plot the fold change of insertions from the two backgrounds in dnrp1 , in WT and dbem1

data_wta=data.loc["wt_a"]
data_wtb=data.loc["wt_b"]

data_bem1a=data.loc["bem1-aid_a"]
data_bem1b=data.loc["bem1-aid_b"]

# +
## Represent the fold change of insertions and reads from the two backgrounds in dnrp1 , in WT and dbem1 



data_insertions_nrp1a=data_wta[data_wta.loc[:,"Gene name"]=="NRP1"]["Insertions"].values[0]/data_wta.loc[:,"Insertions"].sum()
data_insertions_nrp1b=data_wtb[data_wtb.loc[:,"Gene name"]=="NRP1"]["Insertions"].values[0]/data_wtb.loc[:,"Insertions"].sum()

data_insertions_bem1a=data_wta[data_wta.loc[:,"Gene name"]=="BEM1"]["Insertions"].values[0]/data_wta.loc[:,"Insertions"].sum()
data_insertions_bem1b=data_wtb[data_wtb.loc[:,"Gene name"]=="BEM1"]["Insertions"].values[0]/data_wtb.loc[:,"Insertions"].sum()

data_insertions_nrp1_bem1a=data_bem1a[data_bem1a.loc[:,"Gene name"]=="NRP1"]["Insertions"].values[0]/data_bem1a.loc[:,"Insertions"].sum()
data_insertions_nrp1_bem1b=data_bem1b[data_bem1b.loc[:,"Gene name"]=="NRP1"]["Insertions"].values[0]/data_bem1b.loc[:,"Insertions"].sum()

values2plot_mean=[np.mean([data_insertions_nrp1a,data_insertions_nrp1b]),np.mean([data_insertions_bem1a,data_insertions_bem1b]),
np.mean([data_insertions_nrp1_bem1a,data_insertions_nrp1_bem1b])]

=======
# +
x=data_fitness.loc["wt_merged"].loc[:,"fitness_gene"]
y=data_fitness.loc["dnrp1_merged"].loc[:,"fitness_gene"]
z=data_fitness.loc["bem1-aid_merged"].loc[:,"fitness_gene"]

xa=data_fitness.loc["wt_a"].loc[:,"fitness_gene"]
xb=data_fitness.loc["wt_b"].loc[:,"fitness_gene"]

# +


# select the common indexe.setdefault()
common_index = set(x.index).intersection(set(y.index))
common_index=set(common_index).intersection(set(z.index))
common_index_ab=set(xa.index).intersection(set(xb.index))


plt.hist(xa.loc[common_index_ab]-xb.loc[common_index_ab],color="red",alpha=0.4,label="wt_a-wt_b");
plt.hist(y.loc[common_index]-x.loc[common_index],alpha=0.7,label="dnrp1-wt",color="green");
plt.hist(z.loc[common_index]-x.loc[common_index],alpha=0.3,label="bem1-aid-wt",color="blue");
plt.legend()

# +
bem1_aid=list_data_pd.loc["bem1-aid_merged"]
wt_merged=list_data_pd.loc["wt_merged"]
bem1_aid.index=bem1_aid.loc[:,"Gene name"]

bem1_aid.Reads.sum(),wt_merged.Reads.sum(),bem1_aid.Insertions.sum(),wt_merged.Insertions.sum()

# +
## Plotting that BEM1 is highly favored in bem1-aid background


bem1_aid=list_data_pd.loc["bem1-aid_merged"]
bem1_aid.index=bem1_aid.loc[:,"Gene name"]
x=bem1_aid.Insertions.sort_values(ascending=False).head(15)
y=bem1_aid.Reads.sort_values(ascending=False).head(15)
ax,fig=plt.subplots(1,2,figsize=(10,5))
plt.subplots_adjust(wspace=0.3)
plt.subplot(1,2,1)
plt.bar(x.index,x,color="gray")
plt.ylabel("Transposon insertions")
plt.xticks(rotation=90);
plt.subplot(1,2,2)
plt.bar(y.index,y,color="gray")
plt.ylabel("Reads")
plt.xticks(rotation=90);
plt.tight_layout()
plt.savefig("../figures/bem1_aid_bem1_overrepresented.png",dpi=300)

# +
## plot the represssion to bem1d from L Laan eLife paper from evolutionary experiments   

growth_rate_steps=[0.0082,0.001,0.006,0.007]

labels=["WT","dbem1","dbem1dbem3","dbem1dbem3dnrp1"]

plt.figure(figsize=(5,5))

#plt.plot(labels[1:],growth_rate_steps[1:],"o",color="black",markersize=10)
plt.bar(x=labels[1:],height=growth_rate_steps[1:],color="gray",alpha=0.6)
#plt.plot(growth_rate_steps[1:],growth_rate_steps[1:],"--",color="black",markersize=10,alpha=0.3)
plt.hlines(growth_rate_steps[0],0,2,linestyles="dashed",label="Wild type ",color="purple")
plt.xticks(rotation=90);
plt.ylabel("Growth rate (min$^{-1}$)")
plt.grid(linewidth=0.2)
plt.legend(loc="upper right")
plt.ylim(0,0.012)
plt.tight_layout()
plt.savefig("../figures/fitness_laan_evolution.png",dpi=300)

# +
backgrounds=["wt_a","wt_b","dnrp1_1","dnrp1_2","bem1-aid_a","bem1-aid_b"]

data=list_data_pd.loc[backgrounds]

## Plot the fold change of insertions from the two backgrounds in dnrp1 , in WT and dbem1

data_wta=data.loc["wt_a"]
data_wtb=data.loc["wt_b"]

data_bem1a=data.loc["bem1-aid_a"]
data_bem1b=data.loc["bem1-aid_b"]

# +
## Represent the fold change of insertions and reads from the two backgrounds in dnrp1 , in WT and dbem1 



data_insertions_nrp1a=data_wta[data_wta.loc[:,"Gene name"]=="NRP1"]["Insertions"].values[0]/data_wta.loc[:,"Insertions"].sum()
data_insertions_nrp1b=data_wtb[data_wtb.loc[:,"Gene name"]=="NRP1"]["Insertions"].values[0]/data_wtb.loc[:,"Insertions"].sum()

data_insertions_bem1a=data_wta[data_wta.loc[:,"Gene name"]=="BEM1"]["Insertions"].values[0]/data_wta.loc[:,"Insertions"].sum()
data_insertions_bem1b=data_wtb[data_wtb.loc[:,"Gene name"]=="BEM1"]["Insertions"].values[0]/data_wtb.loc[:,"Insertions"].sum()

data_insertions_nrp1_bem1a=data_bem1a[data_bem1a.loc[:,"Gene name"]=="NRP1"]["Insertions"].values[0]/data_bem1a.loc[:,"Insertions"].sum()
data_insertions_nrp1_bem1b=data_bem1b[data_bem1b.loc[:,"Gene name"]=="NRP1"]["Insertions"].values[0]/data_bem1b.loc[:,"Insertions"].sum()

values2plot_mean=[np.mean([data_insertions_nrp1a,data_insertions_nrp1b]),np.mean([data_insertions_bem1a,data_insertions_bem1b]),
np.mean([data_insertions_nrp1_bem1a,data_insertions_nrp1_bem1b])]

>>>>>>> 1c5527a6c3a1e3cda5126ec83fb46e9a8a798f6f
values2plot_std=[np.std([data_insertions_nrp1a,data_insertions_nrp1b]),np.std([data_insertions_bem1a,data_insertions_bem1b]),
np.std([data_insertions_nrp1_bem1a,data_insertions_nrp1_bem1b])]


plt.figure(figsize=(5,5))
plt.bar(np.arange(0,3),values2plot_mean,yerr=values2plot_std,color=["#E8A6BF70","#1EC9E0","#E8A6BF70"],alpha=0.5,capsize=5)
plt.xticks(np.arange(0,3),["dnrp1|wt","dbem1|wt","dbem1|dnrp1"],rotation=90);
plt.ylabel("# of insertions normalized to the library")
plt.grid(linewidth=0.2)
plt.tight_layout()

plt.savefig("../figures/insertions_dnrp1_dbem1.png",dpi=300)

# +
## Plot the fold change of reads from the two backgrounds in dnrp1 , in WT and dbem1

data_reads_nrp1a=data_wta[data_wta.loc[:,"Gene name"]=="NRP1"]["Reads"].values[0]/data_wta.loc[:,"Reads"].sum()
data_reads_nrp1b=data_wtb[data_wtb.loc[:,"Gene name"]=="NRP1"]["Reads"].values[0]/data_wtb.loc[:,"Reads"].sum()
<<<<<<< HEAD

data_reads_bem1a=data_wta[data_wta.loc[:,"Gene name"]=="BEM1"]["Reads"].values[0]/data_wta.loc[:,"Reads"].sum()
data_reads_bem1b=data_wtb[data_wtb.loc[:,"Gene name"]=="BEM1"]["Reads"].values[0]/data_wtb.loc[:,"Reads"].sum()

data_reads_nrp1_bem1a=data_bem1a[data_bem1a.loc[:,"Gene name"]=="NRP1"]["Reads"].values[0]/data_bem1a.loc[:,"Reads"].sum()
data_reads_nrp1_bem1b=data_bem1b[data_bem1b.loc[:,"Gene name"]=="NRP1"]["Reads"].values[0]/data_bem1b.loc[:,"Reads"].sum()

values2plot_mean=[np.mean([data_reads_nrp1a,data_reads_nrp1b]),np.mean([data_reads_bem1a,data_reads_bem1b]),
np.mean([data_reads_nrp1_bem1a,data_reads_nrp1_bem1b])]

=======

data_reads_bem1a=data_wta[data_wta.loc[:,"Gene name"]=="BEM1"]["Reads"].values[0]/data_wta.loc[:,"Reads"].sum()
data_reads_bem1b=data_wtb[data_wtb.loc[:,"Gene name"]=="BEM1"]["Reads"].values[0]/data_wtb.loc[:,"Reads"].sum()

data_reads_nrp1_bem1a=data_bem1a[data_bem1a.loc[:,"Gene name"]=="NRP1"]["Reads"].values[0]/data_bem1a.loc[:,"Reads"].sum()
data_reads_nrp1_bem1b=data_bem1b[data_bem1b.loc[:,"Gene name"]=="NRP1"]["Reads"].values[0]/data_bem1b.loc[:,"Reads"].sum()

values2plot_mean=[np.mean([data_reads_nrp1a,data_reads_nrp1b]),np.mean([data_reads_bem1a,data_reads_bem1b]),
np.mean([data_reads_nrp1_bem1a,data_reads_nrp1_bem1b])]

>>>>>>> 1c5527a6c3a1e3cda5126ec83fb46e9a8a798f6f
values2plot_std=[np.std([data_reads_nrp1a,data_reads_nrp1b]),np.std([data_reads_bem1a,data_reads_bem1b]),
np.std([data_reads_nrp1_bem1a,data_reads_nrp1_bem1b])]

plt.figure(figsize=(5,5))
plt.bar(np.arange(0,3),values2plot_mean,yerr=values2plot_std,color=["#E8A6BF70","#1EC9E0","#E8A6BF70"],alpha=0.5,capsize=5)
plt.xticks(np.arange(0,3),["dnrp1|wt","dbem1|wt","dbem1|dnrp1"],rotation=90);
plt.ylabel("# of reads normalized to the library")

plt.grid(linewidth=0.2)
plt.tight_layout()
plt.savefig("../figures/reads_dnrp1_dbem1.png",dpi=300)
# -

# ## relation with mass spec data and satay data 
<<<<<<< HEAD

# +

 ## import mass spec data and plot all changes in gene expression against the fitness 

data_mass_spec=pd.read_excel("../data/MP_CDR_ES16072021_TMT10plex_v01.xlsx",engine='openpyxl')
# change columns names to a more understanble name like : intensity_strain_nickname_replicatename: intensity_wt_a

data_mass_spec.rename(columns={"Intensity TMT10-126": "i_wt_a", 
                    "Intensity TMT10-127N": "i_wt_b",
                    "Intensity TMT10-127C": "i_dnrp1_a",
                    "Intensity TMT10-128N": "i_dnrp1_b",
                    "Intensity TMT10-128C": "i_dbem1gal_a",
                    "Intensity TMT10-129N": "i_dbem1gal_b",
                    "Intensity TMT10-129C": "i_dbem1dnrp1gal_a",
                    "Intensity TMT10-130N": "i_dbem1dnrp1gal_b",
                    "Intensity TMT10-130C": "i_dbem1galura_a",
                    "Intensity TMT10-131": "i_dbem1galura_b",
                    "Intensity 1(TMT10-126; TMT10-127N)": "i_wt_merged",
                    "Intensity 11(TMT10-127C; TMT10-128N)": "i_dnrp1_merged",
                    "Intensity 32(TMT10-128C; TMT10-129N)": "i_dbem1gal_merged",
                    "Intensity 39(TMT10-129C; TMT10-130N)": "i_dbem1dnrp1gal_merged",
                    "Intensity 40(TMT10-130C; TMT10-131)": "i_dbem1galura_merged",
                    "Ratio 1(TMT10-126; TMT10-127N)": "r_wt_merged",
                    "Ratio 11(TMT10-127C; TMT10-128N)": "r_dnrp1_merged",
                    "Ratio 32(TMT10-128C; TMT10-129N)": "r_dbem1gal_merged",
                    "Ratio 39(TMT10-129C; TMT10-130N)": "r_dbem1dnrp1gal_merged",
                    "Ratio 40(TMT10-130C; TMT10-131)": "r_dbem1galura_merged",
                    "Ratio TMT10-126": "r_wt_a", "Ratio TMT10-127N": "r_wt_b",
                    "Ratio TMT10-127C": "r_dnrp1_a","Ratio TMT10-128N": "r_dnrp1_b",
                    "Ratio TMT10-128C": "r_dbem1gal_a","Ratio TMT10-129N": "r_dbem1gal_b",
                    "Ratio TMT10-129C": "r_dbem1dnrp1gal_a","Ratio TMT10-130N": "r_dbem1dnrp1gal_b",
                    "Ratio TMT10-130C": "r_dbem1galura_a","Ratio TMT10-131": "r_dbem1galura_b"
                    },inplace=True)

data_processed=data_mass_spec
data_processed["p-value"]=10**(-data_processed["Significance"]/10)
data_processed.drop(columns=["-10lgP","Protein Group","Protein ID","PTM"],inplace=True)


data_merged_intensities=data_processed.iloc[:,22:26]

data_merged_intensities.index=data_processed.Accession

data_merged_intensities_relative_2wt=data_merged_intensities.iloc[:,0:4].div(data_merged_intensities.iloc[:,0],axis=0)

data_merged_intensities_dbem1dnrp1_relative_2dbem1=data_merged_intensities.iloc[:,0:4].div(data_merged_intensities.iloc[:,2],axis=0)

data2plot=data_merged_intensities_relative_2wt.iloc[:,0:3].copy()
data2plot["i_dbem1dnrp12dbem1"]=data_merged_intensities_dbem1dnrp1_relative_2dbem1.iloc[:,3]


# +
mapping_accesion=pd.read_csv("../postprocessed-data/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.03.27-15.14.21.67.tsv",
            sep="\t")
accesion_uniprot=mapping_accesion.Entry+"|"+mapping_accesion.loc[:,"Entry Name"]

mapping_accesion["Entry Name"]=accesion_uniprot

mapping_accesion.set_index("Entry Name",inplace=True)

# +
mapping_index=[]

for i in data_processed.Accession:
    if i in mapping_accesion.index:
        mapping_index.append(mapping_accesion.loc[i,"Gene Names (primary)"])
        
    else:
        mapping_index.append(i)
        

data_processed.index=mapping_index
data2plot.index=mapping_index

# +
# accesion_numbers2uniprot=[]

# for i in data_processed.Accession:
#     accesion_numbers2uniprot.append(i.strip("|").split("|")[0]+" "+i.strip("|").split("|")[1].split("_")[1])

# file = open("../postprocessed-data/genes_mass_spec.txt","w")
# for item in mapping_index:
#     if type(item)==str:
#         file.write(item+"\n")
# file.close()

# +

dnrp1_a_fitness=data_fitness.loc["wt_a"].loc["NRP1","fitness_gene"]
dnrp1_b_fitness=data_fitness.loc["wt_b"].loc["NRP1","fitness_gene"]
dnrp1_fitness=(dnrp1_a_fitness+dnrp1_b_fitness)/2
fitness=data_fitness.loc["dnrp1_merged"]
index_mass_spec_fitnes=data2plot.index.intersection(fitness.index)

len(index_mass_spec_fitnes)

## normalize fitness to the median of the growth rate of dnrp1 in wt
fitness_average=(fitness.loc[:,"fitness_gene"])*dnrp1_fitness/np.median(fitness.loc[:,"fitness_gene"])

fitness_domains=(fitness.loc[:,"fitness_domains_corrected"])*dnrp1_fitness/np.median(fitness.loc[:,"fitness_domains_corrected"])



=======

# +

 ## import mass spec data and plot all changes in gene expression against the fitness 

data_mass_spec=pd.read_excel("../data/MP_CDR_ES16072021_TMT10plex_v01.xlsx",engine='openpyxl')
# change columns names to a more understanble name like : intensity_strain_nickname_replicatename: intensity_wt_a

data_mass_spec.rename(columns={"Intensity TMT10-126": "i_wt_a", 
                    "Intensity TMT10-127N": "i_wt_b",
                    "Intensity TMT10-127C": "i_dnrp1_a",
                    "Intensity TMT10-128N": "i_dnrp1_b",
                    "Intensity TMT10-128C": "i_dbem1gal_a",
                    "Intensity TMT10-129N": "i_dbem1gal_b",
                    "Intensity TMT10-129C": "i_dbem1dnrp1gal_a",
                    "Intensity TMT10-130N": "i_dbem1dnrp1gal_b",
                    "Intensity TMT10-130C": "i_dbem1galura_a",
                    "Intensity TMT10-131": "i_dbem1galura_b",
                    "Intensity 1(TMT10-126; TMT10-127N)": "i_wt_merged",
                    "Intensity 11(TMT10-127C; TMT10-128N)": "i_dnrp1_merged",
                    "Intensity 32(TMT10-128C; TMT10-129N)": "i_dbem1gal_merged",
                    "Intensity 39(TMT10-129C; TMT10-130N)": "i_dbem1dnrp1gal_merged",
                    "Intensity 40(TMT10-130C; TMT10-131)": "i_dbem1galura_merged",
                    "Ratio 1(TMT10-126; TMT10-127N)": "r_wt_merged",
                    "Ratio 11(TMT10-127C; TMT10-128N)": "r_dnrp1_merged",
                    "Ratio 32(TMT10-128C; TMT10-129N)": "r_dbem1gal_merged",
                    "Ratio 39(TMT10-129C; TMT10-130N)": "r_dbem1dnrp1gal_merged",
                    "Ratio 40(TMT10-130C; TMT10-131)": "r_dbem1galura_merged",
                    "Ratio TMT10-126": "r_wt_a", "Ratio TMT10-127N": "r_wt_b",
                    "Ratio TMT10-127C": "r_dnrp1_a","Ratio TMT10-128N": "r_dnrp1_b",
                    "Ratio TMT10-128C": "r_dbem1gal_a","Ratio TMT10-129N": "r_dbem1gal_b",
                    "Ratio TMT10-129C": "r_dbem1dnrp1gal_a","Ratio TMT10-130N": "r_dbem1dnrp1gal_b",
                    "Ratio TMT10-130C": "r_dbem1galura_a","Ratio TMT10-131": "r_dbem1galura_b"
                    },inplace=True)

data_processed=data_mass_spec
data_processed["p-value"]=10**(-data_processed["Significance"]/10)
data_processed.drop(columns=["-10lgP","Protein Group","Protein ID","PTM"],inplace=True)


data_merged_intensities=data_processed.iloc[:,22:26]

data_merged_intensities.index=data_processed.Accession

data_merged_intensities_relative_2wt=data_merged_intensities.iloc[:,0:4].div(data_merged_intensities.iloc[:,0],axis=0)

data_merged_intensities_dbem1dnrp1_relative_2dbem1=data_merged_intensities.iloc[:,0:4].div(data_merged_intensities.iloc[:,2],axis=0)

data2plot=data_merged_intensities_relative_2wt.iloc[:,0:3].copy()
data2plot["i_dbem1dnrp12dbem1"]=data_merged_intensities_dbem1dnrp1_relative_2dbem1.iloc[:,3]


# +
mapping_accesion=pd.read_csv("../postprocessed-data/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.03.27-15.14.21.67.tsv",
            sep="\t")
accesion_uniprot=mapping_accesion.Entry+"|"+mapping_accesion.loc[:,"Entry Name"]

mapping_accesion["Entry Name"]=accesion_uniprot

mapping_accesion.set_index("Entry Name",inplace=True)

# +
mapping_index=[]

for i in data_processed.Accession:
    if i in mapping_accesion.index:
        mapping_index.append(mapping_accesion.loc[i,"Gene Names (primary)"])
        
    else:
        mapping_index.append(i)
        

data_processed.index=mapping_index
data2plot.index=mapping_index

# +
# accesion_numbers2uniprot=[]

# for i in data_processed.Accession:
#     accesion_numbers2uniprot.append(i.strip("|").split("|")[0]+" "+i.strip("|").split("|")[1].split("_")[1])

# file = open("../postprocessed-data/genes_mass_spec.txt","w")
# for item in mapping_index:
#     if type(item)==str:
#         file.write(item+"\n")
# file.close()

# +
fitness=data_fitness.loc["dnrp1_merged"]
index_mass_spec_fitnes=data2plot.index.intersection(fitness.index)

len(index_mass_spec_fitnes)

## normalize fitness from 0 to 1
fitness_average=(fitness.loc[:,"fitness_gene"])
#fitness_average=(fitness_average-np.min(fitness_average))/(np.max(fitness_average)-np.min(fitness_average))

fitness_average=fitness_average
fitness_domains=(fitness.loc[:,"fitness_domains_corrected"])

#fitness_domains=(fitness_domains-np.min(fitness_domains))/(np.max(fitness_domains)-np.min(fitness_domains))
fitness_domains=fitness_domains
>>>>>>> 1c5527a6c3a1e3cda5126ec83fb46e9a8a798f6f

# +
x=data2plot.loc[index_mass_spec_fitnes,"i_dnrp1_merged"]
over_index=x[x>1].index
under_index=x[x<1].index
equal_index=x[x==1].index

relative_expression_dnrp1_wt=x
# -

<<<<<<< HEAD
relative_expression_dnrp1_wt["HOM3"],relative_expression_dnrp1_wt["BAT1"]

fitness_average.loc[over_index].sort_values(ascending=False).head(20)

fitness_average.loc[under_index].sort_values(ascending=True)

# +
over_mean=fitness_average.loc[over_index].mean()
over_std=fitness_average.loc[over_index].std()/np.sqrt(len(over_index))

under_mean=fitness_average.loc[under_index].mean()
under_std=fitness_average.loc[under_index].std()/np.sqrt(len(under_index))

print(over_mean,over_std,under_mean,under_std)

# +
plt.figure(figsize=(5,5))
plt.scatter(fitness_average.loc[over_index],data2plot.loc[over_index,"i_dnrp1_merged"],
c="purple",alpha=0.4,s=50,vmin=0,vmax=1.2)

plt.scatter(fitness_average.loc[under_index],data2plot.loc[under_index,"i_dnrp1_merged"],
c="green",alpha=0.4,s=50,vmin=0,vmax=1.2)

plt.scatter(fitness_average.loc[equal_index],data2plot.loc[equal_index,"i_dnrp1_merged"],
c="grey",alpha=0.4,s=100,vmin=0,vmax=1.2)
plt.hlines(1,0,1.75,linestyles="dashed",color="k")

plt.errorbar(over_mean,1.75,xerr=over_std,fmt="o",color="k",capsize=5)
plt.errorbar(under_mean,0.5,xerr=under_std,fmt="o",color="k",capsize=5)

plt.grid(linewidth=0.5)
plt.xlabel("fitness")
plt.ylabel("relative intensity to WT")
plt.tight_layout()
plt.savefig("../figures/fitness_vs_expression_dnrp1.png",dpi=300)

# +
## To see genes over and underexpressed in dbem1dnrp1 
data2plot=data_processed.loc[:, data_processed.columns.str.startswith('r_')]

data2plot=data2plot.loc[:, data2plot.columns.str.endswith('merged')]
data2plot=data2plot.iloc[:,3:5] # to compare with dbem1+URA
data2plot.loc[:,"r_dbem1dnrp1_vs_dbem1"]=data2plot.iloc[:,0]/data2plot.iloc[:,1]
data2plot.index=mapping_index


### Taking the fitness in dbem1
dbem1_a_fitness=data_fitness.loc["wt_a"].loc["BEM1","fitness_gene"]
dbem1_b_fitness=data_fitness.loc["wt_b"].loc["BEM1","fitness_gene"]
dbem1_fitness=(dbem1_a_fitness+dbem1_b_fitness)/2
fitness_dbem1=data_fitness.loc["bem1-aid_merged"]
index_mass_spec_fitnes=data2plot.index.intersection(fitness_dbem1.index)

## normalize fitness to the median of the fitness
fitness_average=(fitness_dbem1.loc[:,"fitness_gene"])*dbem1_fitness/np.median(fitness_dbem1.loc[:,"fitness_gene"])
#fitness_average=(fitness_average-np.min(fitness_average))/(np.max(fitness_average)-np.min(fitness_average))

#fitness_average=fitness_average/np.median(fitness_average)


fitness_domains=(fitness_dbem1.loc[:,"fitness_domains_corrected"])*dbem1_fitness/np.median(fitness_dbem1.loc[:,"fitness_domains_corrected"])

#fitness_domains=(fitness_domains-np.min(fitness_domains))/(np.max(fitness_domains)-np.min(fitness_domains))
#fitness_domains=fitness_domains/np.median(fitness_domains)
### 
x=data2plot.loc[index_mass_spec_fitnes,"r_dbem1dnrp1_vs_dbem1"]
over_index=x[x>1].index
under_index=x[x<1].index
equal_index=x[x==1].index
relative_expression_dnrp1_dbem1=x
# -

x=data2plot.loc[index_mass_spec_fitnes,"r_dbem1dnrp1_vs_dbem1"]
x["HOM3"],x["BAT1"]

fitness_average.loc[over_index].sort_values(ascending=False).head(10),fitness_average.loc[under_index].sort_values(ascending=False).head(20)

fitness_average.loc["MDH1"],x["MDH1"]

# +
over_mean=fitness_average.loc[over_index].mean()
over_std=fitness_average.loc[over_index].std()/np.sqrt(len(over_index))

under_mean=fitness_average.loc[under_index].mean()
under_std=fitness_average.loc[under_index].std()/np.sqrt(len(under_index))

print(over_mean,over_std,under_mean,under_std)

# +
plt.figure(figsize=(5,5))
plt.scatter(fitness_average.loc[over_index],data2plot.loc[over_index,"r_dbem1dnrp1_vs_dbem1"],
c="purple",alpha=0.4,s=50,vmin=0,vmax=1.2)

plt.scatter(fitness_average.loc[under_index],data2plot.loc[under_index,"r_dbem1dnrp1_vs_dbem1"],
c="green",alpha=0.4,s=50,vmin=0,vmax=1.2)

plt.scatter(fitness_average.loc[equal_index],data2plot.loc[equal_index,"r_dbem1dnrp1_vs_dbem1"],
c="grey",alpha=0.4,s=100,vmin=0,vmax=1.2)

plt.errorbar(over_mean,1.3,xerr=over_std,fmt="o",color="black",capsize=5)
plt.errorbar(under_mean,0.7,xerr=under_std,fmt="o",color="black",capsize=5)

plt.hlines(1,-10,20,linestyles="dashed",color="k")

plt.grid(linewidth=0.5)
plt.xlim(-5,20)
plt.xlabel("fitness")
plt.ylabel("relative intensity dbem1")
plt.tight_layout()
plt.savefig("../figures/fitness_vs_expression_dnrp1dbem1_dbem1_URA.png",dpi=300)
=======
x=data2plot.loc[index_mass_spec_fitnes,"i_dnrp1_merged"]
over_index=x[x>1].index
under_index=x[x<1].index
equal_index=x[x==1].index

fitness_average.loc[over_index].sort_values(ascending=False).head(20),fitness_average.loc["RIB4"]

fitness_average.loc[under_index].sort_values(ascending=True).head(10)

# +
over_mean=fitness_average.loc[over_index].mean()
over_std=fitness_average.loc[over_index].std()/np.sqrt(len(over_index))

under_mean=fitness_average.loc[under_index].mean()
under_std=fitness_average.loc[under_index].std()/np.sqrt(len(under_index))

print(over_mean,over_std,under_mean,under_std)

# +
plt.figure(figsize=(5,5))
plt.scatter(fitness_average.loc[over_index],data2plot.loc[over_index,"i_dnrp1_merged"],
c="purple",alpha=0.4,s=50,vmin=0,vmax=1.2)

plt.scatter(fitness_average.loc[under_index],data2plot.loc[under_index,"i_dnrp1_merged"],
c="green",alpha=0.4,s=50,vmin=0,vmax=1.2)

plt.scatter(fitness_average.loc[equal_index],data2plot.loc[equal_index,"i_dnrp1_merged"],
c="grey",alpha=0.4,s=100,vmin=0,vmax=1.2)
plt.hlines(1,0,1.5,linestyles="dashed",color="k")

plt.errorbar(over_mean,1.75,xerr=over_std,fmt="o",color="k",capsize=5)
plt.errorbar(under_mean,0.5,xerr=under_std,fmt="o",color="k",capsize=5)

plt.grid(linewidth=0.5)
plt.xlabel("fitness")
plt.ylabel("relative intensity to WT")
plt.tight_layout()
plt.savefig("../figures/fitness_vs_expression_dnrp1.png",dpi=300)

# +
plt.figure(figsize=(5,5))
plt.scatter(fitness_domains.loc[over_index],data2plot.loc[over_index,"i_dnrp1_merged"],
c="purple",alpha=0.4,s=50,vmin=0,vmax=1.2)

plt.scatter(fitness_domains.loc[under_index],data2plot.loc[under_index,"i_dnrp1_merged"],
c="green",alpha=0.4,s=50,vmin=0,vmax=1.2)

plt.scatter(fitness_domains.loc[equal_index],data2plot.loc[equal_index,"i_dnrp1_merged"],
c="grey",alpha=0.4,s=100,vmin=0,vmax=1.2)
plt.hlines(1,0,1.25,linestyles="dashed",color="k")

plt.grid(linewidth=0.5)
plt.xlabel("fitness")
plt.ylabel("relative intensity dnrp1")
plt.tight_layout()
#plt.savefig("../figures/fitness_vs_expression_dnrp1.png",dpi=300)

# +
## To see genes over and underexpressed in dbem1dnrp1 
data2plot=data_processed.loc[:, data_processed.columns.str.startswith('r_')]

data2plot=data2plot.loc[:, data2plot.columns.str.endswith('merged')]
data2plot=data2plot.iloc[:,3:5] # to compare with dbem1+URA
data2plot.loc[:,"r_dbem1dnrp1_vs_dbem1"]=data2plot.iloc[:,0]/data2plot.iloc[:,1]
data2plot.index=mapping_index


### Taking the fitness in dbem1

fitness_dbem1=data_fitness.loc["bem1-aid_merged"]
index_mass_spec_fitnes=data2plot.index.intersection(fitness_dbem1.index)

## normalize fitness from 0 to 1
fitness_average=(fitness_dbem1.loc[:,"fitness_gene"])
#fitness_average=(fitness_average-np.min(fitness_average))/(np.max(fitness_average)-np.min(fitness_average))

#fitness_average=fitness_average/np.median(fitness_average)


fitness_domains=(fitness_dbem1.loc[:,"fitness_domains_corrected"])

#fitness_domains=(fitness_domains-np.min(fitness_domains))/(np.max(fitness_domains)-np.min(fitness_domains))
#fitness_domains=fitness_domains/np.median(fitness_domains)
### 
x=data2plot.loc[index_mass_spec_fitnes,"r_dbem1dnrp1_vs_dbem1"]
over_index=x[x>1].index
under_index=x[x<1].index
equal_index=x[x==1].index
# -

fitness_average.loc[over_index].sort_values(ascending=False).head(5),fitness_average.loc[under_index].sort_values(ascending=False).head(10)

fitness_average.loc["MDH1"],x["MDH1"]
>>>>>>> 1c5527a6c3a1e3cda5126ec83fb46e9a8a798f6f

# +
over_mean=fitness_average.loc[over_index].mean()
over_std=fitness_average.loc[over_index].std()/np.sqrt(len(over_index))

under_mean=fitness_average.loc[under_index].mean()
under_std=fitness_average.loc[under_index].std()/np.sqrt(len(under_index))

print(over_mean,over_std,under_mean,under_std)

# +
<<<<<<< HEAD
over_mean=fitness_domains.loc[over_index].mean()
over_std=fitness_domains.loc[over_index].std()/np.sqrt(len(over_index))

under_mean=fitness_domains.loc[under_index].mean()
under_std=fitness_domains.loc[under_index].std()/np.sqrt(len(under_index))

print(over_mean,over_std,under_mean,under_std)
# -

x=fitness_average.loc[under_index]
x[x>4]

y=fitness_domains.loc[over_index]
y[y>2]


# +
x=list_data_pd.loc["wt_merged"]
y=list_data_pd.loc["bem1-aid_merged"]
z=list_data_pd.loc["dnrp1_merged"]

nrp1_insertions=x[x.loc[:,"Gene name"]=="NRP1"].loc[:,"Insertions"]
nrp1_insertions_dbem1=y[y.loc[:,"Gene name"]=="NRP1"].loc[:,"Insertions"]

bem1_insertions=x[x.loc[:,"Gene name"]=="BEM1"].loc[:,"Insertions"]
bem1_insertions_dnrp1=z[z.loc[:,"Gene name"]=="BEM1"].loc[:,"Insertions"]


data=(data_fitness.loc["wt_merged"])
data_average=(data.loc[:,"fitness_gene"])
data_norm=(data_average-data_average.min())/(data_average.max()-data_average.min())

mean_bar=data_average.loc["NRP1"]
std_error=data.loc["NRP1","fitness_gene_std"]/np.sqrt(nrp1_insertions) # merging of two samples that contain in total 56 insertions in nrp1,rendering 56 different replicates

data_y=(data_fitness.loc["bem1-aid_merged"])
data_y_average=(data_y.loc[:,"fitness_gene"])*0.5*(data.loc["BEM1","fitness_gene"]+ data.loc["BEM1","fitness_domains_corrected"])/np.median(data_y.loc[:,"fitness_gene"])
data_y_norm=(data_y_average-data_y_average.min())/(data_y_average.max()-data_y_average.min())

mean_bar_y=data_y_average.loc["NRP1"]
std_error_y=data_y.loc["NRP1","fitness_gene_std"]/np.sqrt(nrp1_insertions_dbem1) # merging of two samples that contain in total 56 insertions in nrp1,rendering 56 different replicates

mean_bar_bem1=0.5*(data_average.loc["BEM1"]+data.loc["BEM1","fitness_domains_corrected"])
std_error_bem1=data.loc["BEM1","fitness_gene_std"]/np.sqrt(bem1_insertions) 

data_z=data_fitness.loc["dnrp1_merged"]
mean_bar_z=data_z.loc["BEM1","fitness_gene"]
std_error_z=data_z.loc["BEM1","fitness_gene_std"]/np.sqrt(bem1_insertions_dnrp1)






# +
## Plot fitness of nrp1 taking the fitness average  as an error bar plot



## plotting in a error bar plot

plt.figure(figsize=(8,5))

plt.bar(["fitness nrp1 in WT","fitness bem1 in WT","fitness nrp1 in dbem1"],
[mean_bar,mean_bar_bem1,mean_bar_y],
yerr=[std_error.values[0],std_error_bem1.values[0],std_error_y.values[0]],
color=["#E8A6BF70","#1EC9E0","#E8A6BF70"],capsize=5,alpha=0.5)

# write every mean value on top of every background
# plt.text(0,mean_bar+0.5,"{:.2f}".format(mean_bar),ha="center",va="center",color="black")
# plt.text(1,mean_bar_bem1+0.5,"{:.2f}".format(mean_bar_bem1),ha="center",va="center",color="black")
# plt.text(2,mean_bar_y+0.5,"{:.2f}".format(mean_bar_y),ha="center",va="center",color="black")



plt.ylabel("Fitness from SATAY")
plt.grid(linewidth=0.2)
plt.ylabel("Fitness normalized")
# plt.yticks(np.arange(0,30,step=2));

plt.tight_layout()

plt.savefig("../figures/fitness_dnrp1_dbem1.png",dpi=300)

# +
## take the same index across all libraries

index_arrays=np.array([data_fitness.loc["wt_merged"].index,data_fitness.loc["bem1-aid_merged"].index,
data_fitness.loc["bem1-aid-dbem3_a"].index],
dtype="object")

d1=set.intersection(*map(set,index_arrays))
fitness_bem1_wt=(data_fitness.loc["wt_merged"].loc["BEM1","fitness_gene"]+data_fitness.loc["wt_merged"].loc["BEM1","fitness_domains_corrected"])/2

data_wt=data_fitness.loc["wt_merged"].loc[d1,"fitness_gene"]
data_dbem1=data_fitness.loc["bem1-aid_merged"].loc[d1,"fitness_gene"]

## Normalize the datasetss such as the median values corresponds to the fitness of the knockout of the gene of interest 
data_wt_norm=data_wt

data_dbem1_norm=data_dbem1*fitness_bem1_wt/(np.median(data_dbem1))# because their median is around 0 

# +
## Plot distributions of fitness landscapes in WT and dbem1 

plt.figure(figsize=(6,5))

## Normalize fitness .values()


# plt.hist(data_y_norm/data_y_norm.median(),bins=20,color="#1EC9E0",alpha=0.5,label="$\Delta$bem1");
# plt.hist(data_norm/data_norm.median(),bins=20,color="gray",alpha=0.8,label="WT");
plt.boxplot([data_dbem1_norm,data_wt],labels=["$\Delta$bem1","WT"],
showfliers=True,patch_artist=False)

plt.yticks([-3,0,0.4,1,5]);

plt.ylabel("Fitness normalized")
plt.grid(linewidth=0.5,linestyle="--",axis="y")

plt.ylim(-3,5)


plt.tight_layout()
plt.savefig("../figures/fitness_distributions_dbem1_wt.png",dpi=300)

# +
## Plot insertion map onto gene of interest


background="wt_merged"
data=list_data_pd.loc[background]
## remove ade2 data from the dataset
# data.drop(labels=5805,axis=0,inplace=True)
# ## remove the ura3 data from the dataset
# data.drop(labels=1927,axis=0,inplace=True)


## Extract as lists the gene coordinates, and the locations of every insertion and reads per gene 
reads_locations=[]
insertion_locations=[]
gene_coordinates=[]

for j in data.index:
        coi=from_excel_to_list(data.loc[j]["Reads per insertion location"])
        coi_insertions=from_excel_to_list(data.loc[j]["Insertion locations"])
        gene_coordinates.append([data.loc[j,"Start location"],data.loc[j,"End location"]])
        reads_locations.append(coi)
        insertion_locations.append(coi_insertions)
        


# compute the reads per insertions for each gene along the gene length

gene_parts=np.linspace(0,1,11)
r=np.zeros(shape=(len(insertion_locations),len(gene_parts))) # reads_per_insertion_parts array , every gene in the rows 

for i in np.arange(0,len(insertion_locations)):
    if (insertion_locations[i])!=0:
        g=np.array(insertion_locations[i]) # insertion locations
        f=np.linspace(gene_coordinates[i][0],gene_coordinates[i][1],len(gene_parts)) # gene coordinates
        #f=np.array(gene_coordinates[i][0]+gene_coordinates[i][1]*gene_parts)
        binedges = g.searchsorted(f)

        rngs = [list(range(binedges[k], binedges[k+1])) for k in range(len(binedges)-1)]

        for k in np.arange(0,len(rngs)):
            readsperinsert=[]
            for j in np.arange(0,len(rngs[k])):
                readsperinsert.append(reads_locations[i][rngs[k][j]])
                
            if len(readsperinsert)>1:
                r[i,k]=np.sum(readsperinsert)/(len(readsperinsert)-1)#discarding the insertion with the highest read count
            else:
                r[i,k]=0

    
=======
plt.figure(figsize=(5,5))
plt.scatter(fitness_average.loc[over_index],data2plot.loc[over_index,"r_dbem1dnrp1_vs_dbem1"],
c="purple",alpha=0.4,s=50,vmin=0,vmax=1.2)

plt.scatter(fitness_average.loc[under_index],data2plot.loc[under_index,"r_dbem1dnrp1_vs_dbem1"],
c="green",alpha=0.4,s=50,vmin=0,vmax=1.2)

plt.scatter(fitness_average.loc[equal_index],data2plot.loc[equal_index,"r_dbem1dnrp1_vs_dbem1"],
c="grey",alpha=0.4,s=100,vmin=0,vmax=1.2)

plt.errorbar(over_mean,1.3,xerr=over_std,fmt="o",color="black",capsize=5)
plt.errorbar(under_mean,0.7,xerr=under_std,fmt="o",color="black",capsize=5)
>>>>>>> 1c5527a6c3a1e3cda5126ec83fb46e9a8a798f6f

plt.hlines(1,-10,25,linestyles="dashed",color="k")

<<<<<<< HEAD

# +
nrp1_wt_index=np.where(data.loc[:,"Gene name"]=="NRP1")[0]
nrp1_coordinates=np.linspace(gene_coordinates[nrp1_wt_index[0]][0],gene_coordinates[nrp1_wt_index[0]][1],11)
nrp1_coordinates2text=["{:.0f}".format(i) for i in nrp1_coordinates]



# +
plt.figure(figsize=(10,5))

plt.bar(gene_parts,r[nrp1_wt_index[0],:],width=0.1)
plt.xticks(np.linspace(-0.05,0.95,11),
labels=nrp1_coordinates2text,rotation=90);
# -

bem1d=list_data_pd.loc["bem1-aid_merged"]
bem1d.index=bem1d["Gene name"]
bem1d

# +


reads_nrp1_in_dbem1=from_excel_to_list(bem1d.loc["NRP1","Reads per insertion location"])

plt.bar(np.arange(0,len(reads_nrp1_in_dbem1)),reads_nrp1_in_dbem1)
# -

# ### Digenic interaction for Nrp1

gi_nrp1=digenic_GI(data=data_fitness,goi="NRP1",col_fitness="fitness_gene",backg=["wt_a","wt_b","dnrp1_1","dnrp1_2"],significance=0.05)

# +
gi_pd_fitness_gene=pd.DataFrame.from_dict(gi_nrp1,orient="index")

gi_pd_fitness_gene["gene_names"]=gi_pd_fitness_gene.index

gi_pd_fitness_gene
# -

gi_pd_fitness_gene.loc[:,"fold_change"].hist(bins=20)

gi_pd_fitness_gene.loc["KAR3"],gi_pd_fitness_gene.loc["CIN8"]

# +
## import genes from the cellcycle pathway and mapk pathway from wikipathway

cellcycle=pd.read_csv("../postprocessed-data/Cell cycle and cell division - Saccharomyces cerevisiae default  node.csv",sep="\t")

cellcycle=cellcycle.name.tolist()

mapk=pd.read_csv("../postprocessed-data/MAPK signaling-pathway-Saccharomyces-cerevisiae-default-node.csv",sep="\t")

mapk=mapk.name.tolist()

transcription_initiation=pd.read_csv("../postprocessed-data/Eukaryotic transcription initiation - Saccharomyces cerevisiae default  node.csv",sep="\t")

transcription_initiation=transcription_initiation.name.tolist()

lipid=pd.read_csv("../postprocessed-data/Sphingolipid metabolism - Saccharomyces cerevisiae default  node.csv",sep="\t")
lipid=lipid.name.tolist()

lipid_sugar=pd.read_csv("../postprocessed-data/Lipid-linked oligosaccharide biosynthesis - Saccharomyces cerevisiae default  node.csv",sep=",")
lipid_sugar=lipid_sugar.name.tolist()

mating=pd.read_csv("../postprocessed-data/Mating-pheromone response pathway - Saccharomyces cerevisiae default  node.csv",sep=",")
mating=mating.name.tolist()

glucose_repression=pd.read_csv("../postprocessed-data/Glucose repression - Saccharomyces cerevisiae default  node.csv",sep=",")
glucose_repression=glucose_repression.name.tolist()
# -

all_pathways=cellcycle+mapk+transcription_initiation+lipid+lipid_sugar+mating+glucose_repression
all_pathways=list(set(all_pathways))

data2plot

# +

# make all the gene names uppercase

all_pathways=[i.upper() for i in all_pathways]

all_pathways_gi=gi_pd_fitness_gene[gi_pd_fitness_gene.loc[:,"gene_names"].isin(all_pathways)]

for i in  all_pathways_gi.index:
    if i in data2plot.index:
        all_pathways_gi.loc[i,"expression_dnrp1"]=data2plot.loc[i,"i_dnrp1_merged"]
        all_pathways_gi.loc[i,"expression_dnrp1dbem1"]=data2plot.loc[i,"i_dbem1dnrp12dbem1"]

all_pathways_gi.fillna(0,inplace=True)

all_pathways_gi=all_pathways_gi.loc[:,["gene_names","fold_change","p_statistic","expression_dnrp1","expression_dnrp1dbem1"]]

all_pathways_gi.to_csv("../postprocessed-data/all_pathways_gi_nrp1.csv",sep="\t")
# cellcycle_mapk_gi.to_csv("../postprocessed-data/cellcycle_mapk_gi_nrp1.csv",sep="\t")
#transcription_initiation_gi.to_csv("../postprocessed-data/transcription_initiation_gi_nrp1.csv",sep="\t")
# -

all_pathways_gi.sort_values(by="expression_dnrp1",ascending=False)

# +
under_index_wt=relative_expression_dnrp1_wt[relative_expression_dnrp1_wt<1].index
over_index_wt=relative_expression_dnrp1_wt[relative_expression_dnrp1_wt>1].index



# +
gi_under_index_wt=[]
for i in under_index_wt:
    if i in gi_pd_fitness_gene.index:
        gi_under_index_wt.append(i)

gi_over_index_wt=[]
for i in over_index_wt:
    if i in gi_pd_fitness_gene.index:
        gi_over_index_wt.append(i)

print("The interaction scores for nrp1 for downregulated genes in WT range from:",gi_pd_fitness_gene.loc[gi_under_index_wt].fold_change.mean()
+gi_pd_fitness_gene.loc[gi_under_index_wt].fold_change.std(),"to",gi_pd_fitness_gene.loc[gi_under_index_wt].fold_change.mean()
-gi_pd_fitness_gene.loc[gi_under_index_wt].fold_change.std())

print("The interaction scores for nrp1 for upregulated genes in WT range from:",gi_pd_fitness_gene.loc[gi_over_index_wt].fold_change.mean()
+gi_pd_fitness_gene.loc[gi_over_index_wt].fold_change.std(),"to",gi_pd_fitness_gene.loc[gi_over_index_wt].fold_change.mean()
-gi_pd_fitness_gene.loc[gi_over_index_wt].fold_change.std())

# +
gi_under_index=[]
for i in under_index:
    if i in gi_pd_fitness_gene.index:
        gi_under_index.append(i)

gi_over_index=[]
for i in over_index:
    if i in gi_pd_fitness_gene.index:
        gi_over_index.append(i)

print("The interaction scores for nrp1 for downregulated genes in dbem1dnrp1 range from:",gi_pd_fitness_gene.loc[gi_under_index].fold_change.mean()
+gi_pd_fitness_gene.loc[gi_under_index].fold_change.std(),"to",gi_pd_fitness_gene.loc[gi_under_index].fold_change.mean()
-gi_pd_fitness_gene.loc[gi_under_index].fold_change.std())

print("The interaction scores for nrp1 for upregulated genes in dbem1dnrp1 range from:",gi_pd_fitness_gene.loc[gi_over_index].fold_change.mean()
+gi_pd_fitness_gene.loc[gi_over_index].fold_change.std(),"to",gi_pd_fitness_gene.loc[gi_over_index].fold_change.mean()
-gi_pd_fitness_gene.loc[gi_over_index].fold_change.std())
# -

relative_expression_dnrp1_dbem1[under_index].sort_values(ascending=True)[0:5]

relative_expression_dnrp1_dbem1["GLY1"]

gi_pd_fitness_gene.loc[gi_over_index].sort_values(by="fold_change",ascending=True)

# +
plt.boxplot([gi_pd_fitness_gene.loc[gi_under_index_wt].fold_change,gi_pd_fitness_gene.loc[gi_over_index_wt].fold_change,
gi_pd_fitness_gene.loc[gi_under_index].fold_change,gi_pd_fitness_gene.loc[gi_over_index].fold_change],
labels=["Underexpressed dnrp1","Overexpressed dnrp1","Underexpressed dbem1dnrp1","Overexpressed dbem1dnrp1"]);

plt.xticks(rotation=90);
# -

a=relative_expression_dnrp1_wt[gi_over_index_wt].sort_values(ascending=False)[0:5]
gi_pd_fitness_gene.loc[a.index]

# +
plt.figure(figsize=(8,6))
plt.scatter(relative_expression_dnrp1_wt[gi_under_index_wt],gi_pd_fitness_gene.loc[gi_under_index_wt].fold_change,c=gi_pd_fitness_gene.loc[gi_under_index_wt].p_statistic,
s=100,cmap="Blues_r",vmax=1,vmin=0)
plt.hlines(y=0,xmin=0.4,xmax=1,linestyles="dashed",color="black")


plt.colorbar(label="p-value")

plt.xlabel("Expression ratio dnrp1/wt",fontsize=18)
plt.ylabel("Interaction score with Nrp1",fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.tight_layout()

#plt.savefig("../figures/interaction_score_vs_expression_ratio_dnrp1_wt_underexpressed.png",dpi=300,transparent=True)
# -

gi_pd_fitness_gene.loc[gi_under_index_wt].fold_change.mean()

# +
plt.figure(figsize=(8,6))
plt.scatter(relative_expression_dnrp1_dbem1[gi_under_index],gi_pd_fitness_gene.loc[gi_under_index].fold_change,c=gi_pd_fitness_gene.loc[gi_under_index].p_statistic,
s=100,cmap="Blues_r",vmax=1,vmin=0)
plt.hlines(y=0,xmin=0.7,xmax=1,linestyles="dashed",color="black")

plt.colorbar(label="p-value")

plt.xlabel("Expression ratio bem1dnrp1/bem1")
plt.ylabel("Interaction score eith Nrp1")

plt.tight_layout()

#plt.savefig("../figures/interaction_score_vs_expression_ratio_bem1dnrp1_bem1_underexpressed.png",dpi=300,transparent=True)

# +
plt.figure(figsize=(8,6))
plt.scatter(relative_expression_dnrp1_wt[gi_over_index_wt],gi_pd_fitness_gene.loc[gi_over_index_wt].fold_change,c=gi_pd_fitness_gene.loc[gi_over_index_wt].p_statistic,
s=100,cmap="Blues_r",vmax=1,vmin=0)
plt.hlines(y=0,xmin=1,xmax=2.2,linestyles="dashed",color="black")

plt.colorbar(label="p-value")

plt.xlabel("Expression ratio dnrp1/wt")
plt.ylabel("Interaction score eith Nrp1")

plt.tight_layout()

#plt.savefig("../figures/interaction_score_vs_expression_ratio_dnrp1_wt_overexpressed.png",dpi=300,transparent=True)

# +
plt.figure(figsize=(8,6))
plt.scatter(relative_expression_dnrp1_dbem1[gi_over_index],gi_pd_fitness_gene.loc[gi_over_index].fold_change,c=gi_pd_fitness_gene.loc[gi_over_index].p_statistic,
s=100,cmap="Blues_r",vmax=1,vmin=0)
plt.hlines(y=0,xmin=1,xmax=1.4,linestyles="dashed",color="black")

plt.colorbar(label="p-value")

plt.xlabel("Expression ratio bem1dnrp1/bem1")
plt.ylabel("Interaction score eith Nrp1")

plt.tight_layout()

plt.savefig("../figures/interaction_score_vs_expression_ratio_bem1dnrp1_bem1_overderexpressed.png",dpi=300,transparent=True)

# +
## Boxplot of the normalized fold change
plt.figure(figsize=(5,5))
x=gi_pd_fitness_gene[gi_pd_fitness_gene.loc[:,"significance"]==True]["fold_change"]
x_m=np.abs(np.mean(x))
x_std=np.std(x)
plt.hist(gi_pd_fitness_gene.loc[:,"fold_change"],bins=100,histtype="step",color="black");
# plt.vlines(gi_pd_fitness_gene.loc[:,"fold_change"].mean(),0,600,color="red",linestyles="dashed")
# plt.vlines(gi_pd_fitness_gene.loc[:,"fold_change"].mean()+gi_pd_fitness_gene.loc[:,"fold_change"].std(),0,500,color="gray",linestyles="dashed")
# plt.vlines(gi_pd_fitness_gene.loc[:,"fold_change"].mean()-gi_pd_fitness_gene.loc[:,"fold_change"].std(),0,500,color="gray",linestyles="dashed")
plt.xlabel("fold_change")
plt.xticks([-x_m-x_std,0,x_m+x_std],labels=["-$\mu-\sigma$","0","$\mu+\sigma$"],rotation=90);

plt.ylabel("# of genes")

plt.tight_layout()
#plt.savefig("../figures/fold_change_distribution_NRP1_library.png",dpi=300)
# -

x_m+x_std,-x_m-x_std

colors=sns.diverging_palette(145, 300, s=60, as_cmap=True)
genes=["HOM3","BAT1","PIL1","RDI1","STE3","WHI5","RAX1","BEM3","STE11","FUS3","FAR1","BUD6"]
row_colors = ["darkblue" if gi_pd_fitness_gene.loc[genes[i],"p_statistic"] < 0.1 else "gray"  for i in range (0,len(genes))]
g=sns.clustermap(data=gi_pd_fitness_gene.loc[genes,["e_a","e_b"]],cmap=colors,figsize=(5,5),vmin=-1.5,vmax=1.5,metric="euclidean",method="ward",
                 row_cluster=True,col_cluster=False,xticklabels=True,yticklabels=True,row_colors=row_colors)

nrp1_PI_sig,nrp1_NI_sig,nrp1_PI,nrp1_NI,nrp1_PI_all,nrp1_NI_all=classify_GI(gi_pd_fitness_gene,col="fold_change")

len(nrp1_PI_sig),len(nrp1_NI_sig),len(nrp1_PI),len(nrp1_NI),len(nrp1_PI_all),len(nrp1_NI_all)

for i in gi_pd_fitness_gene.index:
    if i in nrp1_NI_sig:
        gi_pd_fitness_gene.loc[i,"significance4FC"]="neg"
    elif i in nrp1_PI_sig:
        gi_pd_fitness_gene.loc[i,"significance4FC"]="pos"
    else:
        gi_pd_fitness_gene.loc[i,"significance4FC"]="none"

gi_pd_fitness_gene[gi_pd_fitness_gene.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)

# +
# look for the fold change of the standard_essentials genes

gi_standard_essentials=defaultdict(dict)
for i in standard_essentials:
    if i in gi_pd_fitness_gene.index:
        gi_standard_essentials[i]["fold_change"]=gi_pd_fitness_gene.loc[i,"fold_change"]
        gi_standard_essentials[i]["p_statistic"]=gi_pd_fitness_gene.loc[i,"p_statistic"]
        gi_standard_essentials[i]["significance"]=gi_pd_fitness_gene.loc[i,"significance"]
        gi_standard_essentials[i]["significance4FC"]=gi_pd_fitness_gene.loc[i,"significance4FC"]

gi_standard_essentials_pd=pd.DataFrame.from_dict(gi_standard_essentials,orient="index")
gi_standard_essentials_pd[gi_standard_essentials_pd.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)

# +
## interaction of NRp1 with the orthologous in yeast of cut7(kinesin) from fission yeast . It was found
## that nrp1 deletion rescue cut7 mutants in pombe.

cu7_orthologous_1="KIP1"
cut7_orthologous_2="CIN8"

gi_pd_fitness_gene.loc[cu7_orthologous_1,"fold_change"],gi_pd_fitness_gene.loc[cut7_orthologous_2,"fold_change"]
# it seems in yeast that the orthologous of cut7 are not interacting with nrp1 therefore nrp1 is not involved in the mitotic spendly assembly
# -

## Interaction with predicted interactors of nrp1,from the profile similarity network from cell map from the significant nrp1 satay interactors
gene_1="CDC28" # pos
gene_2="SEC6" # pos
gene_3="SEC3" # neg
gene_4="SEC5" # neg
gene_5="EXO70" # pos
gene_6="SEC10" # neg
gene_7="SEC15" # npos
gene_8="EXO84" # neg
gene_9="CHS5"
gene_10="STE2"
gene_11="STE3"
gene_12="STE6"
gene_13="CHS7"
gene_14="KIP3"
gene_15="KAR3"
gene_16="CIN8"
for i in [gene_1,gene_2,gene_3,gene_4,gene_5,gene_6,gene_7,gene_8,gene_9,gene_10,gene_11,gene_12,gene_13,gene_14,gene_15,gene_16]:
    if i in gi_pd_fitness_gene.index:
        print(i,gi_pd_fitness_gene.loc[i,"fold_change"],gi_pd_fitness_gene.loc[i,"significance"])
    else:
        print(i,"not found")


# +
from annotate_volcano import annotate_volcano

volcano_df=gi_pd_fitness_gene

fig=annotate_volcano(volcano_df,figure_title="Interactors of nrp1 in WT")

plt.savefig("../figures/fig_volcano_interactors_nrp1.png",dpi=300,transparent=True)

# +
from annotate_volcano import annotate_volcano   #import annotate_volcano function
volcano_df=gi_pd_fitness_gene
trackgene_list=["SPC29","SUI2","PBY1","ECM25"]
fig=annotate_volcano(volcano_df,figure_title="Interactors of nrp1 in WT",trackgene_list=trackgene_list)

plt.savefig("../figures/fig_volcano_interactors_nrp1_with annotations.png",dpi=300,transparent=True)


# +
## interaction score for polarity genes 

polarity_genes_gi=defaultdict(dict)

for i in polarity_genes.index:
    if i in gi_pd_fitness_gene.index:
        polarity_genes_gi[i]["fold_change"]=gi_pd_fitness_gene.loc[i,"fold_change"]
        polarity_genes_gi[i]["p_statistic"]=gi_pd_fitness_gene.loc[i,"p_statistic"]
        polarity_genes_gi[i]["significance"]=gi_pd_fitness_gene.loc[i,"significance"]
    
polarity_genes_gi_pd=pd.DataFrame.from_dict(polarity_genes_gi,orient="index")




# +
## Plot a heatmap with the GI scores for the polarity_genes

import seaborn as sns

x=polarity_genes_gi_pd.loc[:,["fold_change"]]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,10))
# array2heatmap[~np.isfinite(array2heatmap)] = 0
sns.heatmap(x.sort_values(by="fold_change",ascending=False),cmap="PRGn_r",yticklabels=x.sort_values(by="fold_change",ascending=False).index,
cbar=True,annot=True)


# plt.tight_layout()
plt.savefig("../figures/fig_heatmap_polarity_genes_GI_NRP1.png",dpi=300,transparent=True)

# -

polarity_genes.loc[x.sort_values(by="fold_change",ascending=True).index,:]

polarity_genes_gi_pd.loc[x.sort_values(by="fold_change",ascending=True).index]

# +
plt.figure(figsize=(8,6))

plt.scatter(x.sort_values(by="fold_change",ascending=True).index,x.sort_values(by="fold_change",ascending=True),
c=polarity_genes_gi_pd.loc[x.sort_values(by="fold_change",ascending=True).index,"p_statistic"],s=100,cmap="Blues_r",vmax=1,vmin=0)

#plt.hlines(0,xmin=0,xmax=len(x.sort_values(by="fold_change",ascending=True).index),linestyles="dashed",color="black")
plt.grid(linewidth=0.2,linestyle="dashed")
plt.xticks(rotation=90)
plt.ylabel("GI score")
plt.colorbar(label="p-value")
plt.tight_layout()
plt.savefig("../figures/fig_polarity_genes_GI_NRP1.png",dpi=300,transparent=True)

# +
## export to a txt the gene names of the significant positive and negative regulators of nrp1
nrp1_PI_50=gi_pd_fitness_gene.loc[nrp1_PI,:].sort_values(by="fold_change",ascending=False).iloc[:80,:].index
nrp1_NI_50=gi_pd_fitness_gene.loc[nrp1_NI,:].sort_values(by="fold_change",ascending=True).iloc[:80,:].index
with open("../postprocessed-data/NRP1_positive_satay_genes.txt","w") as f:
    for i in nrp1_PI_50:
        f.write(i+"\n")
f.close()
with open("../postprocessed-data/NRP1_negative_satay_genes.txt","w") as f:
    for i in nrp1_NI_50:
        f.write(i+"\n")
f.close()



# +
import gseapy as gp
from gseapy import barplot, dotplot
type_gi="Significant_positive_GI"
goi=nrp1_PI_sig.tolist()
yeast = gp.get_library_name(organism='Yeast')

sets=['GO_Biological_Process_AutoRIF','GO_Cellular_Component_AutoRIF',
'GO_Molecular_Function_AutoRIF','Pfam_Domains_2019','Phenotype_AutoRIF','WikiPathways_2018' ] 

# #%% enrichment 

enr_data=[]

for i in np.arange(0,len(sets)): 

  enr=gp.enrichr(gene_list=goi,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                
                  outdir='../postprocessed-data/enrich-analysis_dnrp1/'+type_gi+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
  enr_data.append(enr.res2d)
      
# to save your figure, make sure that ``ofname`` is not None
  #ax = dotplot(enr.res2d, title=sets[i],cmap='viridis_r', size=20, figsize=(3,5),ofname=sets[i]+type_gi)
  
    

# +
# #%% barplot

# #%% barplot

fig,ax=plt.subplots(1,2,figsize=(25, 10))

plt.subplots_adjust(wspace=1.2, hspace=0.5)
sets_reset=[0,5]


for i,j in zip(sets_reset,np.arange(0,len(sets_reset))):
    data=enr_data[i]
    data=data[data.loc[:,"Adjusted P-value"]!=0]
    tmp=data.sort_values(by="Adjusted P-value",ascending=True)[0:10]
    
    tmp_axis_x=[]
    for k in np.arange(0,len(tmp)):
        tmp_axis_x.append(tmp.Term.tolist()[k].replace("_"," "))

    ax[j].barh(tmp_axis_x,tmp.loc[:,"Adjusted P-value"],color="purple",alpha=0.5)
    ax[j].set_title(sets[i].replace("_"," "),fontsize=18)
    ax[j].set_xlabel("Adjusted P-value",fontsize=18)
    ax[j].tick_params(axis='both', labelsize=18)


    # ax[j//2,j%2].set_xlim(0,0.025)
    
plt.tight_layout()

plt.savefig("../figures/fig_barplot_enrichment_significant_positive_GI_NRP1.png",dpi=300,transparent=True)

# +
import gseapy as gp
from gseapy import barplot, dotplot
type_gi="Significant_negative_GI"
goi=nrp1_NI_sig.tolist()
yeast = gp.get_library_name(organism='Yeast')



# #%% enrichment 
enr_data=[]
for i in np.arange(0,len(sets)): 

=======
plt.grid(linewidth=0.5)
plt.xlabel("fitness")
plt.ylabel("relative intensity dbem1")
plt.tight_layout()
plt.savefig("../figures/fitness_vs_expression_dnrp1dbem1_dbem1_URA.png",dpi=300)

# +
plt.figure(figsize=(5,5))
plt.scatter(fitness_domains.loc[over_index],data2plot.loc[over_index,"r_dbem1dnrp1_vs_dbem1"],
c="purple",alpha=0.4,s=50,vmin=0,vmax=1.2)

plt.scatter(fitness_domains.loc[under_index],data2plot.loc[under_index,"r_dbem1dnrp1_vs_dbem1"],
c="green",alpha=0.4,s=50,vmin=0,vmax=1.2)

plt.scatter(fitness_domains.loc[equal_index],data2plot.loc[equal_index,"r_dbem1dnrp1_vs_dbem1"],
c="grey",alpha=0.4,s=100,vmin=0,vmax=1.2)
plt.hlines(1,0,2,linestyles="dashed",color="k")

plt.grid(linewidth=0.5)
plt.xlabel("fitness")
plt.ylabel("relative intensity dbem1")
plt.tight_layout()

# +
over_mean=fitness_average.loc[over_index].mean()
over_std=fitness_average.loc[over_index].std()/np.sqrt(len(over_index))

under_mean=fitness_average.loc[under_index].mean()
under_std=fitness_average.loc[under_index].std()/np.sqrt(len(under_index))

print(over_mean,over_std,under_mean,under_std)

# +
over_mean=fitness_domains.loc[over_index].mean()
over_std=fitness_domains.loc[over_index].std()/np.sqrt(len(over_index))

under_mean=fitness_domains.loc[under_index].mean()
under_std=fitness_domains.loc[under_index].std()/np.sqrt(len(under_index))

print(over_mean,over_std,under_mean,under_std)
# -

x=fitness_domains.loc[under_index]
x[x>2]

y=fitness_domains.loc[over_index]
y[y>2]


# +
x=list_data_pd.loc["wt_merged"]
y=list_data_pd.loc["bem1-aid_merged"]
z=list_data_pd.loc["dnrp1_merged"]

nrp1_insertions=x[x.loc[:,"Gene name"]=="NRP1"].loc[:,"Insertions"]
nrp1_insertions_dbem1=y[y.loc[:,"Gene name"]=="NRP1"].loc[:,"Insertions"]

bem1_insertions=x[x.loc[:,"Gene name"]=="BEM1"].loc[:,"Insertions"]
bem1_insertions_dnrp1=z[z.loc[:,"Gene name"]=="BEM1"].loc[:,"Insertions"]


data=(data_fitness.loc["wt_merged"])
data_average=(data.loc[:,"fitness_gene"])
data_norm=(data_average-data_average.min())/(data_average.max()-data_average.min())

mean_bar=data_average.loc["NRP1"]
std_error=data.loc["NRP1","fitness_gene_std"]/np.sqrt(nrp1_insertions) # merging of two samples that contain in total 56 insertions in nrp1,rendering 56 different replicates

data_y=(data_fitness.loc["bem1-aid_merged"])
data_y_average=(data_y.loc[:,"fitness_gene"])
data_y_norm=(data_y_average-data_y_average.min())/(data_y_average.max()-data_y_average.min())

mean_bar_y=data_y_average.loc["NRP1"]
std_error_y=data_y.loc["NRP1","fitness_gene_std"]/np.sqrt(nrp1_insertions_dbem1) # merging of two samples that contain in total 56 insertions in nrp1,rendering 56 different replicates

mean_bar_bem1=data_average.loc["BEM1"]
std_error_bem1=data.loc["BEM1","fitness_gene_std"]/np.sqrt(bem1_insertions) 

data_z=data_fitness.loc["dnrp1_merged"]
mean_bar_z=data_z.loc["BEM1","fitness_gene"]
std_error_z=data_z.loc["BEM1","fitness_gene_std"]/np.sqrt(bem1_insertions_dnrp1)






# +
data_a=data_fitness.loc["wt_a"]
data_b=data_fitness.loc["wt_b"]

data_a.loc["BEM1"],data_b.loc["BEM1"]

# +
# select the the values of nrp1 fitness that are below 0.5

nrp1_essentials=[]

for i in data_z.index:
    if data_z.loc[i,"fitness_gene"]<0.5:
        nrp1_essentials.append(i)




# +
## Plot fitness of nrp1 taking the fitness average  as an error bar plot



## plotting in a error bar plot

plt.figure(figsize=(8,5))

plt.bar(["fitness nrp1 in WT","fitness bem1 in WT","fitness nrp1 in dbem1"],
[mean_bar,mean_bar_bem1,mean_bar_y],
yerr=[std_error.values[0],std_error_bem1.values[0],std_error_y.values[0]],
color=["#E8A6BF70","#1EC9E0","#E8A6BF70"],capsize=5,alpha=0.5)

# write every mean value on top of every background
# plt.text(0,mean_bar+0.5,"{:.2f}".format(mean_bar),ha="center",va="center",color="black")
# plt.text(1,mean_bar_bem1+0.5,"{:.2f}".format(mean_bar_bem1),ha="center",va="center",color="black")
# plt.text(2,mean_bar_y+0.5,"{:.2f}".format(mean_bar_y),ha="center",va="center",color="black")



plt.ylabel("Fitness from SATAY")
plt.grid(linewidth=0.2)
plt.ylabel("Fitness normalized")
# plt.yticks(np.arange(0,30,step=2));

plt.tight_layout()

plt.savefig("../figures/fitness_dnrp1_dbem1.png",dpi=300)

# +
## Plot distributions of fitness landscapes in WT and dbem1 

plt.figure(figsize=(8,5))

## Normalize fitness .values()


# plt.hist(data_y_norm/data_y_norm.median(),bins=20,color="#1EC9E0",alpha=0.5,label="$\Delta$bem1");
# plt.hist(data_norm/data_norm.median(),bins=20,color="gray",alpha=0.8,label="WT");
plt.boxplot([data_y_average,data_average],labels=["$\Delta$bem1","WT"],
showfliers=True,patch_artist=False)

plt.yticks(np.arange(0,60,step=10));

plt.ylabel("Fitness normalized")
plt.grid(linewidth=0.2)

plt.tight_layout()
plt.savefig("../figures/fitness_distributions_dbem1_wt.png",dpi=300)
# -

x=data_y_norm/data_y_norm.median()
x.max()

data_y_norm.loc["NRP1"]/data_y_norm.median(),data_norm.loc["NRP1"]/data_norm.median(),data_norm.loc["BEM1"]/data_norm.median()

# +
# fitness of cla4 and whi3 in WT and dnrp1

# fitness_cla4=data.loc["CLA4","fitness_gene"]
# fitness_cla4_nrp1=data_z.loc["CLA4","fitness_gene"]

# fitness_whi3=data.loc["WHI3","fitness_gene"]
# fitness_whi3_nrp1=data_z.loc["WHI3","fitness_gene"]

# fitness_cla4_std=data.loc["CLA4","fitness_gene_std"]/np.sqrt(x[x.loc[:,"Gene name"]=="CLA4"].loc[:,"Insertions"])
# fitness_cla4_nrp1_std=data_z.loc["CLA4","fitness_gene_std"]/np.sqrt(z[z.loc[:,"Gene name"]=="CLA4"].loc[:,"Insertions"])

# fitness_whi3_std=data.loc["WHI3","fitness_gene_std"]/np.sqrt(x[x.loc[:,"Gene name"]=="WHI3"].loc[:,"Insertions"])
# fitness_whi3_nrp1_std=data_z.loc["WHI3","fitness_gene_std"]/np.sqrt(z[z.loc[:,"Gene name"]=="WHI3"].loc[:,"Insertions"])

# +
# ## bar plot 

# plt.figure(figsize=(12,5))

# plt.bar(["fitness cla4 in WT","fitness cla4 in dnrp1","fitness whi3 in WT","fitness whi3 in dnrp1"],
# [fitness_cla4,fitness_cla4_nrp1,fitness_whi3,fitness_whi3_nrp1],
# yerr=[fitness_cla4_std.values[0],fitness_cla4_nrp1_std.values[0],fitness_whi3_std.values[0],fitness_whi3_nrp1_std.values[0]],
# color="black",capsize=5,alpha=0.5)

# +
## Plot insertion map onto gene of interest


background="wt_merged"
data=list_data_pd.loc[background]
## remove ade2 data from the dataset
# data.drop(labels=5805,axis=0,inplace=True)
# ## remove the ura3 data from the dataset
# data.drop(labels=1927,axis=0,inplace=True)


## Extract as lists the gene coordinates, and the locations of every insertion and reads per gene 
reads_locations=[]
insertion_locations=[]
gene_coordinates=[]

for j in data.index:
        coi=from_excel_to_list(data.loc[j]["Reads per insertion location"])
        coi_insertions=from_excel_to_list(data.loc[j]["Insertion locations"])
        gene_coordinates.append([data.loc[j,"Start location"],data.loc[j,"End location"]])
        reads_locations.append(coi)
        insertion_locations.append(coi_insertions)
        


# compute the reads per insertions for each gene along the gene length

gene_parts=np.linspace(0,1,11)
r=np.zeros(shape=(len(insertion_locations),len(gene_parts))) # reads_per_insertion_parts array , every gene in the rows 

for i in np.arange(0,len(insertion_locations)):
    if (insertion_locations[i])!=0:
        g=np.array(insertion_locations[i]) # insertion locations
        f=np.linspace(gene_coordinates[i][0],gene_coordinates[i][1],len(gene_parts)) # gene coordinates
        #f=np.array(gene_coordinates[i][0]+gene_coordinates[i][1]*gene_parts)
        binedges = g.searchsorted(f)

        rngs = [list(range(binedges[k], binedges[k+1])) for k in range(len(binedges)-1)]

        for k in np.arange(0,len(rngs)):
            readsperinsert=[]
            for j in np.arange(0,len(rngs[k])):
                readsperinsert.append(reads_locations[i][rngs[k][j]])
                
            if len(readsperinsert)>1:
                r[i,k]=np.sum(readsperinsert)/(len(readsperinsert)-1)#discarding the insertion with the highest read count
            else:
                r[i,k]=0

    



# +
nrp1_wt_index=np.where(data.loc[:,"Gene name"]=="NRP1")[0]
nrp1_coordinates=np.linspace(gene_coordinates[nrp1_wt_index[0]][0],gene_coordinates[nrp1_wt_index[0]][1],11)
nrp1_coordinates2text=["{:.0f}".format(i) for i in nrp1_coordinates]



# +
plt.figure(figsize=(10,5))

plt.bar(gene_parts,r[nrp1_wt_index[0],:],width=0.1)
plt.xticks(np.linspace(-0.05,0.95,11),
labels=nrp1_coordinates2text,rotation=90);
# -

data_fitness.loc["bem1-aid_merged","NRP1"]

bem1d=list_data_pd.loc["bem1-aid_merged"]
bem1d.index=bem1d["Gene name"]
bem1d

# +


reads_nrp1_in_dbem1=from_excel_to_list(bem1d.loc["NRP1","Reads per insertion location"])

plt.bar(np.arange(0,len(reads_nrp1_in_dbem1)),reads_nrp1_in_dbem1)

# +
## Computing the interactors of NRP1 using the whole fitness of the gene 
gene="NRP1"
nrp1_f_a=data_fitness.loc["wt_merged",gene]["fitness_gene"] # in wt_A , nrp1 is not found
nrp1_f_b=data_fitness.loc["wt_b",gene]["fitness_gene"]
data_b=data_fitness.loc["wt_b"]
data_a=data_fitness.loc["wt_a"]

intersection_genes=list((set(data_b.index)&set(data_a.index)))

significance_threshold = 0.05 #set significance threshold
gi=defaultdict(dict)
ttest_tval_list = [np.nan]*2 #initialize list for storing t statistics
ttest_pval_list = [np.nan]*2 #initialize list for storing p-values
signif_thres_list = False #initialize boolean list for indicating datapoints with p-value above threshold
fc_list = [np.nan]*2
for gene in intersection_genes :
    geneX=gene
    
    geneX_f_a=data_fitness.loc["wt_a",geneX]["fitness_gene"]
    if geneX in data_fitness.loc["wt_b"].index:
        geneX_f_b=data_fitness.loc["wt_b",geneX]["fitness_gene"]
        if geneX in data_fitness.loc["dnrp1_1"].index and geneX in data_fitness.loc["dnrp1_2"].index:
            geneXnrp1_f_a=data_fitness.loc["dnrp1_1",geneX]["fitness_gene"]
            geneXnrp1_f_b=data_fitness.loc["dnrp1_2",geneX]["fitness_gene"]
            
    variable_a_array=[geneXnrp1_f_a,geneXnrp1_f_b]
    variable_b_array=[geneX_f_a*nrp1_f_a,geneX_f_b*nrp1_f_b]
    
    ttest_val = stats.ttest_ind(variable_a_array, variable_b_array) #T-test
    gi[gene]["p_statistic"]=ttest_val[1]
    ttest_tval_list = ttest_val[0]
    gi[gene]["t_statistic"]=ttest_tval_list
    if not ttest_val[1] == 0: #prevent p=0 to be inputted in log
        ttest_pval_list = -1*np.log10(ttest_val[1])
        
    else:
        ttest_pval_list = 0
    gi[gene]["p_value"]=ttest_pval_list
    if ttest_pval_list >= -1*np.log10(significance_threshold):
        gi[gene]["significance"]=True
    else:
        gi[gene]["significance"]=False
    #DETERMINE FOLD CHANGE PER GENE
    
    fc_list=np.mean(variable_a_array)-np.mean(variable_b_array)
    gi[gene]["fold_change"]=fc_list

        


# +
gi_pd_fitness_gene=pd.DataFrame.from_dict(gi,orient="index")

gi_pd_fitness_gene["gene_names"]=gi_pd_fitness_gene.index

gi_pd_fitness_gene

# +
## normalize the fold change to go from -1 to 1

gi_pd_fitness_gene["fold_change_norm"]=(gi_pd_fitness_gene["fold_change"]-gi_pd_fitness_gene["fold_change"].min())/(gi_pd_fitness_gene["fold_change"].max()-gi_pd_fitness_gene["fold_change"].min())*2-1


## add significance based on the log2FC
gi_pd_significant=gi_pd_fitness_gene[gi_pd_fitness_gene.loc[:,"significance"]==True]

# Finding the std of the fold change of the significant genes, genes whose interaction score is more than 1 std away from the mean are considered as satay genes

std_fc=gi_pd_significant.loc[:,"fold_change_norm"].std()
mean_fc=gi_pd_significant.loc[:,"fold_change_norm"].mean()

neg_satay_signif=gi_pd_significant[gi_pd_significant.loc[:,"fold_change_norm"]<mean_fc-0.5*std_fc].loc[:,"gene_names"].index

pos_satay_signif=gi_pd_significant[gi_pd_significant.loc[:,"fold_change_norm"]>mean_fc+std_fc].loc[:,"gene_names"].index

for i in gi_pd_fitness_gene.index:
    if i in neg_satay_signif:
        gi_pd_fitness_gene.loc[i,"significance4FC"]="neg"
    elif i in pos_satay_signif:
        gi_pd_fitness_gene.loc[i,"significance4FC"]="pos"
    else:
        gi_pd_fitness_gene.loc[i,"significance4FC"]="none"
# -

gi_pd_fitness_gene.loc["EXO70"]

gi_pd_fitness_gene[gi_pd_fitness_gene.loc[:,"significance"]==True].sort_values(by="fold_change_norm",ascending=False)

# +
# look for the fold change of the standard_essentials genes

gi_standard_essentials=defaultdict(dict)
for i in standard_essentials:
    if i in gi_pd_fitness_gene.index:
        gi_standard_essentials[i]["fold_change"]=gi_pd_fitness_gene.loc[i,"fold_change"]
        gi_standard_essentials[i]["p_statistic"]=gi_pd_fitness_gene.loc[i,"p_statistic"]
        gi_standard_essentials[i]["significance"]=gi_pd_fitness_gene.loc[i,"significance"]
        gi_standard_essentials[i]["fold_change_norm"]=gi_pd_fitness_gene.loc[i,"fold_change_norm"]
        gi_standard_essentials[i]["significance4FC"]=gi_pd_fitness_gene.loc[i,"significance4FC"]

gi_standard_essentials_pd=pd.DataFrame.from_dict(gi_standard_essentials,orient="index")
gi_standard_essentials_pd[gi_standard_essentials_pd.loc[:,"significance"]==True].sort_values(by="fold_change_norm",ascending=False)

# +
## interaction of NRp1 with the orthologous in yeast of cut7(kinesin) from fission yeast . It was found
## that nrp1 deletion rescue cut7 mutants in pombe.

cu7_orthologous_1="KIP1"
cut7_orthologous_2="CIN8"

gi_pd_fitness_gene.loc[cu7_orthologous_1,"fold_change_norm"],gi_pd_fitness_gene.loc[cut7_orthologous_2,"fold_change_norm"]
# it seems in yeast that the orthologous of cut7 are not interacting with nrp1 therefore nrp1 is not involved in the mitotic spendly assembly
# -

## Interaction with predicted interactors of nrp1,from the profile similarity network from cell map from the significant nrp1 satay interactors
gene_1="SEC8" # pos
gene_2="DRS2" # pos
gene_3="SEC3" # neg
gene_4="SHE4" # neg
gene_5="EPO1" # pos
gene_6="VAB2" # neg
gene_7="SOM1" # npos
gene_8="SCC2" # neg
gene_9="WHI3"
gene_10="CLN3"
gene_11="WHI1"
gene_12="KAR3"
gene_13="BEM1"
for i in [gene_1,gene_2,gene_3,gene_4,gene_5,gene_6,gene_7,gene_8,gene_9,gene_10,gene_11,gene_12,gene_13]:
    if i in gi_pd_fitness_gene.index:
        print(i,gi_pd_fitness_gene.loc[i,"fold_change_norm"],gi_pd_fitness_gene.loc[i,"significance"])
    else:
        print(i,"not found")


# +
from annotate_volcano import annotate_volcano

volcano_df=gi_pd_fitness_gene

fig=annotate_volcano(volcano_df,figure_title="Interactors of nrp1 in WT")

plt.savefig("../figures/fig_volcano_interactors_nrp1.png",dpi=300,transparent=True)

# +
from annotate_volcano import annotate_volcano   #import annotate_volcano function
volcano_df=gi_pd_fitness_gene
trackgene_list=["COX9","NTF2","PBY1","ECM25"]
fig=annotate_volcano(volcano_df,figure_title="Interactors of nrp1 in WT",trackgene_list=trackgene_list)

plt.savefig("../figures/fig_volcano_interactors_nrp1_with annotations.png",dpi=300,transparent=True)


# +
## Make a histogram of the fold change of the interactors of NRP1

plt.figure(figsize=(5,5))


plt.hist(gi_pd_fitness_gene.loc[gi_pd_fitness_gene.loc[:,"significance"]==True,"fold_change_norm"],bins=20,
color="#E8A6BF70",alpha=0.5);

plt.xlabel("Normalized genetic interaction score")

plt.ylabel("Number of genes")

plt.grid(linewidth=0.2)

plt.tight_layout()

plt.savefig("../figures/fig_histogram_siginificant_interactors_nrp1.png",dpi=300,transparent=True)

# +
## interaction score for polarity genes 

polarity_genes_gi=defaultdict(dict)

for i in polarity_genes.index:
    if i in gi_pd_fitness_gene.index:
        polarity_genes_gi[i]["fold_change_norm"]=gi_pd_fitness_gene.loc[i,"fold_change_norm"]
        polarity_genes_gi[i]["p_statistic"]=gi_pd_fitness_gene.loc[i,"p_statistic"]
        polarity_genes_gi[i]["significance"]=gi_pd_fitness_gene.loc[i,"significance"]
    
polarity_genes_gi_pd=pd.DataFrame.from_dict(polarity_genes_gi,orient="index")

len(polarity_genes_gi_pd)


# +
x=polarity_genes_gi_pd[polarity_genes_gi_pd.loc[:,"significance"]==True]

## elete nrp1 from the rows

x=x.drop("NRP1")
# -

polarity_genes.loc["STE18"]

# +
## Plot a heatmap with the GI scores for the polarity_genes

import seaborn as sns

x=x.loc[:,["fold_change_norm"]]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,5))
# array2heatmap[~np.isfinite(array2heatmap)] = 0
sns.heatmap(x.sort_values(by="fold_change_norm",ascending=False),cmap="PRGn_r",vmin=-1,vmax=1,yticklabels=x.sort_values(by="fold_change_norm",ascending=False).index,
cbar=True,annot=True)


# plt.tight_layout()
# plt.savefig("../figures/fig_heatmap_polarity_genes_GI_NRP1.png",dpi=300,transparent=True)


# +
gi_pd=gi_pd_fitness_gene
gi_pd_significant=gi_pd[gi_pd.loc[:,"significance"]==True]

# Finding the std of the fold change of the significant genes, genes whose interaction score is more than 1 std away from the mean are considered as satay genes

std_fc=gi_pd_significant.loc[:,"fold_change_norm"].std()
mean_fc=gi_pd_significant.loc[:,"fold_change_norm"].mean()

neg_satay_signif=gi_pd_significant[gi_pd_significant.loc[:,"fold_change_norm"]<mean_fc-0.5*std_fc].loc[:,"gene_names"].index

neg_satay=gi_pd[gi_pd.loc[:,"fold_change_norm"]<mean_fc-0.5*std_fc].loc[:,"gene_names"].index # all genes 

neg_satay_all=gi_pd[gi_pd.loc[:,"fold_change_norm"]<0].loc[:,"gene_names"].index # all genes

pos_satay_signif=gi_pd_significant[gi_pd_significant.loc[:,"fold_change_norm"]>mean_fc+1*std_fc].loc[:,"gene_names"].index
pos_satay=gi_pd[gi_pd.loc[:,"fold_change_norm"]>1*std_fc].loc[:,"gene_names"].index # all genes 
pos_satay_all=gi_pd[gi_pd.loc[:,"fold_change_norm"]>0].loc[:,"gene_names"].index # all genes

# -

mean_fc-1*std_fc,mean_fc+1*std_fc,len(neg_satay_signif),len(pos_satay_signif)

# +
## export to a txt the gene names of the significant positive and negative regulators of nrp1

with open("../postprocessed-data/positive_satay_genes.txt","w") as f:
    for i in pos_satay_signif:
        f.write(i+"\n")
f.close()
with open("../postprocessed-data/negative_satay_genes.txt","w") as f:
    for i in neg_satay_signif:
        f.write(i+"\n")
f.close()

## export to a txt the 30 more extreme GI for cellmap subnetwork 
pos_satay_signif_50=gi_pd_significant.loc[pos_satay_signif,:].sort_values(by="fold_change_norm",ascending=False).iloc[:50].index
neg_satay_signif_50=gi_pd_significant.loc[neg_satay_signif,:].sort_values(by="fold_change_norm",ascending=True).iloc[:50].index

with open("../postprocessed-data/positive_satay_genes_30.txt","w") as f:
    for i in pos_satay_signif_50:
        f.write(i+"\n")
f.close()
with open("../postprocessed-data/negative_satay_genes_30.txt","w") as f:
    for i in neg_satay_signif_50:
        f.write(i+"\n")
f.close()

# +
import gseapy as gp
from gseapy import barplot, dotplot
type_gi="Significant_positive_GI"
goi=pos_satay_signif.tolist()
yeast = gp.get_library_name(organism='Yeast')

sets=['GO_Biological_Process_AutoRIF','GO_Cellular_Component_AutoRIF',
'GO_Molecular_Function_AutoRIF','Pfam_Domains_2019','Phenotype_AutoRIF' ] 

# #%% enrichment 

enr_data_pos=[]

for i in np.arange(0,len(sets)): 

>>>>>>> 1c5527a6c3a1e3cda5126ec83fb46e9a8a798f6f
  enr=gp.enrichr(gene_list=goi,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                
                  outdir='../postprocessed-data/enrich-analysis_dnrp1/'+type_gi+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
<<<<<<< HEAD

  enr_data.append(enr.res2d)
=======
  enr_data_pos.append(enr.res2d)
      
# to save your figure, make sure that ``ofname`` is not None
  #ax = dotplot(enr.res2d, title=sets[i],cmap='viridis_r', size=20, figsize=(3,5),ofname=sets[i]+type_gi)
  
    

# +
# #%% barplot

fig,ax=plt.subplots(2,2,figsize=(25, 10))

plt.subplots_adjust(wspace=1.2, hspace=0.5)
sets_reset=[0,1,2,4]


for i,j in zip(sets_reset,np.arange(0,len(sets_reset))):
    ax[j//2,j%2].barh(enr_data_pos[i].Term[0:3],enr_data_pos[i].loc[0:2,"Adjusted P-value"],color="purple",alpha=0.5)
    ax[j//2,j%2].set_title(sets[i])
    ax[j//2,j%2].set_xlabel("Adjusted P-value")
    ax[j//2,j%2].set_xlim(0,0.025)
    
plt.tight_layout()
plt.savefig("../figures/fig_barplot_enrichment_significant_positive_GI_NRP1.png",dpi=300,transparent=True)

# +
import gseapy as gp
from gseapy import barplot, dotplot
type_gi="Significant_negative_GI"
goi=neg_satay_signif.tolist()
yeast = gp.get_library_name(organism='Yeast')

sets=['GO_Biological_Process_AutoRIF','GO_Cellular_Component_AutoRIF',
'GO_Molecular_Function_AutoRIF','Pfam_Domains_2019','Phenotype_AutoRIF' ] 

# #%% enrichment 
enr_data_neg=[]
for i in np.arange(0,len(sets)): 

  enr=gp.enrichr(gene_list=goi,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                
                  outdir='../postprocessed-data/enrich-analysis_dnrp1/'+type_gi+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                )

  enr_data_neg.append(enr.res2d)
>>>>>>> 1c5527a6c3a1e3cda5126ec83fb46e9a8a798f6f
      
# to save your figure, make sure that ``ofname`` is not None
  #ax = dotplot(enr.res2d, title=sets[i],cmap='viridis_r', size=20, figsize=(3,5),ofname=sets[i]+type_gi)
  

    
<<<<<<< HEAD

# +
# #%% barplot

fig,ax=plt.subplots(1,2,figsize=(25, 10))

plt.subplots_adjust(wspace=1.2, hspace=0.5)
sets_reset=[0,5]


for i,j in zip(sets_reset,np.arange(0,len(sets_reset))):
    data=enr_data[i]
    data=data[data.loc[:,"Adjusted P-value"]!=0]
    tmp=data.sort_values(by="Adjusted P-value",ascending=True)[0:10]
    tmp_axis_x=[]
    for k in np.arange(0,len(tmp)):
        tmp_axis_x.append(tmp.Term.tolist()[k].replace("_"," "))

    ax[j].barh(tmp_axis_x,tmp.loc[:,"Adjusted P-value"],color="green",alpha=0.5)
    ax[j].set_title(sets[i].replace("_"," "),fontsize=18)
    ax[j].set_xlabel("Adjusted P-value",fontsize=18)
    ax[j].tick_params(axis='both', labelsize=18)
    # ax[j].set_xscale("log")


    # ax[j//2,j%2].set_xlim(0,0.025)
    
plt.tight_layout()
plt.savefig("../figures/fig_barplot_enrichment_significant_negative_GI_NRP1.png",dpi=300,transparent=True)

# +
import gseapy as gp
from gseapy import barplot, dotplot
type_gi="Negative_mass-spec_GI"
#goi=["RDI1","STE3","WHI5","RAX1"]
#goi=["STE11", "BEM3","STE2","FUS3","FAR1"]
goi=["HOM3","BAT1"]
yeast = gp.get_library_name(organism='Yeast')

sets=['GO_Biological_Process_AutoRIF','GO_Cellular_Component_AutoRIF',
'GO_Molecular_Function_AutoRIF','Pfam_Domains_2019','Phenotype_AutoRIF' ] 

# #%% enrichment 
enr_data=[]
for i in np.arange(0,len(sets)): 

  enr=gp.enrichr(gene_list=goi,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                
                  outdir='../postprocessed-data/enrich-analysis_dnrp1/'+type_gi+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                )

  enr_data.append(enr.res2d)
# -

enr_data[0]

# +
neg_nrp1=pd.read_excel("../postprocessed-data/NRP1_genetic_interactions_filtered_by_negat.xlsx",skiprows=8)
pos_nrp1=pd.read_excel("../postprocessed-data/NRP1_genetic_interactions_filtered_by_pos.xlsx",skiprows=8)

e_map_gi_wilmes=pd.read_excel("../postprocessed-data/e-map-wilmes-2008-data-GI-nrp1.xlsx",sheet_name=[1,2,3],skiprows=0)
# -

neg_nrp1=neg_nrp1[neg_nrp1.Action=="Hit"]
pos_nrp1=pos_nrp1[pos_nrp1.Action=="Hit"]

## interactors from the emap assay 
for i in np.arange(1,len(e_map_gi_wilmes)+1):

    columns_clean=[]

    for j in np.arange(0,len(e_map_gi_wilmes[i].columns)):
        columns_clean.append(e_map_gi_wilmes[i].columns[j].split(" ")[0])

    e_map_gi_wilmes[i].columns=columns_clean


# +
gi_nrp1_emap=pd.DataFrame()

for i in np.arange(1,len(e_map_gi_wilmes)+1):
    gi_nrp1_emap=pd.concat([gi_nrp1_emap,e_map_gi_wilmes[i][e_map_gi_wilmes[i].loc[:,"Gene"]=="NRP1 - DELETION"]])

gi_nrp1_emap.index=np.arange(0,len(gi_nrp1_emap))
gi_nrp1_emap.drop(columns=["Gene"],inplace=True)

# +
## select positive scores as positive and negative scores as negative

positive_emap_nrp1=[]
negative_emap_nrp1=[]
for i in gi_nrp1_emap.columns:
    for j in  gi_nrp1_emap.index:
        if gi_nrp1_emap.loc[j,i]!=np.nan:
            if gi_nrp1_emap.loc[j,i]>0:
                positive_emap_nrp1.append(i)
            elif gi_nrp1_emap.loc[j,i]<0:
                negative_emap_nrp1.append(i)
        
positive_emap_nrp1=set(positive_emap_nrp1)
negative_emap_nrp1=set(negative_emap_nrp1)
# -

len(positive_emap_nrp1),len(negative_emap_nrp1),positive_emap_nrp1.intersection(negative_emap_nrp1)

neg_nrp1_genes=set(neg_nrp1.loc[:,"Interactor.1"])
pos_nrp1_genes=set(pos_nrp1.loc[:,"Interactor.1"])

len(neg_nrp1_genes),len(pos_nrp1_genes)

# +
## Exclude known interactors that were not analyzed in the study because their fitness in satay is annotated as zero, for example :
neg_check=[]
pos_check=[]
neg_emap_check=[]
pos_emap_check=[]
for i in neg_nrp1_genes :
    if i in gi_pd_fitness_gene.index:
        neg_check.append(i)
for i in pos_nrp1_genes:
    if i  in gi_pd_fitness_gene.index:
        pos_check.append(i)

for i in positive_emap_nrp1:
    if i in gi_pd_fitness_gene.index:
        pos_emap_check.append(i)

for i in negative_emap_nrp1:
    if i in gi_pd_fitness_gene.index:
        neg_emap_check.append(i)

len(neg_check),len(pos_check),len(neg_nrp1_genes),len(pos_nrp1_genes),len(neg_emap_check),len(pos_emap_check)
# -

for i in pos_check:
    if i in standard_essentials:
        print(i)

neg_nrp1[neg_nrp1.loc[:,"Interactor.1"].isin(neg_check)].Reference

pos_nrp1[pos_nrp1.loc[:,"Interactor.1"].isin(pos_check)].Reference

pos_check,neg_check

# +
consistent_neg_interactors_nrp1=list(set(nrp1_NI_sig) & set(neg_check) )

consistent_pos_interactors_nrp1=list(set(nrp1_PI_sig) & set(pos_check))

common_neg_interactors_nrp1=list(set(nrp1_NI_all) & set(neg_check) )

common_pos_interactors_nrp1=list(set(nrp1_PI_all) & set(pos_check))

pos_satay_neg_costanzo=list(set(nrp1_PI_all) & set(neg_check))

neg_satay_pos_costanzo=list(set(nrp1_NI_all) & set(pos_check))



# +
print("Number of consistent negative interactors with NRP1 : ",len(consistent_neg_interactors_nrp1))
print("Number of consistent positive interactors with NRP1 : ",len(consistent_pos_interactors_nrp1))


print("Number of existing negative interactors with NRP1 : ",len(neg_check))

print("Ratio of common negative interactors with NRP1 : ",len(common_neg_interactors_nrp1)/len(neg_check))
print("Common negative interactors gene",common_neg_interactors_nrp1)
print("Number of positive interactors with NRP1 in satay and negative in costanzo : ",len(pos_satay_neg_costanzo))
print("Genes : ",pos_satay_neg_costanzo)

print("Number of existing positive interactors with NRP1 : ",len(pos_check))
print("Ratio of common positive interactors with NRP1 : ",len(common_pos_interactors_nrp1)/len(pos_check))
print("Common positive interactors gene",common_pos_interactors_nrp1)


print("Number of negative interactors with NRP1 in satay and positive in costanzo : ",len(neg_satay_pos_costanzo))
print("Genes : ",neg_satay_pos_costanzo)
=======

# +
# #%% barplot

fig,ax=plt.subplots(2,2,figsize=(18, 10))

plt.subplots_adjust(wspace=1.8, hspace=0.5)
sets_reset=[0,1,2,4]


for i,j in zip(sets_reset,np.arange(0,len(sets_reset))):
    ax[j//2,j%2].barh(enr_data_neg[i].Term[0:3],enr_data_neg[i].loc[0:2,"Adjusted P-value"],color="green",alpha=0.5)
    ax[j//2,j%2].set_title(sets[i])
    ax[j//2,j%2].set_xlabel("Adjusted P-value")
    ax[j//2,j%2].set_xscale("log")

plt.tight_layout()
plt.savefig("../figures/fig_barplot_enrichment_significant_negative_GI_NRP1.png",dpi=300,transparent=True)

# +
## Viz by the ClueGo app from cytoscape  

posGI_go=pd.read_excel("../postprocessed-data/ClueGOResultTable-PG-NRP1.xlsx",engine="openpyxl")

negGI_go=pd.read_excel("../postprocessed-data/ClueGOResultTable-NG-NRP1.xlsx",engine="openpyxl")

posGI_go.head()

# +
## Bar plot of the GO terms for different groups 
pos_groups=posGI_go.GOGroups.unique()
posGI_go.sort_values(by="% Associated Genes",ascending=False,inplace=True)
terms=posGI_go.Term.unique()[0:10]
# genes=posGI_go["% Associated Genes"][0:10]
# pvalue=posGI_go["Term PValue Corrected with Bonferroni step down"][0:10]

genes=[]
pvalue=[]
for i in terms:
    
    genes.append(posGI_go[posGI_go.Term==i]["% Associated Genes"].values[0])
    pvalue.append(posGI_go[posGI_go.Term==i]["Term PValue"].values[0])

fig,ax=plt.subplots(1,1,figsize=(12,8))
ax.barh(terms,genes,color="purple",alpha=0.5)
# change the fontsize
for (i,p) in zip(pvalue,ax.patches):
    if i<0.05 and i>0.01:
        ax.text(p.get_x() + p.get_width()+0.3 , p.get_y()+p.get_height()/ 5., '*', ha='center')
    elif i<0.01 and i>0.001:
        ax.text(p.get_x() + p.get_width()+0.5 , p.get_y()+p.get_height()/ 5., '**', ha='center')
    elif i<0.001:
        ax.text(p.get_x() + p.get_width()+1.4 , p.get_y()+p.get_height()/ 5., '***', ha='center')
ax.tick_params(axis='y', labelsize=18)
ax.set_xlabel("% Associated Genes",fontsize=18)
plt.tight_layout()
plt.savefig("../figures/fig_barplot_enrichment_significant_positive_GI_NRP1_cluego.png",dpi=500,transparent=True)

# +

## Bar plot of the GO terms for different groups 

negGI_go.sort_values(by="% Associated Genes",ascending=False,inplace=True)
terms=negGI_go.Term.unique()[0:10]
# genes=posGI_go["% Associated Genes"][0:10]
# pvalue=posGI_go["Term PValue Corrected with Bonferroni step down"][0:10]

genes=[]
pvalue=[]
for i in terms:
    
    genes.append(negGI_go[negGI_go.Term==i]["% Associated Genes"].values[0])
    pvalue.append(negGI_go[negGI_go.Term==i]["Term PValue"].values[0])

fig,ax=plt.subplots(1,1,figsize=(15,8))
ax.barh(terms,genes,color="green",alpha=0.5)
for (i,p) in zip(pvalue,ax.patches):
    if i<0.05 and i>0.01:
        ax.text(p.get_x() + p.get_width()+0.5 , p.get_y()+p.get_height()/ 5., '*', ha='center')
    elif i<0.01 and i>0.001:
        ax.text(p.get_x() + p.get_width()+0.8 , p.get_y()+p.get_height()/ 5., '**', ha='center')
    elif i<0.001:
        ax.text(p.get_x() + p.get_width()+1.4 , p.get_y()+p.get_height()/ 5., '***', ha='center')
# change the fontsize
ax.tick_params(axis='y', labelsize=18)
ax.set_xlabel("% Associated Genes",fontsize=18)
plt.tight_layout()
plt.savefig("../figures/fig_barplot_enrichment_significant_negative_GI_NRP1_cluego.png",dpi=300,transparent=True)
# -

neg_nrp1=pd.read_excel("../postprocessed-data/NRP1_genetic_interactions_filtered_by_negat.xlsx",skiprows=8)
pos_nrp1=pd.read_excel("../postprocessed-data/NRP1_genetic_interactions_filtered_by_pos.xlsx",skiprows=8)

neg_nrp1_genes=neg_nrp1.loc[:,"Interactor.1"]
pos_nrp1_genes=pos_nrp1.loc[:,"Interactor.1"]

neg_nrp1_genes.unique()

len(neg_nrp1_genes.unique()),len(pos_nrp1_genes.unique())

# +
common_neg_signif=[]
for gene in neg_satay_signif:
    if gene in neg_nrp1_genes.unique():
        common_neg_signif.append(gene)


common_neg=[]

for gene in neg_satay_all:
    if gene in neg_nrp1_genes.unique():
        common_neg.append(gene)

common_pos_signif=[]
for gene in pos_satay_signif:
    if gene in pos_nrp1_genes.unique() :
        common_pos_signif.append(gene)

common_pos=[]

for gene in pos_satay_all:
    if gene in pos_nrp1_genes.unique() :
        common_pos.append(gene)


# +
consistent_neg_interactors_nrp1=list(set(neg_satay_signif) & set(neg_nrp1_genes) )

consistent_pos_interactors_nrp1=list(set(pos_satay_signif) & set(pos_nrp1_genes))

common_neg_interactors_nrp1=list(set(neg_satay_all) & set(neg_nrp1_genes) )

common_pos_interactors_nrp1=list(set(pos_satay_all) & set(pos_nrp1_genes))

common_neg_interactors_nrp1,common_pos_interactors_nrp1,consistent_pos_interactors_nrp1,consistent_neg_interactors_nrp1,

# +
from matplotlib_venn import venn2, venn2_circles

plt.figure(figsize=(10,5))
venn2(subsets = (len(pos_satay_all), len(pos_nrp1_genes), len(common_pos_interactors_nrp1)), 
set_labels = ('SATAy', 'existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(pos_satay_all), len(pos_nrp1_genes), len(common_pos_interactors_nrp1)));
plt.title("Consistency in positive interactions")
#
#plt.savefig("../figures/venn_pos_nrp1.png",dpi=300)
# -

plt.figure(figsize=(10,5))
venn2(subsets = (len(neg_satay_all), len(neg_nrp1_genes), len(common_neg_interactors_nrp1)), 
set_labels = ('SATAy', 'existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(neg_satay_all), len(neg_nrp1_genes), len(common_neg_interactors_nrp1)));
plt.title("Consistency in negative interactions")
#plt.savefig("../figures/venn_neg_nrp1.png",dpi=300)

# +
pos_satay_neg_costanzo=list(set(pos_satay_all) & set(neg_nrp1_genes) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(pos_satay_all), len(neg_nrp1_genes), len(pos_satay_neg_costanzo)), 
set_labels = ('PI SATAy', 'NI existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(pos_satay_all), len(neg_nrp1_genes), len(pos_satay_neg_costanzo)));
plt.title("Misclasification I")

#plt.savefig("../figures/venn_misclass1.png",dpi=300)

# +
neg_satay_pos_costanzo=list(set(neg_satay_all) & set(pos_nrp1_genes) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(neg_satay_all), len(pos_nrp1_genes), len(neg_satay_pos_costanzo)), 
set_labels = ('NI SATAy', 'PI existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(neg_satay_all), len(pos_nrp1_genes), len(neg_satay_pos_costanzo)));
plt.title("Misclasification II")

#plt.savefig("../figures/venn_misclass2.png",dpi=300)

# +
# simple plotting function
from gseapy import barplot, dotplot
for i in np.arange(0,len(sets)):
    
# to save your figure, make sure that ``ofname`` is not None
    ax = dotplot(enr.res2d, title=sets[i],cmap='viridis_r', size=20, figsize=(3,5),ofname=sets[i])
# -

keys= ['wt_merged','dnrp1_merged','dbem3_merged','dbem1dbem3_a']

# ## Volcano plots

# +
path_a = r"../data/"
filelist_a = ["wt_a/WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_pergene.txt",
"wt_b/WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_pergene.txt"]
path_b = r"../data/"
filelist_b = ["dnrp1_1/dnrp1-1_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt",
"dnrp1_2/dnrp1-2_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt"]


variable = 'tn_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.001 #set threshold above which p-values are regarded significant
normalize=True

trackgene_list = ["COX9","NTF2","PBY1","SMD1","SBH1","NRP1","ECM25"] # ["cdc42"]


figure_title = "$\Delta$nrp1 vs WT "

fc_interval=[1.5,-1.5]
pv_values=[3,3]

volcano_df_nrp1_wt = volcano(path_a=path_a, filelist_a=filelist_a,
            path_b=path_b, filelist_b=filelist_b, 
            variable=variable, fold_change_interval=fc_interval,p_value_interval=pv_values,
            significance_threshold=significance_threshold,
            normalize=normalize,
            trackgene_list=trackgene_list,
            figure_title=figure_title)

# +
volcano_df_nrp1_wt.sort_values(by=['fold_change'], inplace=True)

dnrp1_genes_positive_enriched=volcano_df_nrp1_wt[volcano_df_nrp1_wt["significance"]==True][0:50]

dnrp1_genes_negative_enriched=volcano_df_nrp1_wt[volcano_df_nrp1_wt["significance"]==True][-50:]

volcano_df_nrp1_wt[volcano_df_nrp1_wt.loc[:,"gene_names"]=="CLA4"]
#dnrp1_genes_positive_enriched[dnrp1_genes_positive_enriched.loc[:,"gene_names"]=="whi3"]

# +

def pie_chart_enrichment(data,process,type,savefig=False):

    
    

    terms=[]
    for i in data.loc[:,"Term"].tolist():
        terms.append(i.split(" (")[0])

    ## Data to plot
    
    data2plot = data.loc[:,"Combined Score"][0:10]
    labels = terms[0:10]

    #define Seaborn color palette to use
    colors = sns.color_palette('pastel')[0:10]


    #create pie chart
    fig,axes=plt.subplots(1,1,figsize=(10,10))
    #patches,texts=plt.pie(data2plot, labels = labels, colors = colors, autopct='%.0f%%',textprops={'fontsize': 25});
    
    explode = list()
    for k in labels:
        explode.append(0.1)
        
    pie = plt.pie(data2plot, explode=explode, shadow=True, autopct='%1.1f%%', colors=colors)
    plt.legend(pie[0], labels, loc="best",fontsize=12)
    
    plt.tight_layout()
    if savefig==True:
        plt.savefig("../postprocessed-data/enrich-analysis_dnrp1/pie_chart_enrichment_"+process+"_"+type+".png")

    return 

# +
## Enrichment analysis with gseapy library of genes that have different enrichment score in the volcano plot:
import gseapy as gp 
type_gi="Negative"
goi=dnrp1_genes_negative_enriched.loc[:,"gene_names"].tolist()
yeast = gp.get_library_name(organism='Yeast')
yeast
sets=[yeast[2],yeast[5],yeast[8],yeast[11],yeast[16] ] #['GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
#'Gene_Interaction_Hubs_BioGRID_2018','Pfam_Domains_2019']
# #%% enrichment 
for i in np.arange(0,len(sets)): 

  enr = gp.enrichr(gene_list=goi,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                  description=type_gi +'_enrichment_dnrp1',
                  outdir='../postprocessed-data/enrich-analysis_dnrp1/'+type_gi+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                )



# +
## Import analyzed data and plot it :
type_gi="Negative"
out_dir='../postprocessed-data/enrich-analysis_dnrp1/'+type_gi+'/'
filename_bio="GO_Biological_Process_2018.Yeast.enrichr.reports.txt"
filename_cell="GO_Cellular_Component_2018.Yeast.enrichr.reports.txt"
filename_mol="GO_Molecular_Function_2018.Yeast.enrichr.reports.txt"

path=out_dir+filename_bio
biological_process=pd.read_csv(path,sep="\t")

path=out_dir+filename_cell
cellular_component=pd.read_csv(path,sep="\t")

path=out_dir+filename_mol
molecular_function=pd.read_csv(path,sep="\t")

pie_chart_enrichment(biological_process,"biological_process",type=type_gi,savefig=True)
>>>>>>> 1c5527a6c3a1e3cda5126ec83fb46e9a8a798f6f
# -

gi_pd_fitness_gene.loc["RPD3"]

gi_pd_fitness_gene.loc[common_pos_interactors_nrp1]



# +
consistent_neg_interactors_nrp1=list(set(nrp1_NI_sig) & set(neg_emap_check) )

consistent_pos_interactors_nrp1=list(set(nrp1_PI_sig) & set(pos_emap_check))

common_neg_interactors_nrp1=list(set(nrp1_NI) & set(neg_emap_check) )

common_pos_interactors_nrp1=list(set(nrp1_PI) & set(pos_emap_check))

common_neg_interactors_nrp1,common_pos_interactors_nrp1,consistent_pos_interactors_nrp1,consistent_neg_interactors_nrp1,len(common_neg_interactors_nrp1),len(common_pos_interactors_nrp1),len(consistent_pos_interactors_nrp1),len(consistent_neg_interactors_nrp1)

# +
from matplotlib_venn import venn2, venn2_circles

plt.figure(figsize=(10,5))
venn2(subsets = (len(nrp1_PI_sig), len(pos_emap_check), len(consistent_pos_interactors_nrp1)), 
set_labels = ('SATAy', 'existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(nrp1_PI_sig), len(pos_emap_check), len(consistent_pos_interactors_nrp1)));
plt.title("Consistency in positive interactions")

# +
from matplotlib_venn import venn2, venn2_circles

plt.figure(figsize=(10,5))
venn2(subsets = (len(nrp1_PI), len(pos_emap_check), len(common_pos_interactors_nrp1)), 
set_labels = ('SATAy', 'existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(nrp1_PI), len(pos_emap_check), len(common_pos_interactors_nrp1)));
plt.title("Consistency in positive interactions")
#
#plt.savefig("../figures/venn_pos_nrp1.png",dpi=300)
# -

plt.figure(figsize=(10,5))
venn2(subsets = (len(nrp1_NI_sig), len(neg_check), len(consistent_neg_interactors_nrp1)), 
set_labels = ('SATAy', 'existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(nrp1_NI_sig), len(neg_check), len(consistent_neg_interactors_nrp1)));
plt.title("Consistency in negative interactions")

plt.figure(figsize=(10,5))
venn2(subsets = (len(nrp1_NI), len(neg_emap_check), len(common_neg_interactors_nrp1)), 
set_labels = ('SATAy', 'existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(nrp1_NI), len(neg_emap_check), len(common_neg_interactors_nrp1)));
plt.title("Consistency in negative interactions")
#plt.savefig("../figures/venn_neg_nrp1.png",dpi=300)

pos_satay_neg_costanzo=list(set(nrp1_PI_sig) & set(neg_emap_check) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(nrp1_PI_sig), len(neg_emap_check), len(pos_satay_neg_costanzo)), 
set_labels = ('PI SATAy', 'NI existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(nrp1_PI_sig), len(neg_emap_check), len(pos_satay_neg_costanzo)));
plt.title("Misclasification I")

# +
pos_satay_neg_costanzo=list(set(nrp1_PI) & set(neg_emap_check) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(nrp1_PI), len(neg_emap_check), len(pos_satay_neg_costanzo)), 
set_labels = ('PI SATAy', 'NI existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(nrp1_PI), len(neg_emap_check), len(pos_satay_neg_costanzo)));
plt.title("Misclasification I")

#plt.savefig("../figures/venn_misclass1.png",dpi=300)
# -

neg_satay_pos_costanzo=list(set(nrp1_NI_sig) & set(pos_emap_check) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(nrp1_NI_sig), len(pos_emap_check), len(neg_satay_pos_costanzo)), 
set_labels = ('NI SATAy', 'PI existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(nrp1_NI_sig), len(pos_emap_check), len(neg_satay_pos_costanzo)));
plt.title("Misclasification II")


# +
neg_satay_pos_costanzo=list(set(nrp1_NI) & set(pos_emap_check) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(nrp1_NI), len(pos_emap_check), len(neg_satay_pos_costanzo)), 
set_labels = ('NI SATAy', 'PI existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(nrp1_NI), len(pos_emap_check), len(neg_satay_pos_costanzo)));
plt.title("Misclasification II")

#plt.savefig("../figures/venn_misclass2.png",dpi=300)

# +
## Heatmap of the SAFE scores for all backgrounds 

data_safe=pd.read_excel("../postprocessed-data/tcm-safe-custom-230516_NRP1.xlsx",sheet_name="SAFE enrichment scores",index_col=1)

data_safe.drop("ORF",axis=1,inplace=True)
data_safe.drop("Feature Qualifier",axis=1,inplace=True)
data_safe.drop("Allele",axis=1,inplace=True)

data_safe.replace(np.nan,None,inplace=True)

data_safe.to_csv("../postprocessed-data/tcm-safe-custom-230516_NRP1_analyze.csv")
# -

data_safe=pd.read_csv("../postprocessed-data/tcm-safe-custom-230516_NRP1_analyze.csv",index_col=0)
data_safe4heatmap=data_safe.iloc[:,1:]
data_safe4heatmap=data_safe4heatmap.astype(float)

# +
from matplotlib.patches import Patch
x=data_safe4heatmap.iloc[:,0:6]

labels = data_safe.Annotations
colors=sns.color_palette("colorblind", len(set(labels)))
lut = dict(zip(set(labels), colors))
row_colors = pd.DataFrame(labels)["Annotations"].map(lut)

g=sns.clustermap(x,cmap=sns.color_palette("viridis"),vmax=0.15,vmin=0.05,figsize=(10,10),row_colors=row_colors,col_cluster=False)

x0, _y0, _w, _h = g.cbar_pos
g.ax_cbar.set_position([x0, 0.8, g.ax_row_dendrogram.get_position().width, 0.2])
g.ax_cbar.set_title('Similarity score')

handles = [Patch(facecolor=lut[name]) for name in lut]
plt.legend(handles, lut, title='Annotations', loc='upper left', bbox_to_anchor=(5, 1.5),fontsize=12,ncol=2)
#sns.heatmap(x,cbar=True,cmap="PuBuGn_r",vmax=x.max().max(),vmin=x.min().min())

plt.savefig("../figures/fig_heatmap_safe_scores_NRP1.png",dpi=300,transparent=True)
# -

# ## Volcano plots

# +
path_a = r"../data/"
filelist_a = ["wt_a/WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_pergene.txt",
"wt_b/WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_pergene.txt"]
path_b = r"../data/"
filelist_b = ["dnrp1_1/dnrp1-1_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt",
"dnrp1_2/dnrp1-2_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt"]


variable = 'read_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.001 #set threshold above which p-values are regarded significant
normalize=True

trackgene_list = ["COX9","NTF2","PBY1","SMD1","SBH1","NRP1","ECM25"] # ["cdc42"]


figure_title = "$\Delta$nrp1 vs WT "

fc_interval=[1.5,-1.5]
pv_values=[3,3]

volcano_df_nrp1_wt = volcano(path_a=path_a, filelist_a=filelist_a,
            path_b=path_b, filelist_b=filelist_b, 
            variable=variable, 
            significance_threshold=significance_threshold,
            normalize=normalize,
            trackgene_list=trackgene_list,
            figure_title=figure_title)

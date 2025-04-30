# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
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

# +
# Importing pergene files 

pergene_files=[]

data_dir="../postprocessed-data/"

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




# +
import pickle
with open("../postprocessed-data/fitness_models_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

fitness_all_pd=pd.concat(b,axis=0,keys=keys)

standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str)


polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")
polarity_genes.fillna(0,inplace=True)

# +
half_life=pd.read_csv("../postprocessed-data/Gene_Protein_Half_Life.tsv",sep="\t")
half_life.index=half_life["Gene"]
half_life.drop("Gene",axis=1,inplace=True)
protein_abundance=pd.read_csv("../postprocessed-data/Gene_ProteinAbundance.tsv",sep="\t")
# remove the p at the end of the protein names

protein_abundance["Protein"]=protein_abundance["Protein"].apply(lambda x: x.split("p")[0])
protein_abundance.index=protein_abundance["Protein"]
protein_abundance.drop("Protein",axis=1,inplace=True)

positive_bem1=np.loadtxt("../postprocessed-data/positive_satay_genes_bem1.txt",dtype=str)
negative_bem1=np.loadtxt("../postprocessed-data/negative_satay_genes_bem1.txt",dtype=str)


# +
## Postprocessing


protein_abundance_selection=protein_abundance[(protein_abundance.Treatment=="untreated") & (protein_abundance.Background.isin(["S288c","W303"])) 
& (protein_abundance.Media.isin(["YEPD"])) & (protein_abundance.Method.isin(["quantitative mass spectrometry evidence"]))]

## make index upper.casefold()
protein_abundance_selection.index=protein_abundance_selection.index.map(str.upper)

# mean_abundances=protein_abundance_selection.groupby("Protein").mean()



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

mapping_index=[]

for i in data_processed.Accession:
    if i in mapping_accesion.index:
        mapping_index.append(mapping_accesion.loc[i,"Gene Names (primary)"])
        
    else:
        mapping_index.append(i)
        

data_processed.index=mapping_index
data2plot.index=mapping_index

# +
data2volcano=data_processed.loc[:, data_processed.columns.str.startswith('r_')]
data2volcano=data2volcano.iloc[:, np.r_[0:2, 4:6]] # WT vs pgal-bem1

for i in range(0,len(data2volcano)):
    data2volcano.loc[i,"log2FC"]=np.log2(data2volcano.iloc[i,0:2].mean()/data2volcano.iloc[i,2:4].mean())
    data2volcano.loc[i,"p-value"]=data_processed.loc[i,"p-value"]
    data2volcano.loc[i,"Significance4volcano"]=data_processed.loc[i,"Significance"]/10
    data2volcano.loc[i,"Accession"]=data_processed.loc[i,"Accession"]
    if data2volcano.loc[i,"log2FC"]<-0.5 :
        data2volcano.loc[i,"Significance"]="downregulated"
    elif data2volcano.loc[i,"log2FC"]>0.5:
        data2volcano.loc[i,"Significance"]="upregulated"
    else:
        data2volcano.loc[i,"Significance"]="not significant"


colors = {"not significant":"#D0D3D6","downregulated":'purple', "upregulated":'green'} # based on p-value significance 
plt.figure(figsize=(8,5))
plt.scatter(data2volcano["log2FC"],data2volcano["Significance4volcano"],alpha=0.6,
c=data2volcano["Significance"].apply(lambda x:colors[x]),s=100)
plt.xlabel("log2FC")
plt.ylabel(" p value")
plt.title("Volcano plot  pgal-bem1 vs WT")
plt.xlim(-2,2)
#plt.vlines(data2volcano.log2FC.mean(),1.2,3.5,linestyles="dashed",color="gray",label="mean")
# -

bem1KO_mass_spec=data2plot.loc[:,"i_dbem1gal_merged"]

# +
# find common proteins between bem1KO and  half_life and protein_abundance_selection

common_proteins=list(set(bem1KO_mass_spec.index).intersection(set(half_life.index)).intersection(set(protein_abundance_selection.index)))

bem1KO_mass_spec=bem1KO_mass_spec.loc[common_proteins]
half_life=half_life.loc[common_proteins]
protein_abundance_selection=protein_abundance_selection.loc[common_proteins]


# +
protein_abundance_selection_mean=[]

for i in protein_abundance_selection.index:
    protein_abundance_selection_mean.append(np.mean(protein_abundance_selection.loc[i,"Abundance"]))

protein_abundance_selection["Mean_abundance"]=protein_abundance_selection_mean
# -

protein_abundance_selection.loc["ACT1","Mean_abundance"][0]

# +
# plot bem1KO vs half_life

plt.figure(figsize=(8,5))

plt.scatter(bem1KO_mass_spec,half_life.loc[:,"Half-life"],alpha=0.6,s=100)
plt.yscale("log")
plt.ylabel("Half-life (hours)")
plt.xlabel("bem1KO_mass_spec")

# +
# plot bem1KO vs protein_abundance_selection

plt.figure(figsize=(8,5))

for i in common_proteins:
    plt.scatter(bem1KO_mass_spec.loc[i],protein_abundance_selection.loc[i,"Mean_abundance"][0],alpha=0.6,s=100,c="gray")
    
 
  
plt.yscale("log")
plt.ylabel("Mean_abundance")
plt.xlabel("bem1KO_mass_spec")

# +
## plot half_life vs protein_abundance_selection

plt.figure(figsize=(8,5))

for i in common_proteins:

    plt.scatter(half_life.loc[i,"Half-life"],protein_abundance_selection.loc[i,"Mean_abundance"][0],alpha=0.6,s=100,c="gray")

plt.yscale("log")
plt.xscale("log")
plt.ylabel("Mean_abundance")
plt.xlabel("Half-life (hours)")

# +
# find common genes between bem1KO_mass_spec and interactors 

common_proteins_PI=list(set(bem1KO_mass_spec.index).intersection(set(positive_bem1)))
common_proteins_NI=list(set(bem1KO_mass_spec.index).intersection(set(negative_bem1)))

len(common_proteins_PI),len(common_proteins_NI)

# +
# find common_proteins between interactors of bem1 and half_life and protein_abundance_selection

common_proteins_PI=list(set(positive_bem1).intersection(set(half_life.index)))
common_proteins_NI=list(set(negative_bem1).intersection(set(half_life.index)))

# +
plt.figure(figsize=(20,5))

plt.plot(half_life.loc[common_proteins_PI,"Half-life"],"*",c="purple",label="positive_bem1")

plt.plot(half_life.loc[common_proteins_NI,"Half-life"],"*",c="green",label="negative_bem1")

plt.xticks(rotation=90);
plt.ylabel("Half-life (hours)")

# +
# find common_proteins between interactors of bem1 and half_life and protein_abundance_selection

common_proteins_PI=list(set(positive_bem1).intersection(set(protein_abundance_selection.index)))

common_proteins_NI=list(set(negative_bem1).intersection(set(protein_abundance_selection.index)))

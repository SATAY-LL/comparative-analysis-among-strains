# ---
# jupyter:
#   jupytext:
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

from functions_interaction_computations import filter_fitness
from functions_interaction_computations import digenic_GI
from functions_interaction_computations import classify_GI


# +

def uncertainty_propagation(f,a,b,sigma_a,sigma_b):
    return f*np.sqrt(np.divide(sigma_a,a)**2+np.divide(sigma_b,b)**2)


# +

# Importing pergene files 

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


# +
import pickle
with open("../postprocessed-data/fitness_models_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

fitness_all_pd=pd.concat(b,axis=0,keys=keys)

standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str)


polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")
polarity_genes.fillna(0,inplace=True)


# +
# with open("../postprocessed-data/polarity_genes_Werner_name.txt","w") as f:
#     for i in polarity_genes.index:
#         f.write(i+"\n")
# f.close()
# -

# ### Fitness for interactions, just eliminate from the datatset those genes with not enough reads, insertions or flanking regions , because we have very little information about them , to incorporate them to the interaction calculation. Only for the domain correction fitness if there are not enough insertions then wwe will assign zero fitness

data_fitness=filter_fitness(fitness_all_pd,backgrounds=keys,goi=["BEM1","BEM3","NRP1"],discard=["Not enough flanking regions"],set2zero=["Not enough reads",
    "Not enough insertions"],cols=["fitness_gene","fitness_domains_corrected"])

len(data_fitness.loc["wt_merged"]),len(data_fitness.loc["bem1-aid_a"]),len(data_fitness.loc["bem1-aid_b"])

data_fitness.loc["wt_merged"].loc["BEM1"]

# +


column="fitness_gene"
arr_dbem1=np.asanyarray([data_fitness.loc["bem1-aid_a"].loc[:,column],data_fitness.loc["bem1-aid_b"].loc[:,column]],dtype="object")
arr_wt=np.asanyarray([data_fitness.loc["wt_a"].loc[:,column],data_fitness.loc["wt_b"].loc[:,column]],dtype="object")
arr_dbem1dbem3=np.asanyarray([data_fitness.loc["dbem1dbem3_a"].loc[:,column],data_fitness.loc["dbem1dbem3_b"].loc[:,column]],dtype="object")
arr_dbem3=np.asanyarray([data_fitness.loc["dbem3_a"].loc[:,column],data_fitness.loc["dbem3_b"].loc[:,column]],dtype="object")

data_fitness_wt=np.mean(arr_wt,axis=0)
data_fitness_dbem1=np.mean(arr_dbem1,axis=0)
data_fitness_dbem1dbem3=np.mean(arr_dbem1dbem3,axis=0)
data_fitness_dbem3=np.mean(arr_dbem3,axis=0)

data=[data_fitness_wt, data_fitness_dbem1, data_fitness_dbem1dbem3, data_fitness_dbem3]



#plt.savefig("../figures/fitness_distribution_dbem1_wt_dbem1dbem3_dbem3.png", dpi=300, bbox_inches="tight")

# +
##Analysis evolutionary trajectory of dbem1 cells


## plot the represssion to bem1d from L Laan eLife paper from evolutionary experiments   

growth_rate_steps=[0.01244,0.0022,0.010,0.01122]
growth_rate_steps2WT=np.divide(growth_rate_steps,growth_rate_steps[0])

growth_rate_steps_std=[0.000284623951135,0.000157371635941,0.000903009013801,0.000527856896616]



growth_rate_steps_std2wt=uncertainty_propagation(f=growth_rate_steps2WT,
    a=growth_rate_steps,b=growth_rate_steps[0],sigma_a=growth_rate_steps_std,sigma_b=growth_rate_steps_std[0])

labels=["WT","dbem1","dbem1dbem3","dbem1dbem3dnrp1"]

plt.figure(figsize=(6,6))

plt.errorbar(labels,growth_rate_steps2WT,yerr=growth_rate_steps_std2wt,fmt="o",color="gray",markersize=15,
    capsize=10,capthick=2)
# #plt.plot(labels[1:],growth_rate_steps[1:],"o",color="black",markersize=10)
# #plt.bar(x=labels[1:],height=growth_rate_steps[1:],color="gray",alpha=0.6)
plt.plot(labels,growth_rate_steps2WT,"--",color="black",alpha=0.3)
# #plt.hlines(growth_rate_steps[0],0,2,linestyles="dashed",label="Wild type ",color="purple")
plt.xticks(rotation=45);
# plt.ylabel("Growth rate compared to WT (min$^{-1}$)")
plt.grid(linewidth=0.2)


plt.tight_layout()
plt.savefig("../figures/evolutionary_trajectory_growth_rate_dots.png",dpi=300)
# -

growth_rate_steps_std2wt

# +
## take the same index across all libraries

index_arrays=np.array([data_fitness.loc["wt_merged"].index,data_fitness.loc["bem1-aid_merged"].index,
data_fitness.loc["bem1-aid-dbem3_a"].index,data_fitness.loc["bem1-aid-dbem3_b"].index,data_fitness.loc["dbem3_merged"].index],
dtype="object")

d1=set.intersection(*map(set,index_arrays))
# because Bem3 will be artficially excluded since its fitness its zero in the dbem1dbem3 because of the strain background 

# +
d1_index=np.array(list(d1))

# index_BEM1=np.where(d1_index=="BEM1")[0][0]
# index_BEM3=np.where(d1_index=="BEM3")[0][0]
index_NRP1=np.where(d1_index=="NRP1")[0][0]


data_wt=data_fitness.loc["wt_merged"].loc[d1,"fitness_gene"]
data_dbem1=data_fitness.loc["bem1-aid_merged"].loc[d1,"fitness_gene"]
data_dbem3=data_fitness.loc["dbem3_merged"].loc[d1,"fitness_gene"]
data_dbem1dbem3=np.mean(np.asanyarray([data_fitness.loc["bem1-aid-dbem3_a"].loc[d1,"fitness_gene"],
data_fitness.loc["bem1-aid-dbem3_b"].loc[d1,"fitness_gene"]],dtype="object"),axis=0)




# +
x=list_data_pd.loc["wt_merged"]
y=list_data_pd.loc["bem1-aid_merged"]
za=list_data_pd.loc["bem1-aid-dbem3_a"]
zb=list_data_pd.loc["bem1-aid-dbem3_b"]


x.index=x.loc[:,"Gene name"]
y.index=y.loc[:,"Gene name"]
za.index=za.loc[:,"Gene name"]
zb.index=zb.loc[:,"Gene name"]



bem1_insertions=x.loc["BEM1","Insertions"]
bem3_insertions=y.loc["BEM3","Insertions"]
nrp1_insertions=np.mean([za.loc["NRP1","Insertions"],zb.loc["NRP1","Insertions"]])

# +


# ## Standarize the data as having mu=0 and std=1 ## NOT A GOOD IDEA IN THIS CASE BECAUSE https://stats.stackexchange.com/questions/19216/variables-are-often-adjusted-e-g-standardised-before-making-a-model-when-is
# data_wt_norm=scaler(data_wt)
# data_dbem1_norm=scaler(data_dbem1)
# data_dbem1dbem3_norm=scaler(data_dbem1dbem3)
fitness_bem1_wt=(data_fitness.loc["wt_merged"].loc["BEM1","fitness_gene"]+data_fitness.loc["wt_merged"].loc["BEM1","fitness_domains_corrected"])/2

fitness_bem3_wt=data_wt.loc["BEM3"]

## Normalize the datasetss such as the median values corresponds to the fitness of the knockout of the gene of interest 
data_wt_norm=data_wt

data_dbem1_norm=data_dbem1*fitness_bem1_wt/(np.median(data_dbem1))# because their median is around 0 
data_dbem3_norm=data_dbem3*fitness_bem3_wt/(np.median(data_dbem3))

fitness_bem3_dbem1=data_dbem1_norm.loc["BEM3"]

data_dbem1dbem3_norm=data_dbem1dbem3*fitness_bem3_dbem1/(np.median(data_dbem1dbem3))







# +
## Reconstitution of evolutionary trajectory by SATAY experiments 

## Get fitness values of bem1 in WT, bem3 in bem1-aid and nrp1 in bem1-aid-dbem3

## get the average fitness from the domain correction and from the whole gene 
# fitness_bem1=np.mean([data_fitness.loc["wt_merged"].loc["BEM1","fitness_gene"],
# data_fitness.loc["wt_merged"].loc["BEM1","fitness_domains_corrected"]])

fitness_bem1=fitness_bem1_wt
fitness_bem1_std= data_fitness.loc["wt_merged"].loc["BEM1","fitness_gene_std"]/np.sqrt(bem1_insertions)


fitness_bem3_dbem1_std= data_fitness.loc["bem1-aid_merged"].loc["BEM3","fitness_gene_std"]/np.sqrt(bem3_insertions)

fitness_nrp1_dbem1_dbem3=data_dbem1dbem3_norm[index_NRP1]

fitness_nrp1_dbem1_dbem3_std= data_fitness.loc["bem1-aid-dbem3_a"].loc["NRP1","fitness_gene_std"]/np.sqrt(nrp1_insertions)

fitness_bem1,fitness_bem3_dbem1,fitness_nrp1_dbem1_dbem3,fitness_bem1_std,fitness_bem3_dbem1_std,fitness_nrp1_dbem1_dbem3_std

# +
##Analysis evolutionary trajectory of dbem1 cells


## Plot the correlation between growth rate values fm the evolutionary experiments and the SATAY experiments   


growth_rate_steps_satay=[np.median(data_wt_norm),fitness_bem1,
fitness_bem3_dbem1,fitness_nrp1_dbem1_dbem3]
growth_rate_steps_satay2WT=np.divide(growth_rate_steps_satay,growth_rate_steps_satay[0])
growth_rate_steps_std=[np.std(data_wt),fitness_bem1_std,fitness_bem3_dbem1_std,fitness_nrp1_dbem1_dbem3_std]

labels=["WT","$\Delta$bem1","$\Delta$bem1$\Delta$bem3","$\Delta$bem1$\Delta$bem3$\Delta$nrp1"]

fig,ax=plt.subplots(figsize=(5,5))




plt.xticks(rotation=90);

ax.errorbar(labels,growth_rate_steps_satay,yerr=growth_rate_steps_std,linestyle="None",color="black",alpha=0.3,capsize=5)
ax.scatter(labels,growth_rate_steps_satay,marker="o",color=["gray","blue","orange","purple"],alpha=0.8,s=250)
ax.plot(labels,growth_rate_steps_satay,linestyle="dashed",color="black",linewidth=0.2)
# ax.annotate("WT",xy=(0,growth_rate_steps_satay[0]),xytext=(0+0.1,growth_rate_steps_satay[0]),fontsize=14)
# ax.annotate("$\Delta$bem1",xy=(1,growth_rate_steps_satay[1]),xytext=(1+0.1,growth_rate_steps_satay[1]),fontsize=14)
# ax.annotate("$\Delta$bem1$\Delta$bem3",xy=(2,growth_rate_steps_satay[2]),xytext=(2+0.1,growth_rate_steps_satay[2]),fontsize=14)
# ax.annotate("$\Delta$bem1$\Delta$bem3\n$\Delta$nrp1",xy=(3,growth_rate_steps_satay[3]),xytext=(3+0.1,growth_rate_steps_satay[3]),fontsize=14)
# # ax.plot([-1,1.5],[0,0.01],linestyle="dashed",color="black",linewidth=0.4)
ax.set_ylabel("Standarized fitness values")
#ax.set_ylabel("Growth rate (min$^{-1}$)")
ax.grid(linewidth=0.2)
ax.set_yscale("log")


#fig.savefig("../figures/fig_comparison_growth_rate_trjectory_vs_satay.png",dpi=300,bbox_inches="tight")


# +
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
growth_rate_steps_satay=[np.median(data_wt_norm),fitness_bem1,
fitness_bem3_dbem1,fitness_nrp1_dbem1_dbem3]

growth_rate_steps_std=[np.std(data_wt),fitness_bem1_std,fitness_bem3_dbem1_std,fitness_nrp1_dbem1_dbem3_std]


fig, ax1 = plt.subplots(figsize=(8,8))

labels=["WT","$\Delta$bem1","$\Delta$bem1$\Delta$bem3"]
labels2scatter=["WT","$\Delta$bem1","$\Delta$bem1$\Delta$bem3","$\Delta$bem1$\Delta$bem3$\Delta$nrp1"]
box=ax1.boxplot([data_wt_norm,data_dbem1_norm,data_dbem1dbem3_norm],labels=labels,
showfliers=False,patch_artist=True,notch=True,boxprops=dict(alpha=.2),whiskerprops=dict(alpha=.2),capprops=dict(alpha=.2));
plt.ylabel("Fitness values")

ax1.scatter(labels2scatter,growth_rate_steps_satay,marker="o",color=["gray","blue","orange","purple"],s=250)
ax1.errorbar(labels2scatter,growth_rate_steps_satay,yerr=growth_rate_steps_std,linestyle="None",color="black",alpha=0.8,capsize=5)
#ax1.plot(labels2scatter,growth_rate_steps_satay,linestyle="dashed",color="black",linewidth=0.2)

plt.grid(linewidth=0.2)
colors = ['black', "blue", 'orange']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)




plt.savefig("../figures/fig_boxplots_fitness_distributions_satay.png",dpi=300,bbox_inches="tight")

# +

## Plot a  3D plot with the three fitness vectors 

from mpl_toolkits.mplot3d import Axes3D


fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')
x=data_wt_norm
y=data_dbem1_norm
z=data_dbem1dbem3_norm

colo = x
  
 
# configuring colorbar
color_map = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
color_map.set_array(colo)

ax.scatter3D(x,y,z, c=colo,s=200,marker="s");


cbar=plt.colorbar(color_map,fraction=0.03, pad=0.1)

cbar.set_label('Fitness WT')

ax.set_xlabel('WT')
ax.set_ylabel('$\Delta$bem1')
ax.set_zlabel('$\Delta$bem1$\Delta$bem3')
# -

gi_bem1=digenic_GI(data=data_fitness,goi="BEM1",col_fitness="fitness_gene",backg=["wt_a","wt_b","bem1-aid_a","bem1-aid_b"])
gi_bem3=digenic_GI(data=data_fitness,goi="BEM3",col_fitness="fitness_gene",backg=["wt_a","wt_b","dbem3_a","dbem3_b"])

# +
gi_pd_fitness_gene=pd.DataFrame.from_dict(gi_bem1,orient="index")

gi_pd_fitness_gene["gene_names"]=gi_pd_fitness_gene.index

gi_pd_fitness_gene_bem1d=gi_pd_fitness_gene
gi_pd_fitness_gene_bem1d

# +
## import genes from the cellcycle pathway and mapk pathway from wikipathway

cellcycle=pd.read_csv("../postprocessed-data/Cell cycle and cell division - Saccharomyces cerevisiae default  node.csv",sep="\t")

cellcycle=cellcycle.name.tolist()

mapk=pd.read_csv("../postprocessed-data/MAPK signaling-pathway-Saccharomyces-cerevisiae-default-node.csv",sep="\t")

mapk=mapk.name.tolist()

lipid=pd.read_csv("../postprocessed-data/Sphingolipid metabolism - Saccharomyces cerevisiae default  node.csv",sep="\t")
lipid=lipid.name.tolist()

glycolisis=pd.read_csv("../postprocessed-data/Glycolysis and gluconeogenesis - Saccharomyces cerevisiae default  node.csv",sep="\t")
glycolisis=glycolisis.name.tolist()

fatty_acid=pd.read_csv("../postprocessed-data/Fatty acid biosynthesis, initial steps - Saccharomyces cerevisiae default  node.csv",sep="\t")
fatty_acid=fatty_acid.name.tolist()

# +
all_pathways=cellcycle+mapk+lipid+glycolisis+fatty_acid

all_pathways=list(set(all_pathways))

# make all the gene names uppercase

all_pathways=[i.upper() for i in all_pathways]

all_pathways_gi=gi_pd_fitness_gene[gi_pd_fitness_gene.loc[:,"gene_names"].isin(all_pathways)]

#all_pathways_gi.to_csv("../postprocessed-data/pathways_bem1_gi.csv",sep="\t")
# -

all_pathways_gi.sort_values(by="fold_change",ascending=False)

x=gi_pd_fitness_gene_bem1d.loc[:,"fold_change"].sort_values(ascending=False)
# np.where(x.index=="BEM3")
gi_pd_fitness_gene_bem1d.loc["BEM3"]

# +
gi_pd_fitness_gene=pd.DataFrame.from_dict(gi_bem3,orient="index")

gi_pd_fitness_gene["gene_names"]=gi_pd_fitness_gene.index

gi_pd_fitness_gene_bem3d=gi_pd_fitness_gene
gi_pd_fitness_gene_bem3d

# +
bem1_PI_sig,bem1_NI_sig,bem1_PI,bem1_NI,bem1_PI_all,bem1_NI_all=classify_GI(gi_pd_fitness_gene_bem1d,col="fold_change")

bem3_PI_sig,bem3_NI_sig,bem3_PI,bem3_NI,bem3_PI_all,bem3_NI_all=classify_GI(gi_pd_fitness_gene_bem3d,col="fold_change")

bem1_PI2export=gi_pd_fitness_gene_bem1d.loc[bem1_PI]
bem1_PI2export=bem1_PI2export[bem1_PI2export["p_statistic"]<0.3]

bem1_NI2export=gi_pd_fitness_gene_bem1d.loc[bem1_NI]
bem1_NI2export=bem1_NI2export[bem1_NI2export["p_statistic"]<0.3]

for i in gi_pd_fitness_gene_bem1d.index:
    if i in bem1_NI_sig or i in bem1_NI2export.index:
        gi_pd_fitness_gene_bem1d.loc[i,"significance4FC"]="neg"
    # elif i in  bem1_NI2export.index:
    #     gi_pd_fitness_gene_bem1d.loc[i,"significance4FC"]="neg"
    elif i in bem1_PI_sig or i in bem1_PI2export.index:
        gi_pd_fitness_gene_bem1d.loc[i,"significance4FC"]="pos"
    # elif i in bem1_PI2export.index:
    #     gi_pd_fitness_gene_bem1d.loc[i,"significance4FC"]="pos"
    else:
        gi_pd_fitness_gene_bem1d.loc[i,"significance4FC"]="none"


for i in gi_pd_fitness_gene_bem3d.index:
    if i in bem3_NI_sig:
        gi_pd_fitness_gene_bem3d.loc[i,"significance4FC"]="neg"
    elif i in bem3_PI_sig:
        gi_pd_fitness_gene_bem3d.loc[i,"significance4FC"]="pos"
    else:
        gi_pd_fitness_gene_bem3d.loc[i,"significance4FC"]="none"
# -

len(bem1_PI_sig),len(bem1_NI_sig),len(bem1_PI),len(bem1_NI),len(bem1_PI_all),len(bem1_NI_all)

# +
gi_pd_fitness_gene_bem1d.loc[:,"fold_change"].hist(bins=50)

plt.vlines(gi_pd_fitness_gene_bem1d.loc[:,"fold_change"].median(),0,300,color="red")


# -

gi_pd_fitness_gene_bem1d[gi_pd_fitness_gene_bem1d.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)

# +
# look for the fold change of the standard_essentials genes

gi_standard_essentials=defaultdict(dict)
gi_pd_fitness_gene=gi_pd_fitness_gene_bem1d
for i in standard_essentials:
    if i in gi_pd_fitness_gene.index:
        gi_standard_essentials[i]["fold_change"]=gi_pd_fitness_gene.loc[i,"fold_change"]
        gi_standard_essentials[i]["p_statistic"]=gi_pd_fitness_gene.loc[i,"p_statistic"]
        gi_standard_essentials[i]["significance"]=gi_pd_fitness_gene.loc[i,"significance"]
        gi_standard_essentials[i]["significance4FC"]=gi_pd_fitness_gene.loc[i,"significance4FC"]
        

gi_standard_essentials_pd=pd.DataFrame.from_dict(gi_standard_essentials,orient="index")
gi_standard_essentials_pd_dbem1=gi_standard_essentials_pd
gi_standard_essentials_pd_dbem1.sort_values(by="fold_change",ascending=False)




# +
from annotate_volcano import annotate_volcano

volcano_df=gi_pd_fitness_gene_bem1d

fig=annotate_volcano(volcano_df,figure_title="Interactors of bem1 in WT")

plt.savefig("../figures/volcano_bem1_wt.png",dpi=300)
# -

gi_pd_fitness_gene_bem1d.loc[bem1_NI2export.index].sort_values(by="fold_change",ascending=False)

# +
from annotate_volcano import annotate_volcano   #import annotate_volcano function
volcano_df=gi_pd_fitness_gene_bem1d
trackgene_list=["BEM3","NRP1"]
fig=annotate_volcano(volcano_df,figure_title="Interactors of bem1 in WT",trackgene_list=trackgene_list)

plt.savefig("../figures/fig_volcano_interactors_bem1_with annotations.png",dpi=300,transparent=True)
# -

gi_pd_fitness_gene_bem3d[gi_pd_fitness_gene_bem3d.loc[:, "p_statistic"]<0.05].sort_values(by="fold_change",ascending=False)

# +
from annotate_volcano import annotate_volcano

volcano_df=gi_pd_fitness_gene_bem3d

fig=annotate_volcano(volcano_df,figure_title="Interactors of bem3 in WT")

plt.savefig("../figures/volcano_bem3_wt.png",dpi=300)

# +
bem1_PI2export=gi_pd_fitness_gene_bem1d.loc[bem1_PI]
bem1_PI2export=bem1_PI2export[bem1_PI2export["p_statistic"]<0.3]

bem1_NI2export=gi_pd_fitness_gene_bem1d.loc[bem1_NI]
bem1_NI2export=bem1_NI2export[bem1_NI2export["p_statistic"]<0.3]


# -

len(bem1_PI2export),len(bem1_NI2export)

gi_pd_fitness_gene_bem1d.loc[bem1_PI].sort_values(by="fold_change",ascending=True)

len(bem1_PI2export),len(bem1_NI2export)

# +
## Export a tab delimited table with the information of the bem1 interactors

table_PI=pd.DataFrame(index=bem1_PI2export.index)
table_PI["score"]=bem1_PI2export.loc[:,"fold_change"]
table_PI["p_value"]=bem1_PI2export.loc[:,"p_statistic"]

table_NI=pd.DataFrame(index=bem1_NI2export.index)
table_NI["score"]=bem1_NI2export.loc[:,"fold_change"]
table_NI["p_value"]=bem1_NI2export.loc[:,"p_statistic"]

table_BEM1_interactors=pd.concat([table_PI,table_NI],axis=0)
table_BEM1_interactors["target_node"]=len(table_PI)*["BEM1"]+len(table_NI)*["BEM1"]
table_BEM1_interactors.to_csv("../postprocessed-data/table_BEM1_interactors.tsv",sep="\t")

# +
## export to a txt the gene names of the significant positive and negative regulators of nrp1

with open("../postprocessed-data/positive_satay_genes_bem1.txt","w") as f:
    for i in bem1_PI2export.index:
        f.write(i+"\n")
f.close()
with open("../postprocessed-data/negative_satay_genes_bem1.txt","w") as f:
    for i in bem1_NI2export.index:
        f.write(i+"\n")
f.close()

with open("../postprocessed-data/positive_satay_genes_bem3.txt","w") as f:
    for i in bem3_PI_sig:
        f.write(i+"\n")
f.close()
with open("../postprocessed-data/negative_satay_genes_bem3.txt","w") as f:
    for i in bem3_NI_sig:
        f.write(i+"\n")
f.close()




# +
## interaction score for polarity genes 

polarity_genes_gi=defaultdict(dict)

gi_pd_fitness_gene=gi_pd_fitness_gene_bem1d
for i in polarity_genes.index:
    if i in gi_pd_fitness_gene.index:
        polarity_genes_gi[i]["p_statistic"]=gi_pd_fitness_gene.loc[i,"p_statistic"]
        polarity_genes_gi[i]["significance"]=gi_pd_fitness_gene.loc[i,"significance"]
        polarity_genes_gi[i]["significance4FC"]=gi_pd_fitness_gene.loc[i,"significance4FC"]
        polarity_genes_gi[i]["fold_change"]=gi_pd_fitness_gene.loc[i,"fold_change"]
    
polarity_genes_gi_pd=pd.DataFrame.from_dict(polarity_genes_gi,orient="index")



polarity_genes_gi_pd_dbem1=polarity_genes_gi_pd

#len(polarity_genes_gi_pd)
polarity_genes_gi_pd_dbem1.sort_values(by="fold_change",ascending=False)

# +
## interaction score for polarity genes 

polarity_genes_gi=defaultdict(dict)

gi_pd_fitness_gene=gi_pd_fitness_gene_bem3d
for i in polarity_genes.index:
    if i in gi_pd_fitness_gene.index:
        polarity_genes_gi[i]["p_statistic"]=gi_pd_fitness_gene.loc[i,"p_statistic"]
        polarity_genes_gi[i]["significance"]=gi_pd_fitness_gene.loc[i,"significance"]
        polarity_genes_gi[i]["significance4FC"]=gi_pd_fitness_gene.loc[i,"significance4FC"]
        polarity_genes_gi[i]["fold_change"]=gi_pd_fitness_gene.loc[i,"fold_change"]
    
polarity_genes_gi_pd=pd.DataFrame.from_dict(polarity_genes_gi,orient="index")



polarity_genes_gi_pd_dbem3=polarity_genes_gi_pd

#len(polarity_genes_gi_pd)
polarity_genes_gi_pd_dbem3.sort_values(by="fold_change",ascending=False)

# +
## Plot a heatmap with the GI scores for the polarity_genes

import seaborn as sns
x=polarity_genes_gi_pd_dbem1[polarity_genes_gi_pd_dbem1.index!="BEM1"]
x=x.loc[:,["fold_change"]]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,5))
# array2heatmap[~np.isfinite(array2heatmap)] = 0
sns.heatmap(x.sort_values(by="fold_change",ascending=False),yticklabels=x.sort_values(by="fold_change",ascending=False).index,
cbar=True,annot=True)
# -

# ### Comparison with existing BEM1 interactors 

# +
## comparison with existing data

sl_bem1=pd.read_excel("../postprocessed-data/BEM1_genetic_interactions_filtered_by_synt.xlsx",skiprows=8)
neg_bem1=pd.read_excel("../postprocessed-data/BEM1_genetic_interactions_filtered_by_negat.xlsx",skiprows=8)
pos_bem1=pd.read_excel("../postprocessed-data/BEM1_genetic_interactions_filtered_by_pos.xlsx",skiprows=8)

## Filtering interactor by Costanzo data 
# x=pos_bem1[pos_bem1.Reference=="Costanzo M"].loc[:,"Interactor.1"]
x=pos_bem1.loc[:,"Interactor.1"] # no filtering by Costanzo data
y=pos_bem1[pos_bem1.Action=="Hit"].loc[:,"Interactor.1"]

z=set(x).intersection(set(y))

pos_bem1_genes=z



# x=neg_bem1[neg_bem1.Reference=="Costanzo M"].loc[:,"Interactor.1"]
x=neg_bem1.loc[:,"Interactor.1"] # no filtering by Costanzo data
y=neg_bem1[neg_bem1.Action=="Hit"].loc[:,"Interactor.1"]
z=set(x).intersection(set(y))
neg_bem1_genes=z


x=sl_bem1[sl_bem1.Assay=="Synthetic Lethality"].loc[:,"Interactor.1"]
y=sl_bem1[sl_bem1.Action=="Hit"].loc[:,"Interactor.1"]
z=set(x).intersection(set(y))

sl_bem1_genes=z



neg_sl_bem1_genes=neg_bem1_genes.union(sl_bem1_genes)

len(neg_sl_bem1_genes),len(pos_bem1_genes)


# +
sl_in_dbem1=[]

for i in z:
    if i in gi_pd_fitness_gene_bem1d.index:
        sl_in_dbem1.append(i)

gi_pd_fitness_gene_bem1d.loc[sl_in_dbem1].sort_values(by="fold_change",ascending=False)
# -

1.47/1.8

gi_pd_fitness_gene_bem1d.loc["MIG3"]

# +
## Exclude known interactors that were not analyzed in the study because their fitness in satay is annotated as zero, for example :
neg_check=[]
pos_check=[]
for i in list(neg_sl_bem1_genes):
    if i in gi_pd_fitness_gene_bem1d.index:
        neg_check.append(i)
for i in list(pos_bem1_genes):
    if i  in gi_pd_fitness_gene_bem1d.index:
        pos_check.append(i)

len(neg_check),len(pos_check),len(neg_sl_bem1_genes),len(pos_bem1_genes)

# +
neg_sl_bem1_genes=neg_check
pos_bem1_genes=pos_check

consistent_neg_interactors=list(set(bem1_NI2export) & set(neg_sl_bem1_genes) )

consistent_pos_interactors=list(set(bem1_PI2export) & set(pos_bem1_genes))


max_consistency_neg=list(set(bem1_NI) & set(neg_sl_bem1_genes) )

max_consistency_pos=list(set(bem1_PI) & set(pos_bem1_genes))



len(consistent_neg_interactors),len(consistent_pos_interactors),len(max_consistency_neg),len(max_consistency_pos),len(pos_bem1_genes)


# -

max_consistency_pos

1/7,len(max_consistency_pos)

len(max_consistency_pos)/len(pos_bem1_genes),len(consistent_neg_interactors)/len(neg_sl_bem1_genes),len(max_consistency_neg)/len(neg_sl_bem1_genes)

# +
from matplotlib_venn import venn2, venn2_circles

plt.figure(figsize=(10,5))
venn2(subsets = (len(bem1_PI_sig), len(pos_bem1_genes), len(consistent_pos_interactors)), 
set_labels = ('SATAY', 'Existing'), set_colors=('purple', 'pink'), alpha = 0.5, normalize_to=1);
# venn2_circles(subsets = (len(pos_satay_signif_special), len(pos_bem1_genes), len(max_consistency_pos)),
#  linestyle='dashed', linewidth=1, color="grey");
plt.title("Consistency in positive interactions")

plt.savefig("../figures/venn_pos_bem1.png",dpi=300)

# +
from matplotlib_venn import venn2, venn2_circles

plt.figure(figsize=(10,5))
venn2(subsets = (len(bem1_PI), len(pos_bem1_genes), len(max_consistency_pos)), 
set_labels = ('SATAY', 'Existing'), set_colors=('purple', 'pink'), alpha = 0.5, normalize_to=1);
# venn2_circles(subsets = (len(pos_satay_signif_special), len(pos_bem1_genes), len(max_consistency_pos)),
#  linestyle='dashed', linewidth=1, color="grey");
plt.title("Consistency in positive interactions")

plt.savefig("../figures/max_consistency_venn_pos_bem1.png",dpi=300)

# +
from matplotlib_venn import venn2, venn2_circles

plt.figure(figsize=(10,5))
venn2(subsets = (len(bem1_NI_sig), len(neg_sl_bem1_genes), len(consistent_neg_interactors)), 
set_labels = ('SATAy', 'Existing'), set_colors=('green', 'cyan'), alpha = 0.5);
#venn2_circles(subsets = (len(pos_satay), len(pos_bem1_genes), len(common_pos_interactors)));
plt.title("Consistency in negative interactions")

plt.savefig("../figures/venn_neg_bem1.png",dpi=300)

# +
from matplotlib_venn import venn2, venn2_circles

plt.figure(figsize=(10,5))
venn2(subsets = (len(bem1_NI), len(neg_sl_bem1_genes), len(max_consistency_neg)), 
set_labels = ('SATAy', 'Existing'), set_colors=('green', 'cyan'), alpha = 0.5);
#venn2_circles(subsets = (len(pos_satay), len(pos_bem1_genes), len(common_pos_interactors)));
plt.title("Consistency in negative interactions")

plt.savefig("../figures/max_consistency_venn_neg_bem1.png",dpi=300)

# +
pos_satay_neg=list(set(bem1_PI_sig) & set(neg_sl_bem1_genes) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(bem1_PI_sig), len(neg_sl_bem1_genes), len(pos_satay_neg)), 
set_labels = ('PI SATAy', 'NI Existing'), set_colors=('purple', 'cyan'), alpha = 0.5);
#venn2_circles(subsets = (len(pos_satay), len(neg_bem1_costanzo_genes), len(pos_satay_neg_costanzo)));
plt.title("Misclasification I")

plt.savefig("../figures/venn_pos_satay_neg_existing_bem1.png",dpi=300)
# -

pos_satay_neg_max=list(set(bem1_PI) & set(neg_sl_bem1_genes) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(bem1_PI), len(neg_sl_bem1_genes), len(pos_satay_neg_max)), 
set_labels = ('PI SATAy', 'NI Existing'), set_colors=('purple', 'cyan'), alpha = 0.5);
#venn2_circles(subsets = (len(pos_satay), len(neg_bem1_costanzo_genes), len(pos_satay_neg_costanzo)));
plt.title("Misclasification I")
pos_satay_neg_max
#plt.savefig("../figures/max_venn_pos_satay_neg_existing_bem1.png",dpi=300)

# +
neg_satay_pos=list(set(bem1_NI_sig) & set(pos_bem1_genes) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(bem1_NI_sig), len(pos_bem1_genes), len(neg_satay_pos)), 
set_labels = ('NI SATAy', 'PI Existing'), set_colors=('green', 'pink'), alpha = 0.5);
#venn2_circles(subsets = (len(neg_satay), len(pos_bem1_costanzo_genes), len(neg_satay_pos_costanzo)));
plt.title("Misclasification II")

plt.savefig("../figures/venn_neg_satay_pos_existing_bem1.png",dpi=300)


# +
neg_satay_pos_max=list(set(bem1_NI) & set(pos_bem1_genes) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(bem1_NI), len(pos_bem1_genes), len(neg_satay_pos_max)), 
set_labels = ('NI SATAy', 'PI Existing'), set_colors=('green', 'pink'), alpha = 0.5);
#venn2_circles(subsets = (len(neg_satay), len(pos_bem1_costanzo_genes), len(neg_satay_pos_costanzo)));
plt.title("Misclasification II")

plt.savefig("../figures/max_venn_neg_satay_pos_existing_bem1.png",dpi=300)

# +

gi_pd_fitness_gene_bem3d.sort_values(by="fold_change",ascending=False)

# +
from annotate_volcano import annotate_volcano   #import annotate_volcano function
volcano_df=gi_pd_fitness_gene_bem3d
trackgene_list=["ARR2","NRP1","PBY1","ECM25","BEM1","FAR1","WHI5"]
fig=annotate_volcano(volcano_df,figure_title="Interactors of dbem3 in WT",trackgene_list=trackgene_list)

plt.savefig("../figures/fig_volcano_interactors_bem3_with annotations.png",dpi=300,transparent=True)
# -

# #### Interactors in the dbem1dbem3 (compared to WT)

# #### Formula trigenic interactions 
#
# $e_{ijk}=e_{i|jk}=f_{ijk}-f_{ij}f_{k}-f_{ik}f_{j}-f_{jk}f_{i}+2f_{i}f_{j}f_{k}$
#
# - In our case i=$\Delta$ genex ,j=$\Delta$ bem1 ,k= $\Delta$ bem3
# - all the double fitness are based on the bem1d or bem3d backgrounds , because that are the libraries we have from SATAY
#
# $e_{x13}=e_{x|13}=f_{x|13}-f_{x|1}f_{3|WT}-f_{x|3}f_{1|WT}-f_{1|3}f_{x|WT}+2f_{x|WT}f_{1|WT}f_{3|WT}$
#

# +
## Normalize all fitness libraries accordingly the fitness of the knockout of interest 

backg=["wt_a","wt_b","bem1-aid_a","bem1-aid_b","dbem3_a","dbem3_b","dbem1dbem3_a","dbem1dbem3_b"]

data_backg=data_fitness.loc[backg]
data_fitness2interact=data_backg[data_backg.loc[:,"fitness_gene"]!=0]

genej="BEM1"
genek="BEM3"
col_fitness="fitness_gene"

## fitness bem1 in wt
genej_f_a=0.5*(data_fitness2interact.loc[backg[0]].loc[genej,col_fitness]+data_fitness2interact.loc[backg[0]].loc[genej,"fitness_domains_corrected"])
genej_f_b=0.5*(data_fitness2interact.loc[backg[1]].loc[genej,col_fitness]+data_fitness2interact.loc[backg[1]].loc[genej,"fitness_domains_corrected"])

## fitness bem3 in wt
genek_f_a=(data_fitness2interact.loc[backg[0]].loc[genek,col_fitness])
genek_f_b=(data_fitness2interact.loc[backg[1]].loc[genek,col_fitness])

### fitness normalization in the bem1-aid library to the median of the fitness of dbem1 in WT , and the dbem1dbem3 library 
# to the median of the fitness of dbem1dbem3 in WT

data_dgenej_a=data_fitness2interact.loc[backg[2]].loc[:,col_fitness]*genej_f_a/2**(np.median(data_backg.loc[backg[2]].loc[:,col_fitness]))
data_dgenej_b=data_fitness2interact.loc[backg[3]].loc[:,col_fitness]*genej_f_b/2**(np.median(data_backg.loc[backg[3]].loc[:,col_fitness]))

fitness_dgenek_dgenej_a=data_dgenej_a.loc[genek] # fitness of bem3 in the bem1-aid_a background
fitness_dgenek_dgenej_b=data_dgenej_b.loc[genek] # fitness of bem3 in the bem1-aid_b background

data_dgenek_a=data_fitness2interact.loc[backg[4]].loc[:,col_fitness]*genek_f_a/(np.median(data_backg.loc[backg[4]].loc[:,col_fitness]))
data_dgenek_b=data_fitness2interact.loc[backg[5]].loc[:,col_fitness]*genek_f_b/(np.median(data_backg.loc[backg[5]].loc[:,col_fitness]))

fitness_dgenej_dgenek_a=data_dgenek_a.loc[genej] # fitness of bem1 in the dbem3_a background
fitness_dgenej_dgenek_b=data_dgenek_b.loc[genej] # fitness of bem1 in the dbem3_b background

# data_dgenej_dgenek_a=data_fitness2interact.loc[backg[6]].loc[:,col_fitness]*0.5*(fitness_dgenek_dgenej_a+fitness_dgenej_dgenek_a)/(np.median(data_fitness.loc[backg[6]].loc[:,col_fitness]))
# data_dgenej_dgenek_b=data_fitness2interact.loc[backg[7]].loc[:,col_fitness]*0.5*(fitness_dgenek_dgenej_b+fitness_dgenej_dgenek_b)/(np.median(data_fitness.loc[backg[7]].loc[:,col_fitness]))
data_dgenej_dgenek_a=data_fitness2interact.loc[backg[6]].loc[:,col_fitness]*0.5*(fitness_dgenek_dgenej_a+fitness_dgenej_dgenek_a)/(np.median(data_backg.loc[backg[6]].loc[:,col_fitness]))
data_dgenej_dgenek_b=data_fitness2interact.loc[backg[7]].loc[:,col_fitness]*0.5*(fitness_dgenek_dgenej_b+fitness_dgenej_dgenek_b)/(np.median(data_backg.loc[backg[7]].loc[:,col_fitness]))

## interaction scores from each gene to bem1 and bem3
e_genex_genej_a=gi_pd_fitness_gene_bem1d.loc[:,"e_a"]
e_genex_genej_b=gi_pd_fitness_gene_bem1d.loc[:,"e_b"]

e_genek_genej_a=e_genex_genej_a.loc[genek]
e_genek_genej_b=e_genex_genej_b.loc[genek]

e_genex_genek_a=gi_pd_fitness_gene_bem3d.loc[:,"e_a"]
e_genex_genek_b=gi_pd_fitness_gene_bem3d.loc[:,"e_b"]

e_genej_genek_a=e_genex_genek_a.loc[genej]
e_genej_genek_b=e_genex_genek_b.loc[genej]


# -

plt.hist(data_dgenej_dgenek_a,bins=70,alpha=0.6);
plt.vlines(0.5*(fitness_dgenek_dgenej_a+fitness_dgenej_dgenek_a),0,600,color="blue")
plt.hist(data_dgenej_dgenek_b,bins=70,alpha=0.2);
plt.vlines(0.5*(fitness_dgenek_dgenej_b+fitness_dgenej_dgenek_b),0,600,color="orange")

# +
## Trigenic interaction scores

intersection_genes=list(set(e_genex_genej_a.index)& set(e_genex_genek_a.index)& set(e_genex_genej_b.index)& set(e_genex_genek_b.index)
& set(data_dgenej_dgenek_a.index)& set(data_dgenej_dgenek_b.index))


significance_threshold = 0.1 #set significance threshold
gi=defaultdict(dict)
ttest_tval_list = [np.nan]*2 #initialize list for storing t statistics
ttest_pval_list = [np.nan]*2 #initialize list for storing p-values
signif_thres_list = False #initialize boolean list for indicating datapoints with p-value above threshold
fc_list = [np.nan]*2
for gene in intersection_genes :
    geneX=gene
    geneX_f_a=data_fitness2interact.loc[backg[0]].loc[geneX,col_fitness]
    geneX_f_b=data_fitness2interact.loc[backg[1]].loc[geneX,col_fitness]

    e_x1a=e_genex_genej_a.loc[geneX]
    e_x1b=e_genex_genej_b.loc[geneX]

    e_x3a=e_genex_genek_a.loc[geneX]
    e_x3b=e_genex_genek_b.loc[geneX]

    genexgenejgenek_f_a=data_dgenej_dgenek_a.loc[geneX]
    genexgenejgenek_f_b=data_dgenej_dgenek_b.loc[geneX]

    
    variable_a_array=[genexgenejgenek_f_a,genexgenejgenek_f_b]
   
    variable_b_array=[geneX_f_a*genej_f_a*genek_f_a+e_x1a*genek_f_a+e_x3a*genej_f_a+0.5*geneX_f_a*(e_genej_genek_a+e_genek_genej_a),
                    geneX_f_b*genej_f_b*genek_f_b+e_x1a*genek_f_b+e_x3b*genej_f_b+0.5*geneX_f_b*(e_genej_genek_b+e_genek_genej_b)] ## ignoring the digenic interaction terms with bem1 and bem3
    
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
    
    fc_list=(np.mean(variable_a_array)-np.mean(variable_b_array))
    gi[gene]["fold_change"]=fc_list

    gi[gene]["multiple-delete-fitness"]=np.mean(variable_a_array)
    gi[gene]["single-delete-fitness-product"]=np.mean(variable_b_array)
    
    gi[gene]["e_a"]=(variable_a_array[0]-variable_b_array[0])
    gi[gene]["e_b"]=(variable_a_array[1]-variable_b_array[1])
    gi[gene]["fold_change_std"]= np.std([gi[gene]["e_a"],gi[gene]["e_b"]])
    

# +
gi_pd_fitness_gene=pd.DataFrame.from_dict(gi,orient="index")

gi_pd_fitness_gene["gene_names"]=gi_pd_fitness_gene.index


gi_pd_fitness_gene_dbem1dbem3=gi_pd_fitness_gene

gi_pd_fitness_gene_dbem1dbem3

# +
all_pathways_gi=gi_pd_fitness_gene_dbem1dbem3[gi_pd_fitness_gene_dbem1dbem3.loc[:,"gene_names"].isin(all_pathways)]


#all_pathways_gi.to_csv("../postprocessed-data/pathways_bem1bem3_gi.csv",sep="\t")
# -

all_pathways_gi.sort_values(by="fold_change",ascending=False)

x=gi_pd_fitness_gene_dbem1dbem3.loc[:,"fold_change"].sort_values(ascending=False)
x[:10]

gi_pd_fitness_gene_bem3d[gi_pd_fitness_gene_bem3d.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)[0:10]

# +

bem1bem3_PI_sig,bem1bem3_NI_sig,bem1bem3_PI,bem1bem3_NI,bem1bem3_PI_all,bem1bem3_NI_all=classify_GI(gi_pd_fitness_gene_dbem1dbem3,col="fold_change")

for i in gi_pd_fitness_gene_dbem1dbem3.index:
    if i in bem1bem3_NI:
        gi_pd_fitness_gene_dbem1dbem3.loc[i,"significance4FC"]="neg"
    elif i in bem1bem3_PI:
        gi_pd_fitness_gene_dbem1dbem3.loc[i,"significance4FC"]="pos"
    else:
        gi_pd_fitness_gene_dbem1dbem3.loc[i,"significance4FC"]="none"

# +
with open("../postprocessed-data/positive_satay_genes_bem1_bem3.txt","w") as f:
    for i in bem1bem3_PI:
        f.write(i+"\n")
f.close()
with open("../postprocessed-data/negative_satay_genes_bem1_bem3.txt","w") as f:
    for i in bem1bem3_NI:
        f.write(i+"\n")
f.close()


# -

gi_pd_fitness_gene_dbem1dbem3[gi_pd_fitness_gene_dbem1dbem3.loc[:,"significance4FC"]=="pos"].sort_values(by="fold_change",ascending=False)

# +
gi_pd_fitness_gene_dbem1dbem3.loc[:,"e_a"].hist(alpha=0.2)
gi_pd_fitness_gene_dbem1dbem3.loc[:,"e_b"].hist(alpha=0.2)
gi_pd_fitness_gene_dbem1dbem3.loc[:,"fold_change"].hist(alpha=0.5)

plt.vlines(gi_pd_fitness_gene_dbem1dbem3.loc[:,"fold_change"].mean()+gi_pd_fitness_gene_dbem1dbem3.loc[:,"fold_change"].std(),
0,600,color="red")
plt.vlines(-gi_pd_fitness_gene_dbem1dbem3.loc[:,"fold_change"].mean()-gi_pd_fitness_gene_dbem1dbem3.loc[:,"fold_change"].std(),
0,600,color="red")
# -

gi_pd_fitness_gene_dbem1dbem3.sort_values(by="fold_change",ascending=False)

# +
# look for the fold change of the standard_essentials genes

gi_standard_essentials=defaultdict(dict)
gi_pd_fitness_gene=gi_pd_fitness_gene_dbem1dbem3
for i in standard_essentials:
    if i in gi_pd_fitness_gene.index:
        gi_standard_essentials[i]["fold_change"]=gi_pd_fitness_gene_dbem1dbem3.loc[i,"fold_change"]
        gi_standard_essentials[i]["p_statistic"]=gi_pd_fitness_gene_dbem1dbem3.loc[i,"p_statistic"]
        gi_standard_essentials[i]["significance"]=gi_pd_fitness_gene_dbem1dbem3.loc[i,"significance"]
        gi_standard_essentials[i]["significance4FC"]=gi_pd_fitness_gene_dbem1dbem3.loc[i,"significance4FC"]

gi_standard_essentials_pd=pd.DataFrame.from_dict(gi_standard_essentials,orient="index")


gi_standard_essentials_pd_dbem1dbem3=gi_standard_essentials_pd
#gi_standard_essentials_pd_dbem1dbem3[gi_standard_essentials_pd_dbem1dbem3.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)
gi_standard_essentials_pd_dbem1dbem3.sort_values(by="fold_change",ascending=False)

# +
## interaction score for polarity genes 

polarity_genes_gi=defaultdict(dict)

gi_pd_fitness_gene=gi_pd_fitness_gene_dbem1dbem3
for i in polarity_genes.index:
    if i in gi_pd_fitness_gene.index:
        polarity_genes_gi[i]["p_statistic"]=gi_pd_fitness_gene.loc[i,"p_statistic"]
        polarity_genes_gi[i]["significance"]=gi_pd_fitness_gene.loc[i,"significance"]
        polarity_genes_gi[i]["significance4FC"]=gi_pd_fitness_gene.loc[i,"significance4FC"]
        polarity_genes_gi[i]["fold_change"]=gi_pd_fitness_gene.loc[i,"fold_change"]
    
polarity_genes_gi_pd=pd.DataFrame.from_dict(polarity_genes_gi,orient="index")



polarity_genes_gi_pd_dbem1dbem3=polarity_genes_gi_pd

#len(polarity_genes_gi_pd)
polarity_genes_gi_pd_dbem1dbem3.sort_values(by="fold_change",ascending=False)

# +
from annotate_volcano import annotate_volcano

volcano_df=gi_pd_fitness_gene_dbem1dbem3

fig=annotate_volcano(volcano_df,figure_title="Interactors of dbem1dbem3 in WT")

plt.savefig("../figures/volcano_bem1bem3_wt.png",dpi=300)

# +
from annotate_volcano import annotate_volcano   #import annotate_volcano function
volcano_df=gi_pd_fitness_gene_dbem1dbem3
trackgene_list=["NRP1","CDC39"]
fig=annotate_volcano(volcano_df,figure_title="Interactors of dbem1dbem3 in WT",trackgene_list=trackgene_list)

plt.savefig("../figures/fig_volcano_interactors_bem1_bem3_with annotations.png",dpi=300,transparent=True)

# +
## Differences in essentiality scores from dbem1 to data_dbem1dbem3

plt.figure(figsize=(8,5))

colors={"pos":"#B55AC5","neg":"#5AAA46","none":"#BBBBBC"}

common_index=[]
for i in gi_pd_fitness_gene_dbem1dbem3.index:
    if i in gi_pd_fitness_gene_bem1d.index:
        common_index.append(i)
plt.scatter(gi_pd_fitness_gene_bem1d.loc[common_index,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[common_index,"fold_change"],
s=50,alpha=0.1,c="blue",label="all genes")

common_index=[]
for i in gi_standard_essentials_pd_dbem1dbem3.index:
    if i in gi_standard_essentials_pd_dbem1.index:
        common_index.append(i)
plt.scatter(gi_standard_essentials_pd_dbem1.loc[common_index,"fold_change"],gi_standard_essentials_pd_dbem1dbem3.loc[common_index,"fold_change"]
,s=50,c="black",alpha=0.6,label="standard essential genes")



plt.grid(linewidth=0.2)
plt.xlabel("dbem1 interaction scores",fontsize=14)
plt.ylabel("dbem1dbem3 interaction scores",fontsize=14)
plt.legend(fontsize=12)
# plt.hlines(0,0,10,linewidth=2,color="purple")
# plt.vlines(0,0,10,linewidth=2,color="purple")

# plt.hlines(0,-8,0,linewidth=2,color="green")
# plt.vlines(0,-8,0,linewidth=2,color="green")

plt.savefig("../figures/fig_scatter_dbem1_vs_dbem1dbem3_GI_Scores.png",dpi=300)

# +
common_index=[]
for i in gi_pd_fitness_gene_dbem1dbem3.index:
    if i in gi_pd_fitness_gene_bem1d.index:
        common_index.append(i)

gi_dbem1=gi_pd_fitness_gene_bem1d.loc[common_index,"fold_change"]
gi_dbem1dbem3=gi_pd_fitness_gene_dbem1dbem3.loc[common_index,"fold_change"]


gi_supression_shift=set(bem1_NI).intersection(bem1bem3_PI)

gi_essentiality_shift=set(bem1_PI).intersection(bem1bem3_NI)

len(gi_essentiality_shift)/len(gi_dbem1dbem3),len(gi_supression_shift)/len(gi_dbem1dbem3)
# -

len(gi_essentiality_shift),len(gi_supression_shift)

# +
## Differences in essentiality scores from dbem1 to data_dbem1dbem3

plt.figure(figsize=(8,5))

colors={"pos":"#B55AC5","neg":"#5AAA46","none":"#BBBBBC"}

common_index=[]
for i in gi_pd_fitness_gene_dbem1dbem3.index:
    if i in gi_pd_fitness_gene_bem1d.index:
        common_index.append(i)


plt.scatter(gi_pd_fitness_gene_bem1d.loc[common_index,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[common_index,"fold_change"],
s=50,alpha=0.1,c="gray")
plt.scatter(gi_pd_fitness_gene_bem1d.loc[gi_supression_shift,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[gi_supression_shift,"fold_change"],
s=50,alpha=1,c="gray",edgecolors="black",linewidth=2)
plt.scatter(gi_pd_fitness_gene_bem1d.loc[gi_essentiality_shift,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[gi_essentiality_shift,"fold_change"],
s=50,alpha=1,c="gray",edgecolors="black",linewidth=2)
plt.xlabel("dbem1 interaction scores",fontsize=14)
plt.ylabel("dbem1dbem3 interaction scores",fontsize=14)

plt.vlines(0,-8,4,color="black",linewidth=2,linestyles="dashed")
plt.hlines(0,-8,25,color="black",linewidth=2,linestyles="dashed")

plt.savefig("../figures/fig_scatter_dbem1_vs_dbem1dbem3_GI_Scores_with_shifts.png",dpi=300,transparent=True)

# +
common_index=[]
for i in gi_pd_fitness_gene_dbem1dbem3.index:
    if i in gi_pd_fitness_gene_bem3d.index:
        common_index.append(i)

gi_dbem3=gi_pd_fitness_gene_bem3d.loc[common_index,"fold_change"]
gi_dbem1dbem3=gi_pd_fitness_gene_dbem1dbem3.loc[common_index,"fold_change"]


gi_supression_shift=set(bem3_NI).intersection(bem1bem3_PI)

gi_essentiality_shift=set(bem3_PI).intersection(bem1bem3_NI)

len(gi_essentiality_shift),len(gi_supression_shift)
# -

gi_dbem1dbem3.loc[gi_essentiality_shift].sort_values(ascending=True)

# +
## Differences in essentiality scores from dbem1 to data_dbem1dbem3

plt.figure(figsize=(8,5))

colors={"pos":"#B55AC5","neg":"#5AAA46","none":"#BBBBBC"}



plt.scatter(gi_pd_fitness_gene_bem3d.loc[common_index,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[common_index,"fold_change"],
s=50,alpha=0.1,c="gray")
plt.scatter(gi_pd_fitness_gene_bem3d.loc[gi_supression_shift,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[gi_supression_shift,"fold_change"],
s=50,alpha=1,c="gray",edgecolors="black",linewidth=2)
plt.scatter(gi_pd_fitness_gene_bem3d.loc[gi_essentiality_shift,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[gi_essentiality_shift,"fold_change"],
s=50,alpha=1,c="gray",edgecolors="black",linewidth=2)
plt.xlabel("dbem3 interaction scores",fontsize=14)
plt.ylabel("dbem1dbem3 interaction scores",fontsize=14)

plt.vlines(0,-8,4,color="black",linewidth=2,linestyles="dashed")
plt.hlines(0,-1.5,1,color="black",linewidth=2,linestyles="dashed")

plt.savefig("../figures/fig_scatter_dbem3_vs_dbem1dbem3_GI_Scores_with_shifts.png",dpi=300,transparent=True)

# +
common_index=[]
for i in gi_pd_fitness_gene_bem3d.index:
    if i in gi_pd_fitness_gene_bem1d.index:
        common_index.append(i)

gi_dbem1=gi_pd_fitness_gene_bem1d.loc[common_index,"fold_change"]
gi_dbem3=gi_pd_fitness_gene_bem3d.loc[common_index,"fold_change"]


gi_supression_shift=set(bem1_NI).intersection(bem3_PI)

gi_essentiality_shift=set(bem1_PI).intersection(bem3_NI)

len(gi_essentiality_shift),len(gi_supression_shift)

# +
## Differences in essentiality scores from dbem1 to data_dbem1dbem3

plt.figure(figsize=(8,5))

colors={"pos":"#B55AC5","neg":"#5AAA46","none":"#BBBBBC"}



plt.scatter(gi_pd_fitness_gene_bem1d.loc[common_index,"fold_change"],gi_pd_fitness_gene_bem3d.loc[common_index,"fold_change"],
s=50,alpha=0.1,c="gray")
plt.scatter(gi_pd_fitness_gene_bem1d.loc[gi_supression_shift,"fold_change"],gi_pd_fitness_gene_bem3d.loc[gi_supression_shift,"fold_change"],
s=100,alpha=1,c="purple",edgecolors="green",linewidth=2)
plt.scatter(gi_pd_fitness_gene_bem1d.loc[gi_essentiality_shift,"fold_change"],gi_pd_fitness_gene_bem3d.loc[gi_essentiality_shift,"fold_change"],
s=100,alpha=1,c="green",edgecolors="blue",linewidth=2)
plt.xlabel("dbem1 interaction scores",fontsize=14)
plt.ylabel("dbem3 interaction scores",fontsize=14)

plt.vlines(0,-1.5,1,color="black",linewidth=2,linestyles="dashed")
plt.hlines(0,-5,25,color="black",linewidth=2,linestyles="dashed")

plt.savefig("../figures/fig_scatter_dbem3_vs_dbem1_GI_Scores_with_shifts.png",dpi=300)
# -

gi_dbem1.loc[gi_essentiality_shift].sort_values(ascending=True)

# +
## wich genes flip their interaction from dbem1 to data_dbem1dbem3

gene_pos_1_neg_13_sig=list(set(bem1_PI_sig).intersection(set(bem1bem3_NI_sig)))

gene_neg_1_pos_13_sig=list(set(bem1_NI_sig).intersection(set(bem1bem3_PI_sig)))

gene_pos_1_neg_13=list(set(bem1_PI).intersection(set(bem1bem3_NI)))

gene_neg_1_pos_13=list(set(bem1_NI).intersection(set(bem1bem3_PI)))

len(gene_pos_1_neg_13),len(gene_neg_1_pos_13),len(gene_pos_1_neg_13_sig),len(gene_neg_1_pos_13_sig)

# -

gi_pd_fitness_gene_bem1d.loc[gene_pos_1_neg_13_sig].sort_values(by="fold_change",ascending=False)

# +
# look for the fold change of the standard_essentials genes

gi_standard_essentials=defaultdict(dict)
gi_pd_fitness_gene=gi_pd_fitness_gene_bem3d
for i in standard_essentials:
    if i in gi_pd_fitness_gene.index:
        gi_standard_essentials[i]["fold_change"]=gi_pd_fitness_gene.loc[i,"fold_change"]
        gi_standard_essentials[i]["p_statistic"]=gi_pd_fitness_gene.loc[i,"p_statistic"]
        gi_standard_essentials[i]["significance"]=gi_pd_fitness_gene.loc[i,"significance"]
        gi_standard_essentials[i]["significance4FC"]=gi_pd_fitness_gene.loc[i,"significance4FC"]

gi_standard_essentials_pd=pd.DataFrame.from_dict(gi_standard_essentials,orient="index")


gi_standard_essentials_pd_bem3d=gi_standard_essentials_pd
gi_standard_essentials_pd_bem3d[gi_standard_essentials_pd_bem3d.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)

# +
## Changes across polarity proteins in different backgrounds 

p_b1=polarity_genes_gi_pd_dbem1.fold_change

p_b3=polarity_genes_gi_pd_dbem3.fold_change

p_b1b3=polarity_genes_gi_pd_dbem1dbem3.fold_change

### Build a dataframe with the fold changes of the polarity genes in the different backgrounds

p=pd.concat([p_b1,p_b3,p_b1b3],axis=1)

p.columns=["dbem1","dbem3","dbem1dbem3"]

p.dropna(inplace=True)


# +

import seaborn as sns
x=p


g=sns.clustermap(x,xticklabels=x.columns,cbar=True,annot=True,cmap="PuBuGn_r",vmax=x.max().max(),vmin=x.min().min(),figsize=(8,8))

g.savefig("../figures/fig_clustermap_polarity_genes_all_backgrounds.png",dpi=300)

# -

def data_safe2structure(data_safe):
   
    data_safe.drop("ORF",axis=1,inplace=True)
    data_safe.drop("Feature Qualifier",axis=1,inplace=True)
    data_safe.drop("Allele",axis=1,inplace=True)

    

    data_safe.reset_index(inplace=True)

    data_safe["Name"].replace("",np.nan,inplace=True)
    data_safe.dropna(inplace=True)

    data_safe.set_index("Name",inplace=True)
    # data_safe.to_csv("../postprocessed-data/tcm-safe-custom-230515-alltogether_analyze.csv")
    return data_safe

# +
## Heatmap of the SAFE scores for all backgrounds 

data_safe_ALL=pd.read_excel("../postprocessed-data/tcm-safe-custom-230515-alltogether.xlsx",sheet_name="SAFE enrichment scores",index_col=1)

data_safe_bem3PI=pd.read_excel("../postprocessed-data/tcm-safe-custom-230515-pi-bem3.xlsx",sheet_name="SAFE enrichment scores",index_col=1)

data_safe_bem3NI=pd.read_excel("../postprocessed-data/tcm-safe-custom-230515-ni-bem3.xlsx",sheet_name="SAFE enrichment scores",index_col=1)

data_safe_bem3_SGD=pd.read_excel("../postprocessed-data/tcm-safe-BEM3-230517_SGA_all.xlsx",sheet_name="SAFE enrichment scores",index_col=1)

data_safe_bem1PI=pd.read_excel("../postprocessed-data/tcm-safe-custom-230515-pi-bem1.xlsx",sheet_name="SAFE enrichment scores",index_col=1)

data_safe_bem1NI=pd.read_excel("../postprocessed-data/tcm-safe-custom-230515-ni-bem1.xlsx",sheet_name="SAFE enrichment scores",index_col=1)

data_safe_bem1_SGD=pd.read_excel("../postprocessed-data/tcm-safe-BEM1-230517_SGA_all.xlsx",sheet_name="SAFE enrichment scores",index_col=1)
data_safe_all_edited=data_safe2structure(data_safe_ALL)
data_safe_bem3PI_edited=data_safe2structure(data_safe_bem3PI)
data_safe_bem3NI_edited=data_safe2structure(data_safe_bem3NI)

data_safe_bem3_SGD_edited=data_safe2structure(data_safe_bem3_SGD)

data_safe_bem1PI_edited=data_safe2structure(data_safe_bem1PI)
data_safe_bem1NI_edited=data_safe2structure(data_safe_bem1NI)

data_safe_bem1_SGD_edited=data_safe2structure(data_safe_bem1_SGD)



# +

data_safe_bem3=pd.concat([data_safe_bem3NI_edited,data_safe_bem3PI_edited],axis=0)
data_safe_bem3.fillna(0,inplace=True)

data_safe_bem1=pd.concat([data_safe_bem1NI_edited,data_safe_bem1PI_edited],axis=0)
data_safe_bem1.fillna(0,inplace=True)


# +
from matplotlib.patches import Patch
def datasafe2heatmap(data_edited,background,savefig=False):

    data_safe4heatmap=data_edited.iloc[:,1:]
    data_safe4heatmap=data_safe4heatmap.astype(float)
    
    data_safe4heatmap.dropna(inplace=True)

### pLOTTING THE HEATMAP

    x=data_safe4heatmap

    labels = data_edited.Annotations
    colors=sns.color_palette("colorblind", len(set(labels)))
    lut = dict(zip(set(labels), colors))
    row_colors = pd.DataFrame(labels)["Annotations"].map(lut)

    g=sns.clustermap(x,cmap=sns.color_palette("blend:#7AB,#EDA"),vmax=x.max()[1:].max(),vmin=x.min()[1:].min(),
    figsize=(15,12),row_colors=row_colors,col_cluster=False,row_cluster=False)



    handles = [Patch(facecolor=lut[name]) for name in lut]
    plt.legend(handles, lut, title='Annotations', loc='upper left', bbox_to_anchor=(2.5, 0.9),fontsize=14,ncol=3)

    if savefig:
        plt.tight_layout()
        plt.savefig("../figures/fig_clustermap_safe_"+ background + "backgrounds.png",dpi=300)

    return g
# -

g=datasafe2heatmap(data_safe_bem1_SGD_edited,"BEM1_SGD",savefig=True)

# +
## Common interactors across backgrounds 

common_negative=set(bem1bem3_NI)&set(bem1_NI)&set(bem3_NI)
common_positive=set(bem1bem3_PI)&set(bem1_PI)&set(bem3_PI)

common_positive_d13_d1=set(bem1bem3_PI)&set(bem1_PI)
common_positive_d13_d3=set(bem1bem3_PI)&set(bem3_PI)

common_negative_d13_d1=set(bem1bem3_NI)&set(bem1_NI)
common_negative_d13_d3=set(bem1bem3_PI)&set(bem1_PI)

len(common_negative),len(common_positive),len(common_positive_d13_d1),len(common_positive_d13_d3),len(common_negative_d13_d1),len(common_negative_d13_d3)

# +
## Common essential genes across backgrounds
essentials_dbem3_pos=gi_standard_essentials_pd_bem3d[gi_standard_essentials_pd_bem3d.loc[:,"significance"]==True & (gi_standard_essentials_pd_bem3d.loc[:,"significance4FC"]=="pos")].index

essentials_dbem1_pos=gi_standard_essentials_pd_dbem1[gi_standard_essentials_pd_dbem1.loc[:,"significance"]==True & (gi_standard_essentials_pd_dbem1.loc[:,"significance4FC"]=="pos")].index

essentials_dbem1dbem3_pos=gi_standard_essentials_pd_dbem1dbem3[gi_standard_essentials_pd_dbem1dbem3.loc[:,"significance"]==True & (gi_standard_essentials_pd_dbem1dbem3.loc[:,"significance4FC"]=="pos")].index

common_essentials_pos=set(essentials_dbem3_pos)&set(essentials_dbem1_pos)&set(essentials_dbem1dbem3_pos)

common_essentials_pos_db13_db1=set(essentials_dbem1dbem3_pos)&set(essentials_dbem1_pos)
common_essentials_pos_db13_db3=set(essentials_dbem1dbem3_pos)&set(essentials_dbem3_pos)
# -

len(common_essentials_pos),len(common_essentials_pos_db13_db1),len(common_essentials_pos_db13_db3)

# ### Fitness for detecting essentiality, making zero regions that do not have enough reads , insertions or flanking regions

data_fitness=filter_fitness(fitness_all_pd,backgrounds=keys,goi=["BEM1","BEM3","NRP1"],discard=["Not enough flanking regions"],set2zero=["Not enough reads",
    "Not enough insertions"],cols=["fitness_gene","fitness_domains_corrected"],essentiality=True)

# +
## take the same index across all libraries

index_arrays=np.array([data_fitness.loc["wt_merged"].index,data_fitness.loc["bem1-aid_merged"].index,
data_fitness.loc["bem1-aid-dbem3_a"].index,data_fitness.loc["bem1-aid-dbem3_b"].index,data_fitness.loc["dbem3_merged"].index],
dtype="object")

d1=set.intersection(*map(set,index_arrays))

d1_index=np.array(list(d1))

# index_BEM1=np.where(d1_index=="BEM1")[0][0]
# index_BEM3=np.where(d1_index=="BEM3")[0][0]
index_NRP1=np.where(d1_index=="NRP1")[0][0]


data_wt=data_fitness.loc["wt_merged"].loc[d1,"fitness_gene"]
data_dbem1=data_fitness.loc["bem1-aid_merged"].loc[d1,"fitness_gene"]
data_dbem3=data_fitness.loc["dbem3_merged"].loc[d1,"fitness_gene"]
data_dbem1dbem3=np.mean(np.asanyarray([data_fitness.loc["bem1-aid-dbem3_a"].loc[d1,"fitness_gene"],
data_fitness.loc["bem1-aid-dbem3_b"].loc[d1,"fitness_gene"]],dtype="object"),axis=0)

fitness_bem1_wt=(data_fitness.loc["wt_merged"].loc["BEM1","fitness_gene"]+data_fitness.loc["wt_merged"].loc["BEM1","fitness_domains_corrected"])/2

fitness_bem3_wt=data_wt.loc["BEM3"]

## Normalize the datasetss such as the median values corresponds to the fitness of the knockout of the gene of interest 
data_wt_norm=data_wt

data_dbem1_norm=data_dbem1*fitness_bem1_wt/np.median(data_dbem1)
data_dbem3_norm=data_dbem3*fitness_bem3_wt/np.median(data_dbem3)

fitness_bem3_dbem1=data_dbem1_norm.loc["BEM3"]

data_dbem1dbem3_norm=data_dbem1dbem3*fitness_bem3_dbem1/np.median(data_dbem1dbem3)


# +
## Predicting number of essential genes based on the corrected fitness values

## Essential genes

# in WT threshold f<0.4 (from domains corrected)

data_fitness_wta=data_fitness.loc["wt_a","fitness_domains_corrected"]
data_fitness_wtb=data_fitness.loc["wt_b","fitness_domains_corrected"]



fitness_dbem3a=data_fitness.loc["wt_a","fitness_gene"].loc["BEM3"]
fitness_dbem3b=data_fitness.loc["wt_b","fitness_gene"].loc["BEM3"]

data_fitness_dbem1a=data_fitness.loc["bem1-aid_a","fitness_domains_corrected"]*fitness_bem1_wt
data_fitness_dbem1b=data_fitness.loc["bem1-aid_b","fitness_domains_corrected"]*fitness_bem1_wt

fitness_bem3d_dbem1a=data_fitness.loc["bem1-aid_a","fitness_gene"].loc["BEM3"]
fitness_bem3d_dbem1b=data_fitness.loc["bem1-aid_b","fitness_gene"].loc["BEM3"]

data_fitness_dbem1dbem3a=data_fitness.loc["dbem1dbem3_a","fitness_domains_corrected"]*fitness_bem3_dbem1
data_fitness_dbem1dbem3b=data_fitness.loc["dbem1dbem3_b","fitness_domains_corrected"]*fitness_bem3_dbem1

data_fitness_dbem3a=data_fitness.loc["dbem3_a","fitness_domains_corrected"]*fitness_bem3_wt
data_fitness_dbem3b=data_fitness.loc["dbem3_b","fitness_domains_corrected"]*fitness_bem3_wt
# in WT threshold f<0.4(mean)-0.28(std) (from domains corrected)
essential_genes_wta=data_fitness_wta[data_fitness_wta<0.12].index
essential_genes_wtb=data_fitness_wtb[data_fitness_wtb<0.12].index

essential_genes_dbem1a=data_fitness_dbem1a[data_fitness_dbem1a<0.12*fitness_bem1_wt].index
essential_genes_dbem1b=data_fitness_dbem1b[data_fitness_dbem1b<0.12*fitness_bem1_wt].index

essential_genes_dbem1dbem3a=data_fitness_dbem1dbem3a[data_fitness_dbem1dbem3a<0.12*fitness_bem3_dbem1].index
essential_genes_dbem1dbem3b=data_fitness_dbem1dbem3b[data_fitness_dbem1dbem3b<0.12*fitness_bem3_dbem1].index

essential_genes_dbem3a=data_fitness_dbem3a[data_fitness_dbem3a<0.12*fitness_bem3_wt].index
essential_genes_dbem3b=data_fitness_dbem3b[data_fitness_dbem3b<0.12*fitness_bem3_wt].index

essential_genes_wt=set(essential_genes_wta)&set(essential_genes_wtb)
essential_genes_dbem1=set(essential_genes_dbem1a)&set(essential_genes_dbem1b)
essential_genes_dbem1dbem3=set(essential_genes_dbem1dbem3a)&set(essential_genes_dbem1dbem3b)
essential_genes_dbem3=set(essential_genes_dbem3a)&set(essential_genes_dbem3b)

# +
number_e_genes=[len(essential_genes_wt),len(essential_genes_dbem1),len(essential_genes_dbem1dbem3),len(essential_genes_dbem3)]

std_number_e_genes=[np.std([len(essential_genes_wta),len(essential_genes_wtb)]),
np.std([len(essential_genes_dbem1a),len(essential_genes_dbem1b)]),
np.std([len(essential_genes_dbem1dbem3a),len(essential_genes_dbem1dbem3b)]),
np.std([len(essential_genes_dbem3a),len(essential_genes_dbem3b)])]


# +
# e_13_pos_1=[]
# for i in essential_genes_dbem1dbem3:
#     if i in gene_pos_1_neg_13:
#         e_13_pos_1.append(i)

# e_13_pos_1

# +
plt.figure(figsize=(8,5))

number_e_genes_refactor=[number_e_genes[1],number_e_genes[3],number_e_genes[2],number_e_genes[0],number_e_genes[1]+number_e_genes[3]]
std_number_e_genes_refactor=[std_number_e_genes[1],std_number_e_genes[3],std_number_e_genes[2],std_number_e_genes[0],std_number_e_genes[1]+std_number_e_genes[3]]
plt.errorbar(["bem1$\Delta$","bem3$\Delta$","bem1$\Delta$bem3$\Delta$","WT","Expected \n in dbem1dbem3"],number_e_genes_refactor,
yerr=std_number_e_genes_refactor,fmt="o",color="black",capsize=10)
plt.ylabel("Predicted number of essential genes")

plt.plot(number_e_genes_refactor,"--",color="black",linewidth=0.5)
plt.yscale("log")

plt.tight_layout()
#plt.savefig("../figures/predicted_number_essential_genes_change_trajectory.png",dpi=300)

# +
fig, ax1 = plt.subplots(figsize=(8,5))
number_e_genes_refactor=[number_e_genes[0],number_e_genes[1],number_e_genes[3],number_e_genes[2],number_e_genes[1]+number_e_genes[3]]
std_number_e_genes_refactor=[std_number_e_genes[0],std_number_e_genes[1],std_number_e_genes[3],std_number_e_genes[2],std_number_e_genes[1]+std_number_e_genes[3]]
ax2 = ax1.twinx()
labels=["WT","$\Delta$bem1","$\Delta$bem3","$\Delta$bem1$\Delta$bem3","Expected"]
ax1.errorbar(labels,number_e_genes_refactor,
yerr=std_number_e_genes_refactor,fmt="x",color="black",capsize=10,markersize=10)
ax1.set_ylabel("Predicted number of essential genes",color="black")


growth_rate_steps=[0.01244,0.0022,0.01161,0.010,0.0022*0.01161]

growth_rate_steps_std=[0.00028,0.00015,0.00073,0.0009]

std_expected=uncertainty_propagation(growth_rate_steps[-1],growth_rate_steps[1],
growth_rate_steps[2],growth_rate_steps_std[1],growth_rate_steps_std[2])

growth_rate_steps_std.append(std_expected)





ax1.plot(number_e_genes_refactor,"--",color="black",linewidth=0.5)

ax2.errorbar(["WT","$\Delta$bem1","$\Delta$bem3","$\Delta$bem1$\Delta$bem3","Expected"],growth_rate_steps,yerr=growth_rate_steps_std,fmt="o",
color="gray",capsize=10,markersize=10)
ax2.plot(growth_rate_steps,"--",color="black",linewidth=0.5)

ax2.set_ylabel("Growth rate (min$^{-1}$)",color="gray")

ax1.tick_params(axis='y', labelcolor="black")
ax2.tick_params(axis='y', labelcolor="gray")

plt.tight_layout()

plt.savefig("../figures/predicted_number_essential_genes_trajectory_growth_rate.png",dpi=300,transparent=True) 
# #plt.savefig("../figures/predicted_number_essential_genes_trajectory.png",dpi=300,transparent=True)

# +
## Scatter plot predcited number of essential genes vs growth rate

plt.figure(figsize=(8,5))

growth_rate_steps=[0.0022,0.010,0.01244] # bem1d,bem1dbem3d,wt
growth_rate_steps_std=[0.00015,0.0009,0.00028]

number_e_genes=[len(essential_genes_dbem1),len(essential_genes_dbem1dbem3),len(essential_genes_wt)]

std_number_e_genes=[np.std([len(essential_genes_dbem1a),len(essential_genes_dbem1b)]),
np.std([len(essential_genes_dbem1dbem3a),len(essential_genes_dbem1dbem3b)]),
np.std([len(essential_genes_wta),len(essential_genes_wtb)])]

plt.errorbar(growth_rate_steps,number_e_genes,yerr=std_number_e_genes,xerr=growth_rate_steps_std,fmt="o",color="black",capsize=5)
plt.plot(growth_rate_steps,number_e_genes,"--",color="black",linewidth=0.5)
# plt.annotate("WT",xy=(growth_rate_steps[2],number_e_genes[2]),xytext=(growth_rate_steps[2]-0.0004,number_e_genes[2]+200),fontsize=14)
# plt.annotate("$\Delta$bem1",xy=(growth_rate_steps[0],number_e_genes[0]),xytext=(growth_rate_steps[0]-0.0005,number_e_genes[0]+500),fontsize=14)
# plt.annotate("$\Delta$bem1$\Delta$bem3",xy=(growth_rate_steps[1],number_e_genes[1]),xytext=(growth_rate_steps[1]-0.0008,
# number_e_genes[1]+400),fontsize=14)

plt.ylim(500,5000)
plt.yticks([1000,2000,3000,4000,5000],fontsize=14)
plt.xticks([0.0022,0.010,0.01244],fontsize=14)
plt.grid(linewidth=0.3,linestyle="dashed")
plt.xlabel("Growth rate (min$^{-1}$)")
plt.ylabel("Predicted number of essential genes")

plt.savefig("../figures/fig_scatter_predicted_number_essential_genes_vs_growth_rate.png",dpi=300)

# +
## How many of the predicted genes are also WT essential gi_standard_essentials

essential_wt_consistent=set(essential_genes_wt)&set(standard_essentials)

essential_dbem1_consistent=set(essential_genes_dbem1)&set(standard_essentials)

essential_dbem1dbem3_consistent=set(essential_genes_dbem1dbem3)&set(standard_essentials)

len(essential_wt_consistent),len(essential_dbem1_consistent),len(essential_dbem1dbem3_consistent)

# +
## Plot how many essential genes are common between backgrounds 

common_wt_bem1_e_genes=set(essential_genes_wt)&set(essential_genes_dbem1)
common_wt_bem1bem3_e_genes=set(essential_genes_wt)&set(essential_genes_dbem1dbem3)
common_bem1_bem1bem3_e_genes=set(essential_genes_dbem1)&set(essential_genes_dbem1dbem3)

common_wt_bem1_bem1bem3_e_genes=set(essential_genes_wt)&set(essential_genes_dbem1)&set(essential_genes_dbem1dbem3)

a=[len(essential_genes_wt), len(essential_genes_dbem1), len(common_wt_bem1_e_genes) ,
len(essential_genes_dbem1dbem3), len(common_wt_bem1bem3_e_genes), len(common_bem1_bem1bem3_e_genes), len(common_wt_bem1_bem1bem3_e_genes)]

from matplotlib_venn import venn3, venn3_circles
# a=(np.round(number_e_genes[0],0).astype(int), np.round(number_e_genes[1],0).astype(int), len(common_wt_bem1_e_genes), number_e_genes[2], 
# len(common_wt_bem1bem3_e_genes), len(common_bem1_bem1bem3_e_genes), len(common_wt_bem1_bem1bem3_e_genes))
venn3(subsets = a, set_labels = ('Essentials WT', 'Essentials $\Delta$bem1','Essentials $\Delta$bem1$\Delta$bem3'), alpha = 0.5,
set_colors=('gray', 'blue', 'yellow'));

plt.savefig("../figures/venn3_common_essentials_dbem1_wt_dbem1dbem3.png", dpi=300, bbox_inches="tight")

# +
## Plot how many of the essentials in different backgrounds are also conserved with existing WT essential genes
a1=(len(essential_wt_consistent)/len(standard_essentials))
a2=(len(essential_dbem1_consistent)/len(standard_essentials))
a4=(len(essential_dbem1dbem3_consistent)/len(standard_essentials))
a3=(len(set(essential_wt_consistent)&set(essential_dbem1_consistent))/(len(essential_wt_consistent)+len(essential_dbem1_consistent)))
a5=(len(set(essential_wt_consistent)&set(essential_dbem1dbem3_consistent))/(len(essential_wt_consistent)+len(essential_dbem1dbem3_consistent)))
a6=(len(set(essential_wt_consistent)&set(essential_dbem1dbem3_consistent))/(len(essential_dbem1_consistent)+len(essential_dbem1dbem3_consistent)))
a7=(len(set(essential_wt_consistent)&set(essential_dbem1_consistent)&set(essential_dbem1dbem3_consistent))/(len(essential_wt_consistent)+len(essential_dbem1_consistent)+len(essential_dbem1dbem3_consistent)))

a=[a1,a2,a3,a4,a5,a6,a7]

venn3(subsets = a, set_labels = ('Essentials WT', 'Essentials $\Delta$bem1','Essentials $\Delta$bem1$\Delta$bem3'), alpha = 0.5,
set_colors=('gray', 'blue', 'yellow'),subset_label_formatter=lambda x: f"{(x):1.0%}");

plt.savefig("../figures/venn3_common_essentials_dbem1_wt_dbem1dbem3_compared2standard_essentials.png", dpi=300, bbox_inches="tight")

# +

data=[data_wt_norm,data_dbem1_norm,data_dbem1dbem3_norm,data_dbem3_norm]
labels=["WT", "$\Delta$bem1", "$\Delta$bem1$\Delta$bem3", "$\Delta$bem3"]
colors=["gray", "blue", "yellow", "green"]
bounds=[0.12,0.12*fitness_bem1_wt,0.12*fitness_bem3_dbem1,0.12*fitness_bem3_wt]
plt.subplots(2,2,figsize=(10,10))
plt.subplots_adjust(hspace=0.5,wspace=0.5)

for i in range(4):
    plt.subplot(2,2,i+1)
    plt.hist(data[i], bins=20, alpha=0.5, label=labels[i],color=colors[i])
    plt.vlines(bounds[i],0,2000,linestyles="dashed",linewidth=1,color="k")
    plt.xlabel("Fitness")
    plt.ylabel("Frequency")
    plt.legend()


plt.savefig("../figures/fitness_distribution_dbem1_wt_dbem1dbem3_dbem3.png", dpi=300, bbox_inches="tight")

# +
import gseapy as gp
from gseapy import barplot, dotplot
type_gi="Positive_GI_SGA"
goi=list(pos_bem1_genes)
yeast = gp.get_library_name(organism='Yeast')

sets=['GO_Biological_Process_AutoRIF','GO_Cellular_Component_AutoRIF',
'GO_Molecular_Function_AutoRIF','Pfam_Domains_2019','Phenotype_AutoRIF' ] 

# #%% enrichment 

enr_data=[]

for i in np.arange(0,len(sets)): 

  enr=gp.enrichr(gene_list=goi,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                
                  outdir='../postprocessed-data/enrich-analysis_dbem1/'+type_gi+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
  enr_data.append(enr.res2d)
# -

pos_bem1_genes

# +
# #%% barplot

fig,ax=plt.subplots(1,2,figsize=(15, 5))

plt.subplots_adjust(wspace=1.2, hspace=0.5)
sets_reset=[0,4]


for i,j in zip(sets_reset,np.arange(0,len(sets_reset))):
    data=enr_data[i]
    data=data[data.loc[:,"Adjusted P-value"]!=0]
    tmp=data.sort_values(by="Adjusted P-value",ascending=True)[0:5]
    
    tmp_axis_x=[]
    for k in np.arange(0,len(tmp)):
        tmp_axis_x.append(tmp.Term.tolist()[k].replace("_"," "))

    ax[j].barh(tmp_axis_x,tmp.loc[:,"Adjusted P-value"],color="purple",alpha=0.5)
    ax[j].set_title(sets[i].replace("_"," "),fontsize=18)
    ax[j].set_xlabel("Adjusted P-value",fontsize=18)
    ax[j].tick_params(axis='both', labelsize=18)


    # ax[j//2,j%2].set_xlim(0,0.025)
    
plt.tight_layout()

plt.savefig("../figures/positive_BEM1_gi_sga_enrichment.png", dpi=300, bbox_inches="tight")

# +
import gseapy as gp
from gseapy import barplot, dotplot
type_gi="Positive_GI_SATAY"
goi=list(bem1_PI)
yeast = gp.get_library_name(organism='Yeast')

sets=['GO_Biological_Process_AutoRIF','GO_Cellular_Component_AutoRIF',
'GO_Molecular_Function_AutoRIF','Pfam_Domains_2019','Phenotype_AutoRIF',"WikiPathways_2018" ] 


# #%% enrichment 

enr_data=[]

for i in np.arange(0,len(sets)): 

  enr=gp.enrichr(gene_list=goi,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                
                  outdir='../postprocessed-data/enrich-analysis_dbem1/'+type_gi+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
  enr_data.append(enr.res2d)

# +
# #%% barplot

fig,ax=plt.subplots(1,2,figsize=(20, 5))

plt.subplots_adjust(wspace=1.2, hspace=0.5)
sets_reset=[0,5]


for i,j in zip(sets_reset,np.arange(0,len(sets_reset))):
    data=enr_data[i]
    data=data[data.loc[:,"Adjusted P-value"]!=0]
    tmp=data.sort_values(by="Adjusted P-value",ascending=True)[0:5]
    
    tmp_axis_x=[]
    for k in np.arange(0,len(tmp)):
        tmp_axis_x.append(tmp.Term.tolist()[k].replace("_"," "))

    ax[j].barh(tmp_axis_x,tmp.loc[:,"Adjusted P-value"],color="purple",alpha=0.5)
    ax[j].set_title(sets[i].replace("_"," "),fontsize=18)
    ax[j].set_xlabel("Adjusted P-value",fontsize=18)
    ax[j].tick_params(axis='both', labelsize=18)


    # ax[j//2,j%2].set_xlim(0,0.025)
    
plt.tight_layout()

#plt.savefig("../figures/positive_BEM1_gi_satay_enrichment.png", dpi=300, bbox_inches="tight")

# +
import gseapy as gp
from gseapy import barplot, dotplot
type_gi="Negative_GI_SGA"
goi=list(neg_sl_bem1_genes)
yeast = gp.get_library_name(organism='Yeast')


# #%% enrichment 

enr_data=[]

for i in np.arange(0,len(sets)): 

  enr=gp.enrichr(gene_list=goi,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                
                  outdir='../postprocessed-data/enrich-analysis_dbem1/'+type_gi+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
  enr_data.append(enr.res2d)

# +
# #%% barplot

fig,ax=plt.subplots(1,2,figsize=(15, 5))

plt.subplots_adjust(wspace=1.2, hspace=0.5)
sets_reset=[0,5]


for i,j in zip(sets_reset,np.arange(0,len(sets_reset))):
    data=enr_data[i]
    data=data[data.loc[:,"Adjusted P-value"]!=0]
    tmp=data.sort_values(by="Adjusted P-value",ascending=True)[0:5]
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

plt.savefig("../figures/negative_BEM1_gi_sga_enrichment.png", dpi=300, bbox_inches="tight")

# +
import gseapy as gp
from gseapy import barplot, dotplot
type_gi="Negative_GI_SATAY"
goi=list(bem1_NI_sig)
yeast = gp.get_library_name(organism='Yeast')


# #%% enrichment 

enr_data=[]

for i in np.arange(0,len(sets)): 

  enr=gp.enrichr(gene_list=goi,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                
                  outdir='../postprocessed-data/enrich-analysis_dbem1/'+type_gi+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
  enr_data.append(enr.res2d)

# +
# #%% barplot

fig,ax=plt.subplots(1,2,figsize=(15, 5))

plt.subplots_adjust(wspace=1.2, hspace=0.5)
sets_reset=[0,5]


for i,j in zip(sets_reset,np.arange(0,len(sets_reset))):
    data=enr_data[i]
    data=data[data.loc[:,"Adjusted P-value"]!=0]
    tmp=data.sort_values(by="Adjusted P-value",ascending=True)[0:5]
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

plt.savefig("../figures/negative_BEM1_gi_satay_enrichment.png", dpi=300, bbox_inches="tight")

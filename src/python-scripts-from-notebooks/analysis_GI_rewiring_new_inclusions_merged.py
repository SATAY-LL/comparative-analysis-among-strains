# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: transposonmapper
#     language: python
#     name: python3
# ---

# +
import scipy
from scipy import stats
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import pickle
import seaborn as sns
from collections import defaultdict
### Genetic interactions 

from functions_interaction_computations import digenic_GI
from functions_interaction_computations import classify_GI

# +
with open("../postprocessed-data/fitness_models_all_backgrounds_enzo_analysis_merged", "rb") as fp:   # Unpickling
    fitness_average = pickle.load(fp)

with open("../postprocessed-data/fitness_models_all_backgrounds_domain_analysis_merged", "rb") as fp:   # Unpickling
    fitness_domains = pickle.load(fp)

with open("../postprocessed-data/fitness_models_all_backgrounds_domain_analysis_corrected", "rb") as fp:   # Unpickling
    fitness_domains_corrected = pickle.load(fp)
standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str) #standard essentials
data_domains=pd.read_excel("../postprocessed-data/genomic-domains-wt.xlsx",index_col="Unnamed: 0")
polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")
polarity_genes.fillna(0,inplace=True)
# -

data_fitness=fitness_average.copy()
data_fitness.columns=['mean', 'sem', 'fitness_gene']

# +
## Boxplots of normalized fitness 

data_fitness_wt=data_fitness.loc["wt_merged","fitness_gene"]

data_fitness_dbem3=data_fitness.loc["dbem3_merged","fitness_gene"]*data_fitness_wt.loc["BEM3"]

data_fitness_dbem1=data_fitness.loc["bem1-aid_merged","fitness_gene"]*data_fitness_wt.loc["BEM1"]


data_fitness_dbem1dbem3=np.mean(pd.merge(data_fitness.loc["dbem1dbem3_a","fitness_gene"],
data_fitness.loc["dbem1dbem3_b","fitness_gene"],left_index=True,right_index=True),axis=1)*data_fitness_dbem1.loc["BEM3"]


## remove nan values from all dataframes 

data_fitness_wt=data_fitness_wt.dropna()
data_fitness_dbem1=data_fitness_dbem1.dropna()
data_fitness_dbem1dbem3=data_fitness_dbem1dbem3.dropna()


# Common indexes among all dataframes

common_index=data_fitness_wt.index.intersection(data_fitness_dbem1.index).intersection(data_fitness_dbem1dbem3.index)

# -

data_fitness_dbem1.loc["BEM3"]

data_fitness_wt.loc["BEM1"]

# +
error_bars_wt=np.std(data_fitness_wt.loc[common_index])/np.sqrt(len(data_fitness_wt.loc[common_index]))

error_bars_dbem1=fitness_average.loc["wt_merged","sem"].loc["BEM1"]

error_bars_dbem1dbem3=fitness_average.loc["bem1-aid_merged","sem"].loc["BEM3"]

error_bars_dbem1dbem3dnrp1=np.mean([fitness_average.loc["dbem1dbem3_a","sem"].loc["NRP1"],
fitness_average.loc["dbem1dbem3_b","sem"].loc["NRP1"]])


# +
fig,ax=plt.subplots(figsize=(5,3))

sns.kdeplot(data_fitness_wt.loc[common_index],ax=ax,color="gray",label="WT",multiple="stack"
            ,fill=True,common_norm=False,alpha=0.5)
sns.kdeplot(data_fitness_dbem1.loc[common_index],ax=ax,color="blue",label="dbem1",multiple="stack",
            fill=True,common_norm=False,alpha=0.5)
sns.kdeplot(data_fitness_dbem1dbem3.loc[common_index],ax=ax,color="orange",label="dbem1dbem3",multiple="stack",
            fill=True,common_norm=False,alpha=0.5)

ax.set_xlabel("Fitness")

#fig.savefig("../figures/fitness_distributions_trajectory.png",dpi=300,bbox_inches="tight")

# +
fig,ax=plt.subplots(figsize=(5,5))

box=ax.boxplot([data_fitness_wt.loc[common_index],
data_fitness_dbem1.loc[common_index],data_fitness_dbem1dbem3.loc[common_index]],labels=["WT","dbem1","dbem1dbem3"],
showfliers=False,patch_artist=True,notch=True,boxprops=dict(alpha=.5),whiskerprops=dict(alpha=.7),capprops=dict(alpha=.7),);



colors = ['black', "black", '#D0D3D6']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)


labels2scatter=["WT","$\Delta$bem1","$\Delta$bem1$\Delta$bem3","$\Delta$bem1$\Delta$bem3$\Delta$nrp1"]

growth_rate_steps_satay=[np.nanmedian(data_fitness_wt.loc[common_index]),
data_fitness_wt.loc["BEM1"],data_fitness_dbem1.loc["BEM3"],data_fitness_dbem1dbem3.loc["NRP1"]]

growth_rate_steps_std=[error_bars_wt,error_bars_dbem1,error_bars_dbem1dbem3,error_bars_dbem1dbem3dnrp1]

ax.scatter(labels2scatter,growth_rate_steps_satay,color=["#464546","#464546","#464546","#D0D3D6"],s=50)
ax.errorbar(labels2scatter,growth_rate_steps_satay,yerr=growth_rate_steps_std,linestyle="None",
color="#464546",alpha=0.6,capsize=5)

ax.set_ylabel("Fitness")

fig.savefig("../figures/fitness_boxplot_trajectory.png",dpi=300,bbox_inches="tight")

# -

gi_bem1=digenic_GI(data=data_fitness,goi="BEM1",col_fitness="fitness_gene",backg=["wt_a","wt_b","bem1-aid_a","bem1-aid_b"])
gi_bem3=digenic_GI(data=data_fitness,goi="BEM3",col_fitness="fitness_gene",backg=["wt_a","wt_b","dbem3_a","dbem3_b"])

# +
gi_pd_fitness_gene=pd.DataFrame.from_dict(gi_bem1,orient="index")

gi_pd_fitness_gene["gene_names"]=gi_pd_fitness_gene.index

gi_pd_fitness_gene_bem1d=gi_pd_fitness_gene
gi_pd_fitness_gene_bem1d.dropna(inplace=True)
gi_pd_fitness_gene_bem1d.sort_values(by="fold_change",ascending=False)
# -

gi_pd_fitness_gene_bem1d[gi_pd_fitness_gene_bem1d.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)

# +
gi_pd_fitness_gene=pd.DataFrame.from_dict(gi_bem3,orient="index")

gi_pd_fitness_gene["gene_names"]=gi_pd_fitness_gene.index

gi_pd_fitness_gene_bem3d=gi_pd_fitness_gene
gi_pd_fitness_gene_bem3d.dropna(inplace=True)
gi_pd_fitness_gene_bem3d.sort_values(by="fold_change",ascending=False)

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

all_pathways=cellcycle+mapk+lipid+glycolisis+fatty_acid

all_pathways=list(set(all_pathways))

# make all the gene names uppercase

all_pathways=[i.upper() for i in all_pathways]

all_pathways_gi=gi_pd_fitness_gene[gi_pd_fitness_gene.loc[:,"gene_names"].isin(all_pathways)]

all_pathways_gi.to_csv("../postprocessed-data/pathways_bem1_gi.csv",sep="\t")

all_pathways_gi=gi_pd_fitness_gene_bem3d[gi_pd_fitness_gene_bem3d.loc[:,"gene_names"].isin(all_pathways)]

all_pathways_gi.to_csv("../postprocessed-data/pathways_bem3_gi.csv",sep="\t")

all_pathways_gi.sort_values(by="fold_change",ascending=False)
# -

np.abs(np.mean(gi_pd_fitness_gene_bem1d.fold_change))+np.std(gi_pd_fitness_gene_bem1d.fold_change)

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

len(bem1_NI2export),len(bem1_PI2export)

# +
gi_pd_fitness_gene_bem1d.loc[:,"fold_change"].hist(bins=50)

plt.vlines(gi_pd_fitness_gene_bem1d.loc[:,"fold_change"].median(),0,650,color="red")

# +
## export to a txt the gene names of the significant positive and negative regulators 

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
from annotate_volcano import annotate_volcano

volcano_df=gi_pd_fitness_gene_bem1d

trackgene_list=["BEM3","NRP1"]

fig=annotate_volcano(volcano_df,figure_title="Interactors of bem1 in WT",trackgene_list=trackgene_list)

#fig.savefig("../figures/volcano_bem1.png",dpi=300,bbox_inches="tight")

# +
volcano_df=gi_pd_fitness_gene_bem3d

fig=annotate_volcano(volcano_df,figure_title="Interactors of bem3 in WT")


# +
## interaction score for polarity genes 
from collections import defaultdict
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
## Existing data from the literature about BEM1 interactors from Costanzo et al. 2010,2016 and Tong et al. 2004

existing_data=pd.read_excel("../postprocessed-data/GI_Bem1_SGD_Costanzo_SL_Tong.xlsx") ## From yeastmine database
existing_data.columns=["query","organism","assay","role_query","interactor","type interaction","source","notes"]
existing_data_interactors_pos=existing_data[existing_data["type interaction"]=="Positive Genetic"].interactor.unique().tolist()
existing_data_interactors_neg=existing_data[existing_data["type interaction"]=="Negative Genetic"].interactor.unique().tolist()
existing_data_interactors_sl=existing_data[existing_data["type interaction"]=="Synthetic Lethality"].interactor.unique().tolist()

existing_data_interactors_neg_all=np.unique(existing_data_interactors_neg+existing_data_interactors_sl)
# -

bem1_PI_sig_SGA=pd.read_csv("../postprocessed-data/bem1_PI_sig_SGA.txt",sep="\t")
bem1_NI_sig_SGA=pd.read_csv("../postprocessed-data/bem1_NI_sig_SGA.txt",sep="\t")
# Take the interactors
bem1_PI_sig_SGA=bem1_PI_sig_SGA["Array allele name"].tolist()
bem1_NI_sig_SGA=bem1_NI_sig_SGA["Array allele name"].tolist()
# Convert to uppercase
bem1_PI_sig_SGA=[i.upper() for i in bem1_PI_sig_SGA]
bem1_NI_sig_SGA=[i.upper() for i in bem1_NI_sig_SGA]


# +
## Create a function to analyze the overlap between the existing data and the predicted data from SATAY 

def overlap_pos_GI(gi_SATAY_sig,gi_SATAY,existing_data_interactors_pos):

    overlap_pos_sig=[]
    overlap_pos=[]
    for i in gi_SATAY_sig:
        if i in existing_data_interactors_pos:
            overlap_pos_sig.append(i)
        elif i in gi_SATAY:
            overlap_pos.append(i)
            

    return overlap_pos_sig,overlap_pos

def overlap_neg_GI(gi_SATAY_sig,gi_SATAY,existing_data_interactors_neg):

    overlap_neg_sig=[]
    overlap_neg=[]
    for i in gi_SATAY_sig:
        if i in existing_data_interactors_neg:
            overlap_neg_sig.append(i)
        elif i in gi_SATAY:
            overlap_neg.append(i)
            

    return overlap_neg_sig,overlap_neg

def mismatch_pos2neg(gi_pos_SATAY_sig,gi_pos_SATAY,existing_data_interactors_neg):

    mismatch_pos2neg_sig=[]
    mismatch_pos2neg=[]
    for i in gi_pos_SATAY_sig:
        if i in existing_data_interactors_neg:
            mismatch_pos2neg_sig.append(i)
    for i in gi_pos_SATAY:
        if i in existing_data_interactors_neg:
            mismatch_pos2neg.append(i)
            

    return mismatch_pos2neg_sig,mismatch_pos2neg

def mismatch_neg2pos(gi_neg_SATAY_sig,gi_neg_SATAY,existing_data_interactors_pos):

    mismatch_neg2pos_sig=[]
    mismatch_neg2pos=[]
    for i in gi_neg_SATAY_sig:
        if i in existing_data_interactors_pos:
            mismatch_neg2pos_sig.append(i)
    for i in gi_neg_SATAY:
        if i in existing_data_interactors_pos:
            mismatch_neg2pos.append(i)
            

    return mismatch_neg2pos_sig,mismatch_neg2pos

# +
overlap_pos_sig,overlap_pos=overlap_pos_GI(bem1_PI_sig,bem1_PI,bem1_PI_sig_SGA)
overlap_neg_sig,overlap_neg=overlap_neg_GI(bem1_NI_sig,bem1_NI,bem1_NI_sig_SGA)

mismatch_pos2neg_sig,mismatch_pos2neg=mismatch_pos2neg(bem1_PI_sig,bem1_PI,bem1_NI_sig_SGA)
mismatch_neg2pos_sig,mismatch_neg2pos=mismatch_neg2pos(bem1_NI_sig,bem1_NI,bem1_PI_sig_SGA)

# -

overlap_pos_sig,overlap_pos,overlap_neg_sig,overlap_neg

mismatch_pos2neg_sig,mismatch_pos2neg,mismatch_neg2pos_sig,mismatch_neg2pos

len(overlap_pos)/len(bem1_PI_sig_SGA), len(overlap_neg)/len(bem1_NI_sig_SGA),len(mismatch_pos2neg)/len(bem1_PI_sig_SGA),len(mismatch_neg2pos)/len(bem1_NI_sig_SGA)

from annotate_volcano import annotate_volcano   #import annotate_volcano function
volcano_df=gi_pd_fitness_gene_bem3d
trackgene_list=["NRP1","IRC4","BEM1","WHI5"]
fig=annotate_volcano(volcano_df,figure_title="Interactors of dbem3 in WT",trackgene_list=trackgene_list)

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
genej_f_a=data_fitness2interact.loc[backg[0]].loc[genej,col_fitness]
genej_f_b=data_fitness2interact.loc[backg[1]].loc[genej,col_fitness]

## fitness bem3 in wt
genek_f_a=(data_fitness2interact.loc[backg[0]].loc[genek,col_fitness])
genek_f_b=(data_fitness2interact.loc[backg[1]].loc[genek,col_fitness])

### fitness normalization in the bem1-aid library to the median of the fitness of dbem1 in WT , and the dbem1dbem3 library 
# to the median of the fitness of dbem1dbem3 in WT

data_dgenej_a=data_fitness2interact.loc[backg[2]].loc[:,col_fitness]*genej_f_a/(np.nanmedian(data_backg.loc[backg[2]].loc[:,col_fitness]))
data_dgenej_b=data_fitness2interact.loc[backg[3]].loc[:,col_fitness]*genej_f_b/(np.nanmedian(data_backg.loc[backg[3]].loc[:,col_fitness]))

fitness_dgenek_dgenej_a=data_dgenej_a.loc[genek] # fitness of bem3 in the bem1-aid_a background
fitness_dgenek_dgenej_b=data_dgenej_b.loc[genek] # fitness of bem3 in the bem1-aid_b background

data_dgenek_a=data_fitness2interact.loc[backg[4]].loc[:,col_fitness]*genek_f_a/(np.nanmedian(data_backg.loc[backg[4]].loc[:,col_fitness]))
data_dgenek_b=data_fitness2interact.loc[backg[5]].loc[:,col_fitness]*genek_f_b/(np.nanmedian(data_backg.loc[backg[5]].loc[:,col_fitness]))

fitness_dgenej_dgenek_a=data_dgenek_a.loc[genej] # fitness of bem1 in the dbem3_a background
fitness_dgenej_dgenek_b=data_dgenek_b.loc[genej] # fitness of bem1 in the dbem3_b background

# data_dgenej_dgenek_a=data_fitness2interact.loc[backg[6]].loc[:,col_fitness]*0.5*(fitness_dgenek_dgenej_a+fitness_dgenej_dgenek_a)/(np.median(data_fitness.loc[backg[6]].loc[:,col_fitness]))
# data_dgenej_dgenek_b=data_fitness2interact.loc[backg[7]].loc[:,col_fitness]*0.5*(fitness_dgenek_dgenej_b+fitness_dgenej_dgenek_b)/(np.median(data_fitness.loc[backg[7]].loc[:,col_fitness]))
data_dgenej_dgenek_a=data_fitness2interact.loc[backg[6]].loc[:,col_fitness]*0.5*(fitness_dgenek_dgenej_a+fitness_dgenej_dgenek_a)/(np.nanmedian(data_backg.loc[backg[6]].loc[:,col_fitness]))
data_dgenej_dgenek_b=data_fitness2interact.loc[backg[7]].loc[:,col_fitness]*0.5*(fitness_dgenek_dgenej_b+fitness_dgenej_dgenek_b)/(np.nanmedian(data_backg.loc[backg[7]].loc[:,col_fitness]))

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
plt.hist(data_dgenej_dgenek_b,bins=70,alpha=0.6);
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
    
    ttest_val = stats.ttest_ind(variable_a_array, variable_b_array,equal_var=False) #T-test
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

gi_pd_fitness_gene_dbem1dbem3.dropna(inplace=True)
gi_pd_fitness_gene_dbem1dbem3.sort_values(by="fold_change",ascending=False)

# +
## sort by significance ==True

a=gi_pd_fitness_gene_dbem1dbem3[gi_pd_fitness_gene_dbem1dbem3.loc[:,"significance"]==True]
mean_fc=a.loc[:,"fold_change"].mean()
std_fc=a.loc[:,"fold_change"].std()

np.abs(mean_fc)+std_fc
# -

gi_pd_fitness_gene_dbem1dbem3[gi_pd_fitness_gene_dbem1dbem3.loc[:,"p_statistic"]<0.3].sort_values(by="fold_change",ascending=False)

# +
bem1bem3_PI_sig,bem1bem3_NI_sig,bem1bem3_PI,bem1bem3_NI,bem1bem3_PI_all,bem1bem3_NI_all=classify_GI(gi_pd_fitness_gene_dbem1dbem3,col="fold_change")

bem1bem3_PI2export=gi_pd_fitness_gene_dbem1dbem3.loc[bem1bem3_PI]
bem1bem3_PI2export=bem1bem3_PI2export[bem1bem3_PI2export["p_statistic"]<0.3]

bem1bem3_NI2export=gi_pd_fitness_gene_dbem1dbem3.loc[bem1bem3_NI]
bem1bem3_NI2export=bem1bem3_NI2export[bem1bem3_NI2export["p_statistic"]<0.3]

for i in gi_pd_fitness_gene_dbem1dbem3.index:
    if i in bem1bem3_NI2export.index:
        gi_pd_fitness_gene_dbem1dbem3.loc[i,"significance4FC"]="neg"
    elif i in bem1bem3_PI2export.index:
        gi_pd_fitness_gene_dbem1dbem3.loc[i,"significance4FC"]="pos"
    else:
        gi_pd_fitness_gene_dbem1dbem3.loc[i,"significance4FC"]="none"
# -

bem1bem3_NI2export.fold_change

# +
### Plot 2D Histograms of the fitness landscapes from wt-dbem1-dbem1dbem3

fig,ax=plt.subplots(1,1,figsize=(5,5))
cmap="seismic"

## select the common indexes among all dataframes

common_index=data_fitness_wt.index.intersection(data_fitness_dbem1.index).intersection(data_fitness_dbem1dbem3.index)

ax.hist2d(data_fitness_wt.loc[common_index],data_fitness_dbem1.loc[common_index],bins=100,cmap=cmap,cmin=1);

ax.set_xlabel("Fitness WT")
ax.set_ylabel("Fitness dbem1")
ax.grid(alpha=0.4)

fig.savefig("../figures/fitness_landscapes_2D_wt-dbem1.png",dpi=300,bbox_inches="tight")

fig,ax=plt.subplots(1,1,figsize=(5,5))
ax.hist2d(data_fitness_dbem1.loc[common_index],data_fitness_dbem1dbem3.loc[common_index],bins=100,cmap=cmap,cmin=1);
ax.set_xlabel("Fitness dbem1")
ax.set_ylabel("Fitness dbem1dbem3")
ax.grid(alpha=0.4)


#fig.savefig("../figures/fitness_landscapes_2D_wt-dbem1-dbem1dbem3.png",dpi=300,bbox_inches="tight")

# +
fig,ax=plt.subplots(1,1,figsize=(5,5))
cmap="seismic"

## select the common indexes among all dataframes

common_index=data_fitness_wt.index.intersection(data_fitness_dbem1.index).intersection(data_fitness_dbem1dbem3.index)

ax.scatter(data_fitness_wt.loc[common_index],data_fitness_dbem1.loc[common_index],color="gray",s=20,alpha=0.1);

ax.set_xlabel("Fitness WT (SATAY)")
ax.set_ylabel("Fitness dbem1 (SATAY)")
ax.grid(alpha=0.4)

for i in bem1_PI2export.index:
    if i in data_fitness_wt.index:
        ax.scatter(data_fitness_wt.loc[i],data_fitness_dbem1.loc[i],color="purple",s=20)
for i in bem1_NI2export.index:
    if i in data_fitness_wt.index:
        ax.scatter(data_fitness_wt.loc[i],data_fitness_dbem1.loc[i],color="green",s=20)

for i in bem1_PI_sig_SGA:
    if i in data_fitness_wt.index and i in data_fitness_dbem1.index:
        ax.scatter(data_fitness_wt.loc[i],data_fitness_dbem1.loc[i],color="darkturquoise",s=20,alpha=0.5)
for i in bem1_NI_sig_SGA:
    if i in data_fitness_wt.index and i in data_fitness_dbem1.index:
        ax.scatter(data_fitness_wt.loc[i],data_fitness_dbem1.loc[i],color="gold",s=20,alpha=0.5)

#fig.savefig("../figures/fitness_landscapes_2D_wt-dbem1_overlap_SGA_SATAY_interactors.png",dpi=300,bbox_inches="tight")
    
# -

len(common_index)*0.67

# +
### Plot 2D Histograms of the fitness landscapes from wt-dbem3 to appendix

fig,ax=plt.subplots(1,1,figsize=(5,5))


data_fitness_dbem3.dropna(inplace=True)
common_index_dbem3=data_fitness_wt.index.intersection(data_fitness_dbem3.index)

ax.hist2d(data_fitness_wt.loc[common_index_dbem3],
          data_fitness_dbem3.loc[common_index_dbem3],bins=100,cmap=cmap,cmin=1);
ax.set_xlabel("Fitness WT")
ax.set_ylabel("Fitness dbem3")
ax.grid(alpha=0.4)
plt.savefig("../figures/fitness_landscapes_2D_wt-dbem3.png",dpi=300,bbox_inches="tight")

# +
### Plot 2D Histograms of the fitness landscapes from dbem3-dbem1dbem3 to appendix

fig,ax=plt.subplots(1,1,figsize=(5,5))


data_fitness_dbem3.dropna(inplace=True)
common_index_dbem3=data_fitness_dbem1dbem3.index.intersection(data_fitness_dbem3.index)

ax.hist2d(data_fitness_dbem3.loc[common_index_dbem3],
          data_fitness_dbem1dbem3.loc[common_index_dbem3],
          bins=100,cmap=cmap,cmin=1);

ax.set_xlabel("Fitness dbem3")
ax.set_ylabel("Fitness dbem1dbem3")
ax.grid(alpha=0.4)

plt.savefig("../figures/fitness_landscapes_2D_dbem3-dbem1dbem3.png",dpi=300,bbox_inches="tight")

# +
fig,ax=plt.subplots(1,1,figsize=(5,5))
cmap="seismic"

## select the common indexes among all dataframes

common_index_interaction=gi_pd_fitness_gene_bem1d.index.intersection(gi_pd_fitness_gene_dbem1dbem3.index)

ax.hist2d(gi_pd_fitness_gene_bem1d.loc[common_index_interaction,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[common_index_interaction,"fold_change"],
          bins=100,cmap=cmap,cmin=1);

ax.set_xlabel("GI BEM1")
ax.set_ylabel("GI BEM1BEM3")
ax.grid(alpha=0.4)

#fig.savefig("../figures/interaction_landscapes_2D_dbem1-dbem1dbem3.png",dpi=300,bbox_inches="tight")

# +
fig,ax=plt.subplots(1,1,figsize=(5,5))
cmap="seismic"

## select the common indexes among all dataframes

common_index_interaction=gi_pd_fitness_gene_bem1d.index.intersection(gi_pd_fitness_gene_dbem1dbem3.index)

ax.hist2d(gi_pd_fitness_gene_bem1d.loc[common_index_interaction,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[common_index_interaction,"fold_change"],
          bins=100,cmap=cmap,cmin=1);

ax.set_xlabel("GI BEM1")
ax.set_ylabel("GI BEM1BEM3")
ax.grid(alpha=0.4)
ax.set_xlim(-5,10)
ax.set_ylim(-4,3)

size_scatter=10

for i in bem1_NI_sig:
    if i in common_index_interaction:
        ax.scatter(gi_pd_fitness_gene_bem1d.loc[i,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[i,"fold_change"],color="gray",s=size_scatter)
for i in bem1_PI_sig:
    if i in common_index_interaction:
        ax.scatter(gi_pd_fitness_gene_bem1d.loc[i,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[i,"fold_change"],color="orange",s=size_scatter)
for i in bem1bem3_NI2export.index:
    if i in common_index_interaction:
        ax.scatter(gi_pd_fitness_gene_bem1d.loc[i,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[i,"fold_change"],color="green",s=size_scatter)
for i in bem1bem3_PI2export.index:
    if i in common_index_interaction:
        ax.scatter(gi_pd_fitness_gene_bem1d.loc[i,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[i,"fold_change"],color="purple",s=size_scatter)

#fig.savefig("../figures/interaction_landscapes_2D_dbem1-dbem1dbem3_significant_genes_highlighted.png",dpi=300,bbox_inches="tight")

# +
## Making distributions of the product of the digenic and trigenic interactions, to know if the distribution of the product of the 
# digenic and trigenic significant interactions is significantly different from the overall distribution of the product of all  digenic and trigenic interactions

fig,ax=plt.subplots(1,2,figsize=(10,5))

common_index_interaction=gi_pd_fitness_gene_bem1d.index.intersection(gi_pd_fitness_gene_dbem1dbem3.index) 
product=gi_pd_fitness_gene_bem1d.loc[common_index_interaction,"fold_change"]*gi_pd_fitness_gene_dbem1dbem3.loc[common_index_interaction,"fold_change"]
ax[0].hist(product,bins=100,color="navy");



## distribution of significant interactions 

significant_interactions=[bem1_PI_sig,bem1_NI_sig,bem1bem3_PI2export.index,bem1bem3_NI2export.index]

## Flaten the list of lists
significant_interactions=[item for sublist in significant_interactions for item in sublist]

#ax[1].hist(gi_pd_fitness_gene_bem1d.loc[significant_interactions,"fold_change"]*gi_pd_fitness_gene_dbem1dbem3.loc[significant_interactions,"fold_change"],bins=10,alpha=0.5,color="orange");
signif_product=gi_pd_fitness_gene_bem1d.loc[significant_interactions,"fold_change"]*gi_pd_fitness_gene_dbem1dbem3.loc[significant_interactions,"fold_change"]
sign_epistasis_signif=signif_product[signif_product<0]
no_sign_epistasis_signif=signif_product[signif_product>0]

bins=np.array([-np.inf,0,np.inf])
colors=["#D0D3D67A", "#FEF0DE7A"]

ax[1].hist(signif_product,bins=10,color="#FEF0DE7A");




## Create a t-test to compare both distributions to see if the one from the significant interactions is significantly different from the overall distribution

t_stat,p_value = stats.ttest_ind(product,signif_product,equal_var=False) #T-test

print(f"T-statistic: {t_stat}, P-value: {p_value}")
fig.savefig("../figures/ttest-sign-epistasis.png",dpi=300,bbox_inches="tight")

# +
fig,ax=plt.subplots(1,1,figsize=(5,5))
cmap="seismic"

## select the common indexes among all dataframes

common_index_interaction=gi_pd_fitness_gene_bem3d.index.intersection(gi_pd_fitness_gene_dbem1dbem3.index)

ax.hist2d(gi_pd_fitness_gene_bem3d.loc[common_index_interaction,"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[common_index_interaction,"fold_change"],
          bins=100,cmap=cmap,cmin=1);

ax.set_xlabel("GI BEM3")
ax.set_ylabel("GI BEM1BEM3")
ax.grid(alpha=0.4)

fig.savefig("../figures/interaction_landscapes_2D_dbem3-dbem1dbem3.png",dpi=300,bbox_inches="tight")

# +
## count how many genes invert their sign from positive to negative and viceversa from bem1 to bem1bem3

e_dbem1=gi_pd_fitness_gene_bem1d.loc[common_index_interaction,"fold_change"]
e_dbem1dbem3=gi_pd_fitness_gene_dbem1dbem3.loc[common_index_interaction,"fold_change"]

sign_epistatsis_genes_pos_dbem1_neg_dbem13=e_dbem1[(e_dbem1>0)].index.intersection(e_dbem1dbem3[(e_dbem1dbem3<0)].index)

sign_epistatsis_genes_neg_dbem1_pos_dbem13=e_dbem1[(e_dbem1<0)].index.intersection(e_dbem1dbem3[(e_dbem1dbem3>0)].index)

conserve_sign_pos=e_dbem1[(e_dbem1>0)].index.intersection(e_dbem1dbem3[(e_dbem1dbem3>0)].index)

conserve_sign_neg=e_dbem1[(e_dbem1<0)].index.intersection(e_dbem1dbem3[(e_dbem1dbem3<0)].index)

len(sign_epistatsis_genes_pos_dbem1_neg_dbem13)/len(common_index),len(sign_epistatsis_genes_neg_dbem1_pos_dbem13)/len(common_index),len(conserve_sign_pos)/len(common_index),len(conserve_sign_neg)/len(common_index)

# +
## count how many genes invert their sign from positive to negative and viceversa from bem1 to bem1bem3 (ONLY SIGNIFICANT GENES)

bem1_significant=bem1_PI_sig.union(bem1_NI_sig)
bem1bem3_significant=bem1bem3_PI2export.index.union(bem1bem3_NI2export.index)
all_significant_genes=bem1_significant.union(bem1bem3_significant)


e_dbem1=gi_pd_fitness_gene_bem1d.loc[all_significant_genes,"fold_change"]
e_dbem1dbem3=gi_pd_fitness_gene_dbem1dbem3.loc[all_significant_genes,"fold_change"]

sign_epistatsis_genes_pos_dbem1_neg_dbem13=e_dbem1[(e_dbem1>0)].index.intersection(e_dbem1dbem3[(e_dbem1dbem3<0)].index)

sign_epistatsis_genes_neg_dbem1_pos_dbem13=e_dbem1[(e_dbem1<0)].index.intersection(e_dbem1dbem3[(e_dbem1dbem3>0)].index)

conserve_sign_pos=e_dbem1[(e_dbem1>0)].index.intersection(e_dbem1dbem3[(e_dbem1dbem3>0)].index)

conserve_sign_neg=e_dbem1[(e_dbem1<0)].index.intersection(e_dbem1dbem3[(e_dbem1dbem3<0)].index)

len(sign_epistatsis_genes_pos_dbem1_neg_dbem13)/len(common_index),len(sign_epistatsis_genes_neg_dbem1_pos_dbem13)/len(common_index),len(conserve_sign_pos)/len(common_index),len(conserve_sign_neg)/len(common_index)
# -

len(all_significant_genes), len(sign_epistatsis_genes_pos_dbem1_neg_dbem13)

# +
## make a pie plot of the sign of the epistasis from bem1 to bem1bem3

labels=["Pos2Pos","Neg2Pos","Neg2Neg","Pos2Neg"]
sizes=[len(conserve_sign_pos)/len(all_significant_genes),
       len(sign_epistatsis_genes_neg_dbem1_pos_dbem13)/len(all_significant_genes),
       len(conserve_sign_neg)/len(all_significant_genes),
       len(sign_epistatsis_genes_pos_dbem1_neg_dbem13)/len(all_significant_genes),]

fig,ax=plt.subplots(1,1,figsize=(5,5))
colors=["#D0D3D67A", "#FEF0DE7A","#D0D3D67A","#FEF0DE7A"]

#ax.pie(sizes,labels=labels,autopct='%1.1f%%', startangle=90,colors=colors);
g=ax.pie(sizes,labels=labels, startangle=90,colors=colors);

fig.savefig("../figures/interaction_signs_pie_dbem1-dbem1dbem3_significant_genes.png",dpi=300,bbox_inches="tight")
# -

sign_epistatsis_genes_neg_dbem1_pos_dbem13

# +
common_index_sig_PI=np.intersect1d(bem1_PI2export.index,bem1bem3_PI2export.index)
common_index_sig_NI=np.intersect1d(bem1_NI2export.index,bem1bem3_NI2export.index)
common_index_sig_bem1_PI_bem1bem3_NI=np.intersect1d(bem1_PI2export.index,bem1bem3_NI2export.index)
common_index_sig_bem1_NI_bem1bem3_PI=np.intersect1d(bem1_NI2export.index,bem1bem3_PI2export.index)

common_index_sig_NI,common_index_sig_bem1_PI_bem1bem3_NI,common_index_sig_bem1_NI_bem1bem3_PI
# -

with open("../postprocessed-data/positive_satay_genes_bem1bem3.txt","w") as f:
    for i in bem1bem3_PI2export.index:
        f.write(i+"\n")
f.close()
with open("../postprocessed-data/negative_satay_genes_bem1bem3.txt","w") as f:
    for i in bem1bem3_NI2export.index:
        f.write(i+"\n")
f.close()

bem1_PI2export.index,bem1bem3_NI2export.index

# +
## Make a hierchical clustering of the most correlated genes with bem1 , bem3 and dbem1dbem3

from scipy.cluster import hierarchy
from scipy.spatial import distance

intersection_genes=list(set(gi_pd_fitness_gene_dbem1dbem3.index)& set(gi_pd_fitness_gene_bem1d)& set(gi_pd_fitness_gene_bem3d.index))

## make a dataframe with the intersection genes and the interaction scores of bem1 and bem3 and in dbem1dbem3
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import linkage, dendrogram

# Create a DataFrame from the gene data
#genes_data=[bem1_PI,bem1_NI,bem3_PI,bem3_NI,bem1bem3_PI,bem1bem3_NI]
genes_data=[gi_pd_fitness_gene_bem1d.index,gi_pd_fitness_gene_bem3d.index,gi_pd_fitness_gene_dbem1dbem3.index]
genes_data= np.unique(np.concatenate(genes_data))


## create a dtaframe with those genes and the interaction scores of bem1 and bem3 and in dbem1dbem3

similarity_profile_pd=pd.DataFrame(index=genes_data,columns=["fold_change_dbem1","fold_change_dbem3","fold_change_dbem1dbem3"])



for i in genes_data:
    if i in gi_pd_fitness_gene_bem1d.index:
        similarity_profile_pd.loc[i,"fold_change_dbem1"]=gi_pd_fitness_gene_bem1d.loc[i,"fold_change"]
    if i in gi_pd_fitness_gene_bem3d.index:
        similarity_profile_pd.loc[i,"fold_change_dbem3"]=gi_pd_fitness_gene_bem3d.loc[i,"fold_change"]
    if i in gi_pd_fitness_gene_dbem1dbem3.index:
        similarity_profile_pd.loc[i,"fold_change_dbem1dbem3"]=gi_pd_fitness_gene_dbem1dbem3.loc[i,"fold_change"]

similarity_profile_pd.dropna(inplace=True)

# Calculate the correlation matrix
correlation_matrix = similarity_profile_pd.T.corr(method='pearson')


# +
from scipy.spatial import distance
from scipy.cluster import hierarchy

cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)

correlation_matrix[correlation_matrix<0]=0
correlations_array = correlation_matrix

row_linkage = hierarchy.linkage(
    distance.pdist(correlations_array), method='average')

col_linkage = hierarchy.linkage(
    distance.pdist(correlations_array.T), method='average')

g=sns.clustermap(correlations_array, row_linkage=row_linkage, col_linkage=col_linkage,figsize=(13, 13), cmap=cmap,vmin=0)
g.savefig("../figures/correlation_matrix_bem1_bem3_dbem1dbem3.png",dpi=300,bbox_inches="tight")

# +
## To extract the clusters from the correlation matrix

# from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

# # Perform hierarchical clustering
# linkage_matrix = linkage(g.data2d.corr(), method='average')

# # Cut the dendrogram to obtain a specific number of clusters
# num_clusters = 20  # Change this to the desired number of clusters
# row_clusters = fcluster(linkage_matrix, num_clusters, criterion='maxclust')

# # Create a DataFrame with the cluster assignments

# df = pd.DataFrame({'Row': g.data2d.index, 'Cluster': row_clusters})

# x=df[df["Cluster"]==4].Row.tolist()

# sns.heatmap(correlation_matrix.loc[x,x],cmap=cmap)

# +
mutations_laan=pd.read_excel("../postprocessed-data/elife-09638-supp1-v2.xlsx")
mutations_laan.gene.unique()

mutations_laan2save=["ERV29","BEN3","PKH3","NRP1",
                     "BEM2","LRG1","SOL4","YOX1","PLC1",
                     "YKL091C","GPR1","MPH3","BEM2","TCB1",
                     "IRA1","GAL83"]
# -

with open("../postprocessed-data/mutations-laan-evolutionary-trajectory.txt","w") as f:
    for i in mutations_laan2save:
        f.write(i+"\n")
f.close()

mutations_laan[mutations_laan.loc[:,"gene"]=="GPR1"]

# +
## Plot a barplot of the interaction scores with BEM1 and BEM1BEM3 of specific genes related to the evolutionary trajectory

genes=['ERV29', 'PKH3', 'NRP1', 'BEM2','LRG1','SOL4',
       'YOX1', 'PLC1', 'YKL091C', 'GPR1', 'MPH3', 'TCB1',
       'IRA1', 'GAL83']

scores_bem1=[]
p_bem1=[]
scores_bem1bem3=[]
p_bem1bem3=[]
categories=[]
for gene in genes:
    if gene in gi_pd_fitness_gene_bem1d.index and gene in gi_pd_fitness_gene_dbem1dbem3.index:
        scores_bem1.append(gi_pd_fitness_gene_bem1d.loc[gene,"fold_change"])
        p_bem1.append(gi_pd_fitness_gene_bem1d.loc[gene,"p_statistic"])
        scores_bem1bem3.append(gi_pd_fitness_gene_dbem1dbem3.loc[gene,"fold_change"])
        p_bem1bem3.append(gi_pd_fitness_gene_dbem1dbem3.loc[gene,"p_statistic"])
        categories.append(gene)
    

plt.figure(figsize=(10,5))

scores_bem1.append(gi_pd_fitness_gene_bem1d.loc["BEM3","fold_change"])
scores_bem1bem3.append(0)

categories=categories+["BEM3"]
variable1 = scores_bem1

variable2 = scores_bem1bem3

bar_width = 0.35  # width of the bars

fig, ax = plt.subplots()

# Plotting the bars
bar1 = np.arange(len(variable1))
bar2 = [x + bar_width for x in bar1]
ax.bar(bar1, variable1, bar_width, label='BEM1',color="magenta")
ax.bar(bar2, variable2, bar_width, label='BEM1BEM3',color="gold")

ax.set_xticks([r + bar_width/2 for r in range(len(variable1))])
ax.set_xticklabels(categories)
ax.legend()


plt.ylabel('Genetic interaction scores')
# plt.title('Bar Plot with Two Variables')
plt.tight_layout()

plt.savefig("../figures/non-fixed-mutations-trajectory-satay.png",dpi=300,bbox_inches="tight")
# -

gi_pd_fitness_gene_bem1d.loc[categories,"p_statistic"],gi_pd_fitness_gene_dbem1dbem3.loc[categories[0:-1],"p_statistic"]

# +
all_pathways_gi=gi_pd_fitness_gene_dbem1dbem3[gi_pd_fitness_gene_dbem1dbem3.loc[:,"gene_names"].isin(all_pathways)]


all_pathways_gi.to_csv("../postprocessed-data/pathways_bem1bem3_gi.csv",sep="\t")

# +
volcano_df=gi_pd_fitness_gene_dbem1dbem3

fig=annotate_volcano(volcano_df,figure_title="Interactors of dbem1dbem3 in WT")

fig.savefig("../figures/volcano_dbem1dbem3.png",dpi=300,bbox_inches="tight")
# -

gi_pd_fitness_gene_dbem1dbem3[gi_pd_fitness_gene_dbem1dbem3.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)

volcano_df=gi_pd_fitness_gene_dbem1dbem3
trackgene_list=["NRP1","BUD5","ERG8","FLO11"]
fig=annotate_volcano(volcano_df,figure_title="Interactors of dbem1dbem3 in WT",trackgene_list=trackgene_list)

# +
common_index=[]
for i in gi_pd_fitness_gene_dbem1dbem3.index:
    if i in gi_pd_fitness_gene_bem1d.index:
        common_index.append(i)

gi_dbem1=gi_pd_fitness_gene_bem1d.loc[common_index,"fold_change"]
gi_dbem1dbem3=gi_pd_fitness_gene_dbem1dbem3.loc[common_index,"fold_change"]


gi_supression_shift=set(bem1_NI).intersection(bem1bem3_PI)

gi_essentiality_shift=set(bem1_PI).intersection(bem1bem3_NI)

len(gi_essentiality_shift),len(gi_supression_shift)

# +
## Plot differences in genetic interaction scores
from matplotlib_venn import venn3, venn3_circles

a1=len(set(bem1_PI))
a2=len(set(bem1bem3_PI))
a3=len(set(bem1_PI).intersection(set(bem1bem3_PI)))
a4=len(set(bem1_NI))
a5=len(set(bem1_PI).intersection(set(bem1_NI)))
a6=len(set(bem1bem3_PI).intersection(set(bem1_NI)))
a7=len(set(bem1_PI).intersection(set(bem1bem3_PI)).intersection(set(bem1_NI)))

a=[a1,a2,a3,a4,a5,a6,a7]
plt.figure(figsize=(10,8))
g=venn3(subsets = a, set_labels = ('BEM1 PI', 'BEM1BEM3 PI','BEM1 NI'), alpha = 0.5,
set_colors=('purple', 'yellow', 'green'));

for x in range(len(g.subset_labels)):
    if g.subset_labels[x] is not None:
        g.subset_labels[x].set_fontsize(15)
#plt.savefig("../figures/fig_venn_bem1_bem1bem3_PI_relationships.png",dpi=300)

# +
from matplotlib_venn import venn3, venn3_circles

a1=len(set(bem1_PI))
a2=len(set(bem1bem3_NI))
a3=len(set(bem1_PI).intersection(set(bem1bem3_NI)))
a4=len(set(bem1_NI))
a5=len(set(bem1_PI).intersection(set(bem1_NI)))
a6=len(set(bem1bem3_NI).intersection(set(bem1_NI)))
a7=len(set(bem1_PI).intersection(set(bem1bem3_NI)).intersection(set(bem1_NI)))

a=[a1,a2,a3,a4,a5,a6,a7]
plt.figure(figsize=(10,8))
g=venn3(subsets = a, set_labels = ('BEM1 PI', 'BEM1BEM3 NI','BEM1 NI'), alpha = 0.5,
set_colors=('purple', 'yellow', 'green'));
for x in range(len(g.subset_labels)):
    if g.subset_labels[x] is not None:
        g.subset_labels[x].set_fontsize(15)

#plt.savefig("../figures/fig_venn_bem1_bem1bem3_NI_relationships.png",dpi=300)

# +
## Plot differences in genetic interaction scores
from matplotlib_venn import venn3, venn3_circles

a1=len(set(bem3_PI))
a2=len(set(bem1bem3_PI))
a3=len(set(bem3_PI).intersection(set(bem1bem3_PI)))
a4=len(set(bem3_NI))
a5=len(set(bem3_PI).intersection(set(bem3_NI)))
a6=len(set(bem1bem3_PI).intersection(set(bem3_NI)))
a7=len(set(bem3_PI).intersection(set(bem1bem3_PI)).intersection(set(bem3_NI)))

a=[a1,a2,a3,a4,a5,a6,a7]
plt.figure(figsize=(10,8))
g=venn3(subsets = a, set_labels = ('BEM3 PI', 'BEM1BEM3 PI','BEM3 NI'), alpha = 0.5,
set_colors=('purple', 'yellow', 'green'));

for x in range(len(g.subset_labels)):
    if g.subset_labels[x] is not None:
        g.subset_labels[x].set_fontsize(15)

#plt.savefig("../figures/fig_venn_bem3_bem1bem3_PI_relationships.png",dpi=300)

# +
## Plot differences in genetic interaction scores
from matplotlib_venn import venn3, venn3_circles

a1=len(set(bem3_PI))
a2=len(set(bem1bem3_NI))
a3=len(set(bem3_PI).intersection(set(bem1bem3_NI)))
a4=len(set(bem3_NI))
a5=len(set(bem3_PI).intersection(set(bem3_NI)))
a6=len(set(bem1bem3_NI).intersection(set(bem3_NI)))
a7=len(set(bem3_PI).intersection(set(bem1bem3_NI)).intersection(set(bem3_NI)))

a=[a1,a2,a3,a4,a5,a6,a7]
plt.figure(figsize=(10,8))
g=venn3(subsets = a, set_labels = ('BEM3 PI', 'BEM1BEM3 NI','BEM3 NI'), alpha = 0.5,
set_colors=('purple', 'yellow', 'green'));

for x in range(len(g.subset_labels)):
    if g.subset_labels[x] is not None:
        g.subset_labels[x].set_fontsize(15)

#plt.savefig("../figures/fig_venn_bem3_bem1bem3_NI_relationships.png",dpi=300)

# +
## how many interactions are lost from dbem1 to dbem1dbem3 and from dbem3 to dbem1dbem3

bem1bem3_all=bem1bem3_NI_sig.union(bem1bem3_PI_sig)
bem1_all=bem1_PI_all.union(bem1_NI_all)
bem3_all=bem3_PI_all.union(bem3_NI_all)
# Gained interactors for bem1bem3
bem1bem3_gained_dbem1=bem1bem3_all.difference(bem1_all)
bem1bem3_gained_dbem3=bem1bem3_all.difference(bem3_all)
# Lost interactors for bem1bem3
bem1bem3_lost_dbem1=bem1_all.difference(bem1bem3_all)
bem1bem3_lost_dbem3=bem3_all.difference(bem1bem3_all)
# -

len(bem1bem3_lost_dbem3),len(bem1bem3_gained_dbem3),len(bem1bem3_lost_dbem1),len(bem1bem3_gained_dbem1)

list(gi_supression_shift)

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
plt.scatter(gi_pd_fitness_gene_bem1d.loc[list(gi_supression_shift),"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[list(gi_supression_shift),"fold_change"],
s=50,alpha=1,c="gray",edgecolors="black",linewidth=2)
plt.scatter(gi_pd_fitness_gene_bem1d.loc[list(gi_essentiality_shift),"fold_change"],gi_pd_fitness_gene_dbem1dbem3.loc[list(gi_essentiality_shift),"fold_change"],
s=50,alpha=1,c="gray",edgecolors="black",linewidth=2)
plt.xlabel("dbem1 interaction scores",fontsize=14)
plt.ylabel("dbem1dbem3 interaction scores",fontsize=14)

plt.vlines(0,-8,4,color="black",linewidth=2,linestyles="dashed")
plt.hlines(0,-8,25,color="black",linewidth=2,linestyles="dashed")

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

# +
## wich genes flip their interaction from dbem1 to data_dbem1dbem3

gene_pos_1_neg_13_sig=list(set(bem1_PI_sig).intersection(set(bem1bem3_NI_sig)))

gene_neg_1_pos_13_sig=list(set(bem1_NI_sig).intersection(set(bem1bem3_PI_sig)))

gene_pos_1_neg_13=list(set(bem1_PI).intersection(set(bem1bem3_NI)))

gene_neg_1_pos_13=list(set(bem1_NI).intersection(set(bem1bem3_PI)))

len(gene_pos_1_neg_13),len(gene_neg_1_pos_13),len(gene_pos_1_neg_13_sig),len(gene_neg_1_pos_13_sig)
# -

gi_pd_fitness_gene_bem1d.loc[gene_pos_1_neg_13].sort_values(by="fold_change",ascending=False)

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
## Changes across polarity proteins in different backgrounds 

p_b1=polarity_genes_gi_pd_dbem1.fold_change

p_b3=polarity_genes_gi_pd_dbem3.fold_change

p_b1b3=polarity_genes_gi_pd_dbem1dbem3.fold_change

p_value_b1=polarity_genes_gi_pd_dbem1.p_statistic

p_value_b3=polarity_genes_gi_pd_dbem3.p_statistic

p_value_b1b3=polarity_genes_gi_pd_dbem1dbem3.p_statistic

### Build a dataframe with the fold changes of the polarity genes in the different backgrounds

p=pd.concat([p_b1,p_b3,p_b1b3],axis=1)

p_value_polarity=pd.concat([p_value_b1,p_value_b3,p_value_b1b3],axis=1)

p.columns=["dbem1","dbem3","dbem1dbem3"]

p_value_polarity.columns=["dbem1","dbem3","dbem1dbem3"]

p.dropna(inplace=True)

p_value_polarity.dropna(inplace=True)
# -

p

# +
cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
fig,ax=plt.subplots(1,1,figsize=(10,10))

p_new=p.drop(columns=["dbem3"])
p_new=p_new.sort_values(by="dbem1dbem3",ascending=False)
plt.imshow(p_new,cmap=cmap,vmin=-2,vmax=2)
plt.colorbar()
ax.set_yticks(np.arange(len(p.index)));
ax.set_xticks(np.arange(len(p_new.columns)));
ax.set_xticklabels(p_new.columns,fontsize=14,rotation=90)
ax.set_yticklabels(p.index,fontsize=14);


fig.savefig("../figures/polarity_genes_heatmap_dbem1_dbem3_dbem1dbem3.png",dpi=300,bbox_inches="tight")

# +
## Plot a 3d scatter plot of the interaction scores of the polarity genes in the different backgrounds

from mpl_toolkits.mplot3d import Axes3D

fig,ax = plt.subplots(1,1,figsize=(5,5))

X=p.iloc[:,0]
Y=p.iloc[:,1]
g=ax.scatter(X,Y,c=p.iloc[:,2],cmap=cmap,vmin=-2.5,vmax=2.5,s=50)
ax.set_xlabel("dbem1")
ax.set_ylabel("dbem3")
ax.set_title("Interaction score polarity genes")

for i in p_value_polarity.index:
    if p_value_polarity.loc[i,"dbem1dbem3"]<0.2:
        ax.scatter(p.loc[i,"dbem1"],p.loc[i,"dbem3"],c="black",marker="*",s=50)
        ax.annotate(i,(p.loc[i,"dbem1"],p.loc[i,"dbem3"]))

cbar = plt.colorbar(g)
cbar.ax.set_ylabel('dbem1dbem3', rotation=270)

# +
from mpl_toolkits.mplot3d import Axes3D

fig,ax = plt.subplots(1,1,figsize=(5,5))

X=p.iloc[:,0]
Y=p.iloc[:,2]
g=ax.scatter(X,Y,s=50,c="black",alpha=0.5)
ax.set_xlabel("dbem1")
ax.set_ylabel("dbem1dbem3")
ax.set_title("Interaction score polarity genes")

for i in p_value_polarity.index:
    if p_value_polarity.loc[i,"dbem1dbem3"]<0.3 and p_value_polarity.loc[i,"dbem1"]<0.3:
        ax.scatter(p.loc[i,"dbem1"],p.loc[i,"dbem1dbem3"],c="black",marker="*",s=50)
        ax.annotate(i,(p.loc[i,"dbem1"],p.loc[i,"dbem1dbem3"]))

fig.savefig("../figures/polarity_genes_scatter_dbem1_dbem1dbem3.png",dpi=300,bbox_inches="tight")
# -

polarity_genes



# +
import seaborn as sns
x=p

cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)

g=sns.clustermap(x,xticklabels=x.columns,cbar=True,annot=False,cmap=cmap,vmax=x.max().max(),vmin=x.min().min(),figsize=(8,15))

g.savefig("../figures/fig_clustermap_polarity_genes_all_backgrounds.png",dpi=300)

# +
## Esentiality prediction across backgrounds

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

fitness_bem1_wt=data_fitness.loc["wt_merged"].loc["BEM1","fitness_gene"]

fitness_bem3_wt=data_wt.loc["BEM3"]

## Normalize the datasetss such as the median values corresponds to the fitness of the knockout of the gene of interest 
data_wt_norm=data_wt

data_dbem1_norm=data_dbem1*fitness_bem1_wt/np.nanmedian(data_dbem1)
data_dbem3_norm=data_dbem3*fitness_bem3_wt/np.nanmedian(data_dbem3)

fitness_bem3_dbem1=data_dbem1_norm.loc["BEM3"]

data_dbem1dbem3_norm=data_dbem1dbem3*fitness_bem3_dbem1/np.nanmedian(data_dbem1dbem3.tolist())


# +

data_fitness_wta=data_fitness.loc["wt_a","fitness_gene"]
data_fitness_wtb=data_fitness.loc["wt_b","fitness_gene"]


data_fitness_dbem1a=data_fitness.loc["bem1-aid_a","fitness_gene"]
data_fitness_dbem1b=data_fitness.loc["bem1-aid_b","fitness_gene"]


data_fitness_dbem1dbem3a=data_fitness.loc["dbem1dbem3_a","fitness_gene"]
data_fitness_dbem1dbem3b=data_fitness.loc["dbem1dbem3_b","fitness_gene"]

data_fitness_dbem3a=data_fitness.loc["dbem3_a","fitness_gene"]
data_fitness_dbem3b=data_fitness.loc["dbem3_b","fitness_gene"]

# +
data_fitness_wta_essentials=data_fitness_wta[data_fitness_wta.index.isin(standard_essentials)]

data_fitness_wtb_essentials=data_fitness_wtb[data_fitness_wtb.index.isin(standard_essentials)]

data_fitness_wt_essentials=pd.concat([data_fitness_wta_essentials,data_fitness_wtb_essentials],axis=1)
data_fitness_wt_essentials.fillna(0,inplace=True)
# -

data_fitness_wt_essentials["mean"]=data_fitness_wt_essentials.mean(axis=1)
cutoff=np.round(np.mean(data_fitness_wt_essentials["mean"]),3)
# cutoff=0.44 # from the mean - std of the merged wt dataset
cutoff

# +
plt.figure(figsize=(5,5))


plt.hist(data_fitness_wt_essentials["mean"],bins=40,color="black",alpha=0.6);
plt.vlines(cutoff,0,100,color="gray",linestyles="dashed")#
plt.annotate(cutoff,xy=(cutoff,100),xytext=(cutoff+0.01,100),color="gray")
plt.xlabel("Fitness")
plt.ylabel("Number of essential genes")
plt.tight_layout()


# +
essential_genes_wta=data_fitness_wta[data_fitness_wta<cutoff].index
essential_genes_wtb=data_fitness_wtb[data_fitness_wtb<cutoff].index

essential_genes_dbem1a=data_fitness_dbem1a[data_fitness_dbem1a<cutoff].index
essential_genes_dbem1b=data_fitness_dbem1b[data_fitness_dbem1b<cutoff].index

essential_genes_dbem1dbem3a=data_fitness_dbem1dbem3a[data_fitness_dbem1dbem3a<cutoff].index
essential_genes_dbem1dbem3b=data_fitness_dbem1dbem3b[data_fitness_dbem1dbem3b<cutoff].index

essential_genes_dbem3a=data_fitness_dbem3a[data_fitness_dbem3a<cutoff].index
essential_genes_dbem3b=data_fitness_dbem3b[data_fitness_dbem3b<cutoff].index


essential_genes_wt=set(essential_genes_wta)&set(essential_genes_wtb)
essential_genes_dbem1=set(essential_genes_dbem1a)&set(essential_genes_dbem1b)
essential_genes_dbem1dbem3=set(essential_genes_dbem1dbem3a)&set(essential_genes_dbem1dbem3b)
essential_genes_dbem3=set(essential_genes_dbem3a)&set(essential_genes_dbem3b)


# +
number_e_genes=[len(essential_genes_wt),len(essential_genes_dbem1),len(essential_genes_dbem1dbem3),
len(essential_genes_dbem3)]

std_number_e_genes=[np.std([len(essential_genes_wta),len(essential_genes_wtb)]),
np.std([len(essential_genes_dbem1a),len(essential_genes_dbem1b)]),
np.std([len(essential_genes_dbem1dbem3a),len(essential_genes_dbem1dbem3b)]),
np.std([len(essential_genes_dbem3a),len(essential_genes_dbem3b)])]

# +
plt.figure(figsize=(8,5))

plt.bar(np.arange(len(number_e_genes)),number_e_genes,yerr=std_number_e_genes,color="black",alpha=0.6)

plt.xticks(np.arange(len(number_e_genes)),["WT","dbem1","dbem1dbem3","dbem3"]);

# +
plt.figure(figsize=(8,5))

number_e_genes_refactor=[number_e_genes[1],number_e_genes[2],number_e_genes[0],number_e_genes[1]+number_e_genes[3]]
std_number_e_genes_refactor=[std_number_e_genes[1],std_number_e_genes[2],std_number_e_genes[0],std_number_e_genes[1]+std_number_e_genes[3]]
plt.errorbar(["bem1$\Delta$","bem1$\Delta$bem3$\Delta$","WT","Expected \n in dbem1dbem3"],number_e_genes_refactor,
yerr=std_number_e_genes_refactor,fmt="o",color="black",capsize=10)
plt.ylabel("Predicted number of essential genes")

plt.plot(number_e_genes_refactor,"--",color="black",linewidth=0.5)
plt.yscale("log")

plt.tight_layout()

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

plt.ylim(0,1000)
# plt.yticks([1000,2000,3000],fontsize=14)
plt.xticks([0.0022,0.010,0.01244],fontsize=14)
plt.grid(linewidth=0.3,linestyle="dashed")
plt.xlabel("Growth rate (min$^{-1}$)")
plt.ylabel("Predicted number of essential genes")


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

# ## This notebook will implement the major influencing features for essentiality prediction using ML , which are :
#
# - Neighboorhood index : Number of transposon insertions within the ORF, normalized by the length of the ORF and the surrounding 10kbp. 
# - Freedom index :  Length of the largest insertion-free region in the ORF (Open Reading Frame), divided by the ORF’s length.
#
# From paper : Levitan A, et al. (2020) Comparing the utility of in vivo transposon mutagenesis approaches in yeast species to infer gene essentiality.

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
from functions_scores_essentiality import write_ones_if_essential

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

standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str)

#######  Import fitness from SATAY #########
import pickle
with open("../postprocessed-data/fitness_models_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

satay_fitness=pd.concat(b,axis=0,keys=keys)
# -

polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")
polarity_genes.fillna(0,inplace=True)

data_norm_pd=pd.read_excel("../postprocessed-data/data_norm_linear_transformation_per_background.xlsx",
engine='openpyxl',index_col="background")
data_norm_pd.drop(columns=["Unnamed: 0","Unnamed: 1"],inplace=True)


# +
list_data_pd_wt=list_data_pd.loc["wt_merged"]

ni_essentials=[] 
ni=data_norm_pd.loc["wt_merged"]["tr_normalized_windows"].values # neighborhood index

ho=np.where(list_data_pd_wt.loc[:,"Gene name"]=="HO")[0][0]

ni_ho=data_norm_pd.loc["wt_merged"]["tr_normalized_windows"].values[ho]
ni_norm=ni/ni_ho


genes=list_data_pd_wt.loc[:,"Gene name"]
genes.reset_index(drop=True,inplace=True)

for i in standard_essentials:
    x=np.where(genes==i)[0]
    
    ni_essentials.append(ni_norm[x])


ni_essentials=np.array(ni_essentials)

ni_essentials_bellow_ho=ni_essentials[(ni_essentials<np.array(0.5))]
ni_bellow_ho=ni_norm[(ni_norm<np.array(0.5))]

d=len(ni_essentials_bellow_ho)/len(standard_essentials)
d_ni=len(ni_bellow_ho)/len(ni)

# -

len(ni_bellow_ho)
d_ni

# +
fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

axes.hist(ni/ni_ho,bins=2000,alpha=0.5,color="gray",label="essential genes");
axes.vlines(0.5,0,300,color="red",linestyle="dashed",linewidth=2,label="HO/2");
axes.annotate(f"{d_ni*100:.2f}" +"%",xy=(0,270),fontsize=16)
axes.set_title("Distribution of NI all genes",fontsize=16)
axes.tick_params(axis="both",labelsize=16)
axes.set_xlabel("Relative Neighborhood Index",fontsize=16)
axes.set_xlim(0,5)
axes.set_ylabel("Count",fontsize=16)
axes.legend(loc="best",fontsize=16)


# These are in unitless percentages of the figure size. (0,0 is bottom left)
left, bottom, width, height = [0.5, 0.3, 0.4, 0.4]
ax2 = fig.add_axes([left, bottom, width, height])

data = [len(ni_essentials_bellow_ho),len(standard_essentials)-len(ni_essentials_bellow_ho)]
labels = ['Truly \n essentials'," "]

colors=["pink","gray"]



ax2.pie(data, labels = labels, colors = colors, autopct='%.0f%%',textprops={'fontsize': 14});
plt.tight_layout()

#fig.savefig("../figures/figures_thesis_chapter_2/fig_distribution_NI_for_all_pi_inset_essentials",dpi=300)

# +
## Distributio of NI for annotated essential genes

# fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

# axes.hist(np.concatenate(ni_essentials),bins=200,alpha=0.5,color="gray",label="essential genes");
# axes.vlines(0.5,0,200,color="red",linestyle="dashed",linewidth=2,label="HO/2");
# axes.annotate(f"{d*100:.2f}" +"%",xy=(0.1,175),fontsize=16)
# axes.set_title("Distribution of normalized insertions for annotated essential genes",fontsize=16)
# axes.tick_params(axis="both",labelsize=16)
# axes.set_xlabel("Normalized insertions by a 10kB windows compared to HO locus",fontsize=16)
# axes.set_xlim(0,3)
# axes.set_ylabel("Count",fontsize=16)
# axes.legend(loc="best")
# plt.tight_layout()

#fig.savefig("../figures/figures_thesis_chapter_2/fig_distribution_of_normalized_insertions_for_annotated_essential_genes.png",dpi=400)

# -

def Free_index(gene,data,background):
    """ Computation of the free index parameter for a given gene,
    that represents the longest interval free of transposons over the 
    length of the gene. This is an important feature for the prediction of essentials
    genes. 

    Parameters
    ----------
    gene : str
        Gene of interest
    data : pd.DataFrame
        data containing the insertions locations per gene
    background : str
        genetic background of the library

    Returns
    -------
    _type_
        _description_
    """
    
    data_background=data.loc[background]
    data_background.index=data_background["Gene name"]

    fi=defaultdict(dict)
    insertion_float=from_excel_to_list(data_background.loc[gene,"Insertion locations"])     

    l=data_background.loc[gene,"End location"]-data_background.loc[gene,"Start location"]

    if (type(insertion_float)!=int) and (len(insertion_float)>1): 
        x=np.array(insertion_float)

        L=np.max(x[1:] - x[:-1])

        fi["free index"][gene]=L/l
    else:
        fi["free index"][gene]=0

    return fi


# +
backgrounds=["wt_merged","dnrp1_merged"]
fi_all=defaultdict(dict)
fi_all_per_background=[]
for background in backgrounds:
    for gene in list_data_pd.loc[background]["Gene name"]:
        fi=Free_index(gene,list_data_pd,background)
        # tmp=pd.DataFrame.from_dict(fi)
        fi_all["free index"][gene]=fi["free index"][gene]
    fi_all_per_background.append(pd.DataFrame(fi_all))
    


# -

fi_all_pd=pd.concat(fi_all_per_background,keys=backgrounds)
fi_all_pd.head(4)
fi_wt=fi_all_pd.loc["wt_merged"]
fi_wt

standard_essentials= [x for x in standard_essentials if pd.isnull(x) == False]

# +
fi_essentials=[]
fi_non_essentials=[]
for i in standard_essentials: 
    if i in fi_wt.index:
        tmp=fi_wt.loc[i][0]
        fi_essentials.append(tmp)
  
ref=np.median(fi_essentials) #reference for essentials , FI close to 1 is essential 

# +

fi_non_essentials=fi_wt.loc[~fi_wt.index.isin(standard_essentials)]


# +

plt.hist(fi_non_essentials,bins=100,alpha=0.5,color="black",label="non essential genes");
plt.hist(fi_essentials,bins=100,alpha=0.5,color="white",label="essential genes");
# -

plt.errorbar(0,np.median(fi_essentials),yerr=np.std(fi_essentials),fmt="o",color="red",label="median of essentials",
capsize=10)
plt.errorbar(0.5,np.median(fi_non_essentials),yerr=np.std(fi_non_essentials),fmt="o",color="black",label="median of non essentials",
capsize=10)
plt.legend(loc="best")

# +
fi_wt=fi_all_pd.loc["wt_merged"]
fi_ho=fi_wt.loc["HO"]

fi_essentials_norm=[]


for i in standard_essentials: 
    if i in fi_wt.index:
        tmp=fi_wt.loc[i][0]
        fi_essentials_norm.append(tmp/fi_ho)

n=1.8
fi_essentials_norm=np.concatenate(np.array(fi_essentials_norm))
fi_essentials_norm_above_ho=fi_essentials_norm[(fi_essentials_norm>np.array(n))]
d=len(fi_essentials_norm_above_ho)/len(standard_essentials)

## Distribution of FI for all genes compared to HO

values=fi_wt/fi_ho
values=values["free index"].values

values_above_ho=values[(values>np.array(n))]
d_all=len(values_above_ho)/len(values)

# +
fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

axes.hist(values,bins=100,alpha=0.5,color="gray",label="all genes");
axes.vlines(n,0,350,color="red",linestyle="dashed",linewidth=2,label="2*HO");
axes.annotate(f"{d_all*100:.2f}" +"%",xy=(3,300),fontsize=16)
axes.set_title("Distribution of relative FI for all genes",fontsize=16)
axes.tick_params(axis="both",labelsize=16)
axes.set_xlabel("Relative Free Index",fontsize=16)
# axes.set_xlim(0,3)
axes.set_ylabel("Count",fontsize=16)
axes.legend(loc="best",fontsize=16)


# These are in unitless percentages of the figure size. (0,0 is bottom left)
left, bottom, width, height = [0.5, 0.3, 0.4, 0.4]
ax2 = fig.add_axes([left, bottom, width, height])

data = [len(fi_essentials_norm_above_ho),len(standard_essentials)-len(fi_essentials_norm_above_ho)]
labels = ['Truly \n essentials'," "]

colors=["pink","gray"]



ax2.pie(data, labels = labels, colors = colors, autopct='%.0f%%',textprops={'fontsize': 14});

plt.tight_layout()

#fig.savefig("../figures/figures_thesis_chapter_2/fig_distribution_of_relative_FI_pie_inset_essentials.png",dpi=400)
# -

# ## Evaluating the FI ad NI as predictors for essentiality

# +
## Compute the ROC curve for the NI and FI to predict essentiality

### tasks:
#- Convert the FI and NI scores into probabilities scores to the true class:
# for example , normalize the FI scores from 0 to 1 and the NI scores from 0 to 1, the probality of being essential is 
# for the case of FI= normalization , and for the NI = 1-normalization.  
# z=xi-min(x)/max(x)-min(x)

## Example

from sklearn import metrics
import numpy as np
y = np.array([1, 1, 2, 2])
scores = np.array([0.1, 0.4, 0.35, 0.8])
fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=2)
area=metrics.auc(fpr,tpr)
plt.plot(fpr, tpr)
print(area)
# -

fi2roc_values=fi_wt.loc[:,"free index"].values/fi_ho.values
ni2roc_values=ni/ni_ho

fi2roc_values=fi_wt.loc[:,"free index"].values
ni2roc_values=ni

# +
fi_all_pd_2roc=(fi2roc_values-np.min(fi2roc_values))/(np.max(fi2roc_values)-np.min(fi2roc_values))
y_pd=write_ones_if_essential(fi_all_pd,"wt_merged",standard_essentials)
y=y_pd.loc[:,"true essential"].values

ni2roc=1-(ni2roc_values-np.min(ni2roc_values))/(np.max(ni2roc_values)-np.min(ni2roc_values)) # for the ni the highest the value the least probable is , that is why I need to substract that value to 1

# +
figure,ax=plt.subplots(nrows=1,ncols=2,figsize=(10,5))
plt.subplots_adjust(wspace=0.3)

ax[1].set_xlabel("Probability of being essential",fontsize=16)
ax[1].tick_params(axis="both",labelsize=16)
ax[1].set_ylabel("Counts",fontsize=16)
ax[1].hist(ni2roc[y==0],bins=300,alpha=0.5,label="Non essential genes",color="gray");
ax[1].hist(ni2roc[y==1],bins=80,label="Essential genes",color="pink");
ax[1].legend(fontsize=16)
ax[1].set_title("Neighborhood Index",fontsize=16)
ax[1].set_xlim(0.7,1)


ax[0].set_xlabel("Probability of being essential",fontsize=16)
ax[0].tick_params(axis="both",labelsize=16)
ax[0].set_ylabel("Counts",fontsize=16)
ax[0].hist(fi_all_pd_2roc[y==0],bins=100,alpha=0.5,label="Non essential genes",color="gray");
ax[0].hist(fi_all_pd_2roc[y==1],bins=100,label="Essential genes",color="pink");
ax[0].legend(fontsize=16)
ax[0].set_title("Free Index",fontsize=16)

plt.tight_layout()

figure.savefig("../figures/figures_thesis_chapter_2/fig_prob_distributions_fi_ni.png",dpi=400)


# +
from sklearn import metrics
import numpy as np

fpr, tpr, thresholds = metrics.roc_curve(y, fi_all_pd_2roc)
fpr_ni,tpr_ni,thresholds_ni=metrics.roc_curve(y,ni2roc)

area=metrics.auc(fpr,tpr)
area_ni=metrics.auc(fpr_ni,tpr_ni)

figure,ax=plt.subplots(nrows=1,ncols=2,figsize=(10,5))

plt.subplots_adjust(wspace=0.3)

ax[0].plot(fpr, tpr ,label=f"AUC={area:.2f}",color="darkorange",lw=2)
ax[0].set_title("ROC Free Index",fontsize=16)

ax[1].plot(fpr_ni,tpr_ni,color="darkorange",lw=2,label=f"AUC={area_ni:.2f}")  
ax[1].set_title("ROC Neighborhood Index",fontsize=16)

for axes in ax:

    axes.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
    axes.set_xlim([0.0, 1.0])
    axes.set_ylim([0.0, 1.05])
    axes.set_xlabel("False Positive Rate",fontsize=16)
    axes.set_ylabel("True Positive Rate",fontsize=16)
    
    axes.legend(loc="lower right",fontsize=16)
  
#figure.savefig("../figures/figures_thesis_chapter_2/fig_roc_curve_fi_ni.png",dpi=400)
# -



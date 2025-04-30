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

import scipy
from scipy import stats
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import pickle

# +
with open("../postprocessed-data/fitness_models_all_backgrounds_enzo_analysis_merged", "rb") as fp:   # Unpickling
    fitness_average = pickle.load(fp)

with open("../postprocessed-data/fitness_models_all_backgrounds_domain_analysis_merged", "rb") as fp:   # Unpickling
    fitness_domains = pickle.load(fp)

with open("../postprocessed-data/fitness_models_all_backgrounds_domain_analysis_corrected", "rb") as fp:   # Unpickling
    fitness_domains_corrected = pickle.load(fp)
standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str) #standard essentials
data_domains=pd.read_excel("../postprocessed-data/genomic-domains-wt.xlsx",index_col="Unnamed: 0")
## data from yeastmine
domains_names=pd.read_csv('../data/Domains_all_genes_protein_coordinates_yeastmine.tsv',sep="\t")
domains_names.index=domains_names["Gene Name"]

# +
# how many essential genes do not have domains annotated 

data_domains_essential_index=data_domains.index.intersection(standard_essentials)
data_domains_essential=data_domains.loc[data_domains_essential_index]

essential_without_domains=[]
essential_with_domains=[]
for i in data_domains_essential.index:
    if len(data_domains_essential.loc[i,"domain coordinates"]) == 2:
        essential_without_domains.append(i)
    else:
        essential_with_domains.append(i)
        


# +
## Determine the fitness for essential genes  using both fitness types 
background="wt_merged"

fitness_domains_corrected_essentials_index=fitness_domains_corrected.loc[background].index.intersection(essential_with_domains)
fitness_domainless_essentials_index=fitness_domains_corrected.loc[background].index.intersection(essential_without_domains)
fitness_average_essentials_index=fitness_average.loc[background].index.intersection(standard_essentials)

fitness_domains_corrected_non_essentials_index=fitness_domains_corrected.loc[background].index.difference(standard_essentials)
fitness_average_non_essentials_index=fitness_average.loc[background].index.difference(standard_essentials)


fitness_domains_corrected_essentials=fitness_domains_corrected.loc[background].loc[fitness_domains_corrected_essentials_index]
fitness_average_essentials=fitness_average.loc[background].loc[fitness_average_essentials_index]


fitness_domains_corrected_non_essentials=fitness_domains_corrected.loc[background].loc[fitness_domains_corrected_non_essentials_index]
fitness_average_non_essentials=fitness_average.loc[background].loc[fitness_average_non_essentials_index]

fitness_domainless_essentials=fitness_domains_corrected.loc[background].loc[fitness_domainless_essentials_index]

# +
## Plot histograms of each fitness highlighting essential genes 


figure, axes = plt.subplots(2, 1, figsize=(8, 5),sharex=True, sharey=True)

fitness_average_non_essentials["mean_norm"].hist(ax=axes[0], bins=100, color="grey", alpha=0.5, label="Non essential");
fitness_average_essentials["mean_norm"].hist(ax=axes[0], bins=100, color="red", alpha=0.5, label="Essential");

fitness_domains_corrected_non_essentials["fitness_domain_corrected"].hist(ax=axes[1], bins=100, color="grey", alpha=0.5, label="Non esential");
fitness_domains_corrected_essentials["fitness_domain_corrected"].hist(ax=axes[1], bins=100, color="red", alpha=0.5, label="Essential");
#fitness_domainless_essentials["fitness_domain_corrected"].hist(ax=axes[1], bins=100, color="blue", alpha=0.5, label="Essential without domains");

for i in axes: 
    i.set_ylabel("Number of genes")
    # i.legend(loc="upper left")
    i.grid(False)

axes[1].axvline(x=fitness_domains_corrected_non_essentials["fitness_domain_corrected"].median(), color="gray", linestyle="--")
axes[1].axvline(x=fitness_domains_corrected_essentials["fitness_domain_corrected"].median(), color="red", linestyle="--")
axes[0].axvline(x=fitness_average_non_essentials["mean_norm"].median(), color="gray", linestyle="--")
axes[0].axvline(x=fitness_average_essentials["mean_norm"].median(), color="red", linestyle="--")
axes[1].set_xlabel("Normalized Fitness")
axes[0].set_title("Average fitness")
axes[1].set_title("Fitness corrected by domain")

plt.tight_layout()
plt.savefig("../figures/fitness_domains_mean_non_essentials_vs_essentials.png", dpi=300)
# -

domains_names.loc["PAN1"]

# +
## std from both distributions 

fitness_average.loc["wt_merged"]["mean_norm"].std(), fitness_domains_corrected["fitness_domain_corrected"].std()

# +
## Values below zero from fitness corrected by domain

fitness_domains_corrected[fitness_domains_corrected["fitness_domain_corrected"]<=0].loc["wt_merged"]
# -

fitness_average[fitness_average["mean_norm"]<=0].loc["wt_merged"]

# +
## Fitness values above 1 from both methods 

fitness_domains_corrected[fitness_domains_corrected["fitness_domain_corrected"]>1].loc["wt_merged"]
# -

fitness_average[fitness_average["mean_norm"]>1].loc["wt_merged"]

fitness_domains_corrected_essentials[fitness_domains_corrected_essentials["fitness_domain_corrected"]<=0]

fitness_average_essentials[fitness_average_essentials["mean_norm"]<=0]

# +
## ROC curve to evaluate the performance of each method to predict existing essential genes 

# +
from sklearn.metrics import roc_curve, auc

data_wt_non_corrected=fitness_average.loc["wt_merged"]["mean_norm"]
data_wt_corrected=fitness_domains_corrected.loc["wt_merged"]["fitness_domain_corrected"]



essentials_pred=data_wt_non_corrected[data_wt_non_corrected<fitness_average_essentials["mean_norm"].median()].index.tolist()
y_corrected=np.zeros(len(data_wt_corrected)) ## labels for ROC curve
y_pred_corrected=np.zeros(len(data_wt_corrected)) ## predicted values for ROC curve
essentials_pred_corrected=data_wt_corrected[data_wt_corrected<fitness_domains_corrected_essentials["fitness_domain_corrected"].median()].index.tolist()
# -

fitness_average_essentials["mean_norm"].median(),fitness_domains_corrected_essentials["fitness_domain_corrected"].median()

# +
# assign 1 to the index that correspond to essential genes 

y=np.zeros(len(data_wt_non_corrected)) ## labels for ROC curve
y_pred=np.zeros(len(data_wt_non_corrected)) ## predicted values for ROC curve

j=0
for i in data_wt_non_corrected.index:
    if i in standard_essentials.tolist():
        y[j]=1
    j+=1
#
j=0
for i in data_wt_non_corrected.index:
    if i in essentials_pred:
        y_pred[j]=1
    j+=1


j=0
for i in data_wt_corrected.index:
    if i in standard_essentials.tolist():
        y_corrected[j]=1
    j+=1
#
j=0
for i in data_wt_corrected.index:
    if i in essentials_pred_corrected:
        y_pred_corrected[j]=1
    j+=1

# +
## to construc a confusiion matrix 
from sklearn.metrics import confusion_matrix
import seaborn as sns

# confusion_matrix(y, y_pred)/len(y)*100


cm = confusion_matrix(y, y_pred)

ax= plt.subplot()
sns.heatmap(cm, annot=True, ax = ax,fmt="d",cmap="Greys",annot_kws={"fontsize":17}); #annot=True to annotate cells

ax.set_xlabel('Predicted labels');ax.set_ylabel('True labels');
ax.set_title('Confusion Matrix');
ax.xaxis.set_ticklabels(['Non-essential', 'Essential']); ax.yaxis.set_ticklabels(['Non-essential', 'Essential']);

plt.savefig("../figures/confusion_matrix_average_fitness2predict_essentials.png", dpi=300)

# +
## to construc a confusiion matrix 
from sklearn.metrics import confusion_matrix

# confusion_matrix(y, y_pred)/len(y)*100


cm = confusion_matrix(y_corrected, y_pred_corrected)

ax= plt.subplot()
sns.heatmap(cm, annot=True, ax = ax,fmt="d",cmap="Greys",annot_kws={"fontsize":17}); #annot=True to annotate cells

ax.set_xlabel('Predicted labels');ax.set_ylabel('True labels');
ax.set_title('Confusion Matrix');
ax.xaxis.set_ticklabels(['Non-essential', 'Essential']); ax.yaxis.set_ticklabels(['Non-essential', 'Essential']);

plt.savefig("../figures/confusion_matrix_domain_corrected_fitness2predict_essentials.png", dpi=300)

# +
from sklearn import metrics

#fpr, tpr, thresholds = metrics.roc_curve(y, fitness2rocprob)
fpr, tpr, thresholds = metrics.roc_curve(y, y_pred)
area=metrics.auc(fpr,tpr)

figure,ax=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate",fontsize=16)
plt.ylabel("True Positive Rate",fontsize=16)
#plt.title("ROC of the  corrected fitness to predict essentiality",fontsize=16)

plt.plot(fpr, tpr ,label=f"AUC={area:.2f}",color="darkorange",lw=2)
ax.tick_params(axis="both",labelsize=16)
ax.legend(loc="lower right",fontsize=16)

figure.savefig("../figures/fig_non_corrected_fitness_ROC_curve.png",dpi=300)

# +
from sklearn import metrics

#fpr, tpr, thresholds = metrics.roc_curve(y, fitness2rocprob)
fpr, tpr, thresholds = metrics.roc_curve(y_corrected, y_pred_corrected)
area=metrics.auc(fpr,tpr)

figure,ax=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate",fontsize=16)
plt.ylabel("True Positive Rate",fontsize=16)
#plt.title("ROC of the  corrected fitness to predict essentiality",fontsize=16)

plt.plot(fpr, tpr ,label=f"AUC={area:.2f}",color="darkorange",lw=2)
ax.tick_params(axis="both",labelsize=16)
ax.legend(loc="lower right",fontsize=16)

figure.savefig("../figures/fig_corrected_fitness_ROC_curve.png",dpi=300)

# +
## Reliability of the fitness models 

wt_a=fitness_average.loc["wt_a"]["mean_norm"]
wt_b=fitness_average.loc["wt_b"]["mean_norm"]

wt_a_domains=fitness_domains_corrected.loc["wt_a"]["fitness_domain_corrected"]
wt_b_domains=fitness_domains_corrected.loc["wt_b"]["fitness_domain_corrected"]

M_1=fitness_average.loc["dnrp1_1"]["mean_norm"]
M_2=fitness_average.loc["dnrp1_2"]["mean_norm"]

M_1_domains=fitness_domains_corrected.loc["dnrp1_1"]["fitness_domain_corrected"]
M_2_domains=fitness_domains_corrected.loc["dnrp1_2"]["fitness_domain_corrected"]


# +
# Pearsom correlation coefficient
### remove nan 
wt_a=wt_a.dropna()
wt_b=wt_b.dropna()

wt_a_domains=wt_a_domains.dropna()
wt_b_domains=wt_b_domains.dropna()


M_1=M_1.dropna()
M_2=M_2.dropna()

M_1_domains=M_1_domains.dropna()
M_2_domains=M_2_domains.dropna()

## take the intersection of the indexes
wt_common=wt_a.index.intersection(wt_b.index)

wt_domains_common=wt_a_domains.index.intersection(wt_b_domains.index)

M_common=M_1.index.intersection(M_2.index)

M_domains_common=M_1_domains.index.intersection(M_2_domains.index)

# -

scipy.stats.pearsonr(wt_a_domains[wt_domains_common],wt_b_domains[wt_domains_common])

scipy.stats.pearsonr(wt_a.loc[wt_common],wt_b.loc[wt_common])

# +
fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

axes.scatter(wt_a.loc[wt_common],wt_b.loc[wt_common],color="grey",alpha=0.4,label="WT")
# axes.vlines(x=1,ymin=-1.5,ymax=2,color="black",linestyle="--")
# axes.hlines(y=1,xmin=-1.5,xmax=2,color="black",linestyle="--")

# plt.plot([-1.5,2.0],[-1.5,2.0],color="black",linestyle="--")

plt.xlabel("Fitness WT1_a",fontsize=16)
plt.ylabel("Fitness WT1_b",fontsize=16)

plt.text(0.1, 0.9, f"$\sigma$={scipy.stats.pearsonr(wt_a.loc[wt_common],wt_b.loc[wt_common])[0]:.2f}", horizontalalignment='center',verticalalignment='center', transform=axes.transAxes,fontsize=16)
plt.ylim(-1.5,2.5)
plt.xlim(-1.5,2.5)
plt.yticks([-1,0,1,2])
plt.xticks([-1,0,1,2])
plt.tight_layout()
plt.savefig("../figures/fig_fitness_reliability_average_technical_replicates.png",dpi=300)

# +
fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

axes.scatter(wt_a_domains[wt_domains_common],wt_b_domains[wt_domains_common],color="grey",alpha=0.4,label="WT")
# axes.vlines(x=1,ymin=-1.5,ymax=2,color="black",linestyle="--")
# axes.hlines(y=1,xmin=-1.5,xmax=2,color="black",linestyle="--")

# plt.plot([-1.5,2.0],[-1.5,2.0],color="black",linestyle="--")

plt.xlabel("Fitness WT1_a",fontsize=16)
plt.ylabel("Fitness WT1_b",fontsize=16)

plt.text(0.1, 0.9, f"$\sigma$={scipy.stats.pearsonr(wt_a_domains[wt_domains_common],wt_b_domains[wt_domains_common])[0]:.2f}", horizontalalignment='center',verticalalignment='center', transform=axes.transAxes,fontsize=16)
plt.ylim(-1.5,2.5)
plt.xlim(-1.5,2.5)
plt.yticks([-1,0,1,2])
plt.xticks([-1,0,1,2])
plt.tight_layout()
plt.savefig("../figures/fig_fitness_reliability_domains_technical_replicates.png",dpi=300)

# +
fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

axes.scatter(M_1.loc[M_common],M_2.loc[M_common],color="grey",alpha=0.4,label="Mutant")

plt.xlabel("Fitness M1",fontsize=16)
plt.ylabel("Fitness M2",fontsize=16)

plt.text(0.1, 0.9, f"$\sigma$={scipy.stats.pearsonr(M_1.loc[M_common],M_2.loc[M_common])[0]:.2f}", horizontalalignment='center',verticalalignment='center', transform=axes.transAxes,fontsize=16)

plt.ylim(-1.5,2.5)
plt.xlim(-1.5,2.5)
plt.yticks([-1,0,1,2])
plt.xticks([-1,0,1,2])
plt.tight_layout()

plt.savefig("../figures/fig_fitness_reliability_average_biological_replicates.png",dpi=300)

# +
fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

axes.scatter(M_1_domains.loc[M_domains_common],M_2_domains.loc[M_domains_common],
             color="grey",alpha=0.4,label="Mutant")

plt.xlabel("Fitness M1",fontsize=16)
plt.ylabel("Fitness M2",fontsize=16)

plt.text(0.1, 0.9, f"$\sigma$={scipy.stats.pearsonr(M_1_domains.loc[M_domains_common],M_2_domains.loc[M_domains_common])[0]:.2f}", horizontalalignment='center',verticalalignment='center', transform=axes.transAxes,fontsize=16)
plt.ylim(-1.5,2.5)
plt.xlim(-1.5,2.5)
plt.yticks([-1,0,1,2])
plt.xticks([-1,0,1,2])
plt.tight_layout()
plt.savefig("../figures/fig_fitness_reliability_domains_biological_replicates.png",dpi=300)

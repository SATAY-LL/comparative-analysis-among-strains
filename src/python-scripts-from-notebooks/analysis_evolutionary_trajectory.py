# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3.8.10 ('satay-dev')
#     language: python
#     name: python3
# ---

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
import scipy.stats as stats
import seaborn as sns
from collections import defaultdict

## standard for plots
plt.rc('font', family='serif',size=14)
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)


# +
## Plot the number of genes with less than 5 transposons in all mutants 

pergene_files=[]
#data_dir= "../satay/data_files/data_unmerged/"
#data_dir="../transposonmapper/data_files/files4test/"
data_dir="../postprocessed-data/"
#data_dir="../transposonmapper/data_files/"
for root, dirs, files in os.walk(data_dir):
    for file in files:
        if file.endswith("pergene_insertions.xlsx"):
            pergene_files.append(os.path.join(root, file))
# -

list_data=[]
for i in pergene_files:
    list_data.append(pd.read_excel(i,engine='openpyxl',index_col="Unnamed: 0"))

keys=[]
for i in np.arange(0,len(pergene_files)):
    keys.append(pergene_files[i].split("/")[-1].split("_")[0]+"_"+pergene_files[i].split("/")[-1].split("_")[1])

# import essential genes used in transposonmapper
standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str)


keys

backgrounds=["wt_a","wt_b",'wt_merged', 'bem1-aid_a', 'bem1-aid_b', "bem1-aid_merged",
 'dbem1dbem3_a', 'dbem1dbem3_b','dnrp1_merged', 'dbem3_merged']

backgrounds=keys

list_data_pd=pd.concat(list_data,axis=0,keys=keys)

list_data_pd=list_data_pd.loc[backgrounds]

# +
# import excel file with the normalized data

# data_norm_pd=pd.read_excel("../postprocessed-data/data_norm_linear_transformation_per_background.xlsx",
# engine='openpyxl',index_col="background")
# data_norm_pd.drop(columns=["Unnamed: 0","Unnamed: 1"],inplace=True)

# +
# data_norm_pd.index.unique()
# -

polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")
polarity_genes.fillna(0,inplace=True)

# +
# data_norm_pd_wt=data_norm_pd.loc["wt_merged"]
# data_norm_pd_wt.columns
# x_norm=data_norm_pd_wt.loc[:,"reads_normalized_windows"]/data_norm_pd_wt.loc[:,"tr_normalized_windows"]
# x=data_norm_pd_wt.loc[:,"Reads"]/data_norm_pd_wt.loc[:,"Insertions"]
# # plt.subplot(1,2,1)
# # plt.hist(x_norm,bins=100);
# # plt.subplot(1,2,2)
# # plt.hist(x,bins=100);
# x_norm.std()/x_norm.mean(),x.std()/x.mean()

# +
# Number of genes with less than 2 transposons and 2 reads in all mutants
Q=[]
for i in np.arange(0,len(backgrounds)):
    tmp=(list_data_pd.loc[backgrounds[i]])
    #L.append(len(tmp[(tmp.loc[:,"Insertions"]<5) & (tmp.loc[:,"Reads"]<2)]))
    Q.append(len(tmp[(tmp.loc[:,"Insertions"]<5)]))



# +
Q_dict=dict(zip(backgrounds,Q))

A={k: v for k, v in sorted(Q_dict.items(), key=lambda item: item[1])}

fig = plt.figure(figsize=(10, 7))
plt.bar(A.keys(),A.values())
plt.xticks(rotation=90);

plt.ylabel("# of genes with less than 5 insertions ",fontsize=9)
plt.tight_layout(pad=3)
#plt.savefig("../figures/fig_number_genes_with_less_than_5_reads.png",dpi=300)



# +
# Number of transposons in each library 
L=[]
R=[]
for i in np.arange(0,len(backgrounds)):
    tmp=(list_data_pd.loc[backgrounds[i]])
    L.append(np.sum(tmp.loc[:,"Insertions"]))
    R.append(np.sum(tmp.loc[:,"Reads"]))

    
L_dict=dict(zip(backgrounds,L))
R_dict=dict(zip(backgrounds,R))

A={k: v for k, v in sorted(L_dict.items(), key=lambda item: item[1])}
B={k: v for k, v in sorted(R_dict.items(), key=lambda item: item[1])}

fig,axes = plt.subplots(nrows=1,ncols=2,figsize=(15, 7))
axes[0].bar(A.keys(),A.values())
axes[0].set_xticklabels(labels=A.keys(),rotation=90);
axes[0].set_ylabel("# transposons per library",fontsize=9)

axes[1].bar(B.keys(),B.values())
axes[1].set_xticklabels(labels=B.keys(),rotation=90);
axes[1].set_ylabel("# reads per library",fontsize=9)

#plt.savefig("../figures/fig_number_insertions_and_reads_for_all_libraries.png",dpi=300)

plt.tight_layout(pad=3)

# +
## Plot length of genes vs number of transposons
from scipy.optimize import curve_fit

def func(x, a,b):
    return a*x +b

fig,axes = plt.subplots(nrows=1,ncols=1,figsize=(8,5))

data=list_data_pd.loc["wt_merged"]
xaxis=(data.loc[:,"End location"]-data.loc[:,"Start location"])
sns.regplot(x=xaxis,y=data.loc[:,"Insertions"],ax=axes,color="gray",scatter_kws={"s":80,"alpha":0.2},label="WT")

# plt.plot(xaxis,data.loc[:,"Insertions"],"o",markersize=10,alpha=0.5)

# x = np.array(xaxis)
# popt, pcov = curve_fit(func, x, np.array(data.loc[:,"Insertions"]))
# plt.plot(x, func(x, *popt), 'ko', label="Fitted Curve->" + str(popt[0]) + '*x+'+str(popt[1]),alpha=0.6)

plt.ylabel("gene length",fontsize=18)
# plt.xscale("log")
# plt.yscale("log")
plt.xlabel("# transposons",fontsize=18)
plt.legend()
#plt.savefig("../figures/figures_thesis_chapter_2/supp_fig_length_vs_number_of_transposons_linear.png",dpi=300)
# -



# +
## Plot the number of genes with less than 5 transposons in all mutants over the total number of transposons per library
from scipy.optimize import curve_fit

keys2plot=["wt_a","wt_b", 'bem1-aid_a', 'bem1-aid_b', 'dbem1dbem3_a', 'dbem1dbem3_b']
L=[]
R=[]
for i in np.arange(0,len(keys2plot)):
    tmp=(list_data_pd.loc[keys2plot[i]])
    L.append(np.sum(tmp.loc[:,"Insertions"]))
    R.append(np.sum(tmp.loc[:,"Reads"]))

    
L_dict=dict(zip(keys2plot,L))
R_dict=dict(zip(keys2plot,R))

Q=[]
for i in np.arange(0,len(backgrounds)):
    tmp=(list_data_pd.loc[backgrounds[i]])
    #L.append(len(tmp[(tmp.loc[:,"Insertions"]<5) & (tmp.loc[:,"Reads"]<2)]))
    Q.append(len(tmp[(tmp.loc[:,"Insertions"]<5)]))

Q_dict=dict(zip(backgrounds,Q))

fig,axes = plt.subplots(nrows=1,ncols=1,figsize=(10, 10))
l_array=[]
q_array=[]
for keys in L_dict.keys():
    axes.plot(L_dict[keys],Q_dict[keys],'*r',markersize=20,alpha=0.5)
    axes.annotate(keys,(L_dict[keys],Q_dict[keys]))
    l_array.append(L_dict[keys])
    q_array.append(Q_dict[keys])

# def func(x, a,b):
#     return b+a / x

# x = np.array(l_array)

# popt, pcov = curve_fit(func, x, np.array(q_array))
# axes.plot(x, func(x, *popt), 'ko', label="Fitted Curve (b+a/x)",alpha=0.6,markersize=10)

axes.set_xlabel("# of insertions in that library",fontsize=16)
axes.set_ylabel("# of genes with less than 5 insertions",fontdict={"size":16})
#axes.legend()
# axes.set_xscale("log")
axes.set_yscale("log")
plt.grid(axis="both",which="both",ls="--",alpha=0.5)
plt.tick_params(axis='both', which='both', labelsize=16)

fig.savefig("../figures/number_of_insertions_less_than_5_vs_total_number_of_reads.png",dpi=300)

# +
# from module_intergenic_model import adding_features2dataframe,getting_r

# list_data_extended=[]
# for i in np.arange(0,len(backgrounds)):
#     tmp=(list_data_pd.loc[backgrounds[i]])
#     list_data_extended.append(adding_features2dataframe(tmp))



# +
# list_data_extended_pd=pd.concat(list_data_extended,axis=0,keys=backgrounds)

# +
# Number of reads per transposons per library
# L=[]
# R=[]
# for i in np.arange(0,len(backgrounds)):
#     tmp=(list_data_extended_pd.loc[backgrounds[i]])
#     L.append(tmp.loc[:,"reads-per-tr"].mean())
#     R.append(tmp.loc[:,"reads-per-tr"].std())

# L_dict=dict(zip(backgrounds,L))
# R_dict=dict(zip(backgrounds,R))

# A={k: v for k, v in sorted(L_dict.items(), key=lambda item: item[1])}
# B={k: v for k, v in sorted(R_dict.items(), key=lambda item: item[1])}

# fig,axes = plt.subplots(nrows=1,ncols=2,figsize=(15, 7))
# axes[0].bar(A.keys(),A.values())
# axes[0].set_xticklabels(labels=A.keys(),rotation=90);
# axes[0].set_ylabel("Mean number of reads per transposons per library",fontsize=9)

# axes[1].bar(B.keys(),B.values())
# axes[1].set_xticklabels(labels=B.keys(),rotation=90);
# axes[1].set_ylabel("Standard deviation of number of reads per transposons per library",fontsize=9)

# plt.tight_layout(pad=3)
#plt.savefig("../figures/fig_mean_and_std_number_reads_per_transposons_per_library.png",dpi=300)
# -

keys=['dbem3_b',
 'WT_1-Benoit',
 'dnrp1_2',
 'bem1-aid_a',
 'dbem1dbem3_b',
 'wt_merged',
 'dbem3_a-trimmed',
 'dbem1dbem3_a',
 'bem1-aid-dbem3_a',
 'bem1-aid_merged',
 'bem1-aid-dbem3_b',
 'dnrp1_1',
 'wt_b',
 'wt_a',
 'dnrp1_merged',
 'WT_2-Benoit',
 'bem1-aid_b',
 'dbem3_merged',
 'dbem3_a']

# +
#######  Import fitness from SATAY #########
import pickle
with open("../postprocessed-data/fitness_models_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

satay_fitness=pd.concat(b,axis=0,keys=keys)

standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str)

# +
satay_wt=satay_fitness.loc["wt_merged"]


satay_wt2compare=satay_wt.loc[:,"fitness_domains_corrected"]
#satay_wt2compare=satay_wt.loc[:,"fitness_gene"]


satay_wt2compare=pd.DataFrame(satay_wt2compare)
satay_wt2compare.columns=["fitness"]

satay_wt2compare.index.name="Gene"
satay_wt2compare.dropna(inplace=True)

satay_wt2compare=satay_wt2compare[satay_wt2compare.loc[:,"fitness"]!= "Not enough flanking regions"]

## asign zero to the fitness that are not enough reads or insertions

satay_wt2compare.loc[satay_wt2compare.loc[:,"fitness"]=="Not enough reads","fitness"]=0
satay_wt2compare.loc[satay_wt2compare.loc[:,"fitness"]=="Not enough insertions","fitness"]=0

# satay_wt2compare=satay_wt2compare[satay_wt2compare.loc[:,"fitness"]!= "Not enough reads"]
# satay_wt2compare=satay_wt2compare[satay_wt2compare.loc[:,"fitness"]!= "Not enough insertions"]



# +
essential_satay=satay_wt2compare[satay_wt2compare.loc[:,"fitness"]<0.5]
len(essential_satay),len(standard_essentials)

# how many essentials for satay are also essentials in S288C 

common_essentials=len(np.intersect1d(essential_satay.index,standard_essentials))

# essentials that are unique for satay_fitness

unique_essentials_satay=np.setdiff1d(essential_satay.index,standard_essentials)

# essentials that are unique for S288C

unique_essentials_S288C=np.setdiff1d(standard_essentials,essential_satay.index)

print("common_essentials",common_essentials/len(standard_essentials)*100,"% for reference")
print("common_essentials",common_essentials/len(essential_satay)*100,"%")
print("unique_essentials_satay",len(unique_essentials_satay)/len(essential_satay)*100,"%")
print("unique_essentials_S288C",len(unique_essentials_S288C)/len(standard_essentials)*100,"%")


# +

values=satay_wt2compare.astype(float).values

fig, axes = plt.subplots(2, 1,  gridspec_kw={"height_ratios":(.10, .30)}, figsize = (8, 5))


sns.violinplot(values, ax=axes[0],color="gray",orient="h",inner="quartile")

g=sns.histplot(values,bins=200,color="gray",ax=axes[1],stat="percent",label="WT",kde=True,element="step")

axes[1].set_xlabel("Corrected fitness values",fontsize=16)
axes[1].set_ylabel("Percent",fontsize=16)
axes[1].tick_params(axis="both",labelsize=16)
axes[0].tick_params(axis="both",labelsize=16)

#axes[1].vlines(x=1,ymin=0,ymax=2.5,color="red",linestyle="--",linewidth=2,label="HO")
#axes[1].vlines(x=np.mean(values),ymin=0,ymax=2.5,color="black",linestyle="--",linewidth=1)

axes[1].legend()
#fig.savefig("../figures/figures_thesis_chapter_2/fig_fitness_distribution_normalized_wt2HO.png",dpi=400)

# +
high=[1.1,np.max(values)]
low=[high[0]*0.5,high[0]*0.7]
neutral=[low[1],high[0]]
essential=low[0]

high,low,neutral,essential

# +
## Normalization of values from 0 to 1

valuesnorm=(values-np.min(values))/(np.max(values)-np.min(values))


# +
## make a pie chart with the mean fitness values compared to HO

high=[1.1,np.max(values)]
low=[high[0]*0.5,high[0]*0.7]
neutral=[low[1],high[0]]
essential=low[0]

values_neutral=values[np.where((values>neutral[0]) & (values<neutral[1]))]
values_high=values[np.where((values>high[0]) & (values<high[1]))]
values_low=values[np.where((values<low[1]) & (values>low[0]))]
values_essential=values[np.where(values<essential)]


# high=[0.8,1]
# neutral=[0.65,0.8]
# low=[0.3,0.65]
# essential=[0,0.3]

# values_neutral=valuesnorm[np.where((valuesnorm>neutral[0]) & (valuesnorm<neutral[1]))]
# values_high=valuesnorm[np.where((valuesnorm>high[0]) & (valuesnorm<high[1]))]
# values_low=valuesnorm[np.where((valuesnorm<low[1]) & (valuesnorm>low[0]))]
# values_essential=valuesnorm[np.where(valuesnorm<essential[1])]


fig,axes=plt.subplots(1,1,figsize=(8,8))
#define data
data = [len(values_neutral),len(values_high),len(values_low),len(values_essential)]
labels = ['Neutral', 'Advantageous', 'Disadvantageous', 'Possible essentials']

#define Seaborn color palette to use
#colors = sns.color_palette('pastel')[0:4]
colors=["blue","green","red","pink"]

#create pie chart
plt.pie(data, labels = labels, colors = colors, autopct='%.0f%%',textprops={'fontsize': 18});
plt.tight_layout()
#fig.savefig("../figures/figures_thesis_chapter_2/fig_fitness_pie_normalized_wt2HO.png",dpi=400)
# -

# ### ROC curve to assess the performance of the low fitness values with essentiality.
#

# +

# fitness2roc=1-fitness
# fitness2rocprob=(fitness2roc-np.min(fitness2roc))/(np.max(fitness2roc)-np.min(fitness2roc))



# # for i in np.arange(0,len(standard_essentials)):
# #     if standard_essentials[i] in satay_wt2compare.index:
# #         x=np.where(satay_wt2compare.index==standard_essentials[i])[0][0]
# #         y[x]=1 # write 1 in the position of essential genes
# #         print(x)

fitness=satay_wt2compare.loc[:,"fitness"].astype(float)
y=np.zeros(len(fitness)) ## labels for ROC curve
y_pred=np.zeros(len(fitness)) ## predicted values for ROC curve
j=0
for i in satay_wt2compare.index:
    if i in standard_essentials.tolist():
        y[j]=1
   
        
    j+=1
#
j=0
for i in satay_wt2compare.index:
    if i in essential_satay.index.tolist():
        y_pred[j]=1
        
   
        
    j+=1





# +
## to construc a confusiion matrix 
from sklearn.metrics import confusion_matrix

# confusion_matrix(y, y_pred)/len(y)*100

# plot a heatmap of the confusion matrix
import seaborn as sns
import matplotlib.pyplot as plt

cm = confusion_matrix(y, y_pred)

ax= plt.subplot()
sns.heatmap(cm, annot=True, ax = ax,fmt="d",cmap="Blues"); #annot=True to annotate cells

# labels, title and ticks
# set the labels in bold 

ax.set_xlabel('Predicted labels');ax.set_ylabel('True labels');
ax.set_title('Confusion Matrix');
ax.xaxis.set_ticklabels(['Non-essential', 'Essential']); ax.yaxis.set_ticklabels(['Non-essential', 'Essential']);

#plt.savefig("../figures/fig_confusion_matrix_corrected_fitness2essentiality.png",dpi=400)
# -

print(4247/(len(y)-len(standard_essentials))*100,  "percentage of non-essential genes")
print(474/len(y)*100, "percentage of false negatives")
print(1229/len(y)*100, "percentage of false positives")
print(636/len(standard_essentials)*100, "percentage of true positives")

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
plt.title("ROC of the  corrected fitness to predict essentiality",fontsize=16)

plt.plot(fpr, tpr ,label=f"AUC={area:.2f}",color="darkorange",lw=2)
ax.tick_params(axis="both",labelsize=16)
ax.legend(loc="lower right",fontsize=16)

#figure.savefig("../figures/fig_corrected_fitness_ROC_curve.png",dpi=300)

# +
## subplots of the confusion matrix and ROC curve

figure,ax=plt.subplots(nrows=1,ncols=2,figsize=(10,4))

plt.subplots_adjust(wspace=0.5)
plt.subplot(1,2,1)
# plot a ROC curve

plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate",fontsize=16)
plt.ylabel("True Positive Rate",fontsize=16)
plt.title("ROC curve",fontsize=16)

plt.plot(fpr, tpr ,label=f"AUC={area:.2f}",color="darkorange",lw=2)
ax[0].tick_params(axis="both",labelsize=14)
ax[0].legend(loc="lower right",fontsize=14)


# plot a heatmap of the confusion matrix
plt.subplot(1,2,2)

cm = confusion_matrix(y, y_pred)

## normalize the confusion matrix
cm = cm.astype('float') / cm.sum(axis=0)

sns.heatmap(cm, annot=True, ax = ax[1],cmap="Blues"); #annot=True to annotate cells

# labels, title and ticks
# set the labels in bold 

ax[1].set_xlabel('Predicted labels');ax[1].set_ylabel('True labels');
ax[1].set_title('Confusion Matrix');
ax[1].xaxis.set_ticklabels(['Non-essential', 'Essential']); ax[1].yaxis.set_ticklabels(['Non-essential', 'Essential']);
plt.tight_layout()

plt.savefig("../figures/fig_corrected_fitness_ROC_confusion_matrix.png",dpi=400)
# -



# ## Heat map polarity genes 

# +
## heatmap of fitnes values across backgrounds
#list_data_pd_wt=list_data_pd.loc["wt_merged"]

index_polarity_genes=[]
chosen_polarity_genes=[]
for i in np.arange(0,len(polarity_genes)):
    tmp=np.where(satay_wt2compare.index==polarity_genes.index[i])
    if tmp[0].size!=0:
        index_polarity_genes.append(tmp[0][0])
        chosen_polarity_genes.append(polarity_genes.index[i])
# -

backgrounds_heatmap=["wt_a","wt_b",'wt_merged', 'dnrp1_1','dnrp1_2','dnrp1_merged', 'dbem3_a',
'dbem3_b','dbem3_merged','bem1-aid_a', 'bem1-aid_b', "bem1-aid_merged",'dbem1dbem3_a', 'dbem1dbem3_b']

# +

satay_fitness_corrected=satay_fitness[satay_fitness.loc[:,"fitness_domains_corrected"]!= "Not enough flanking regions"]
satay_fitness_corrected=satay_fitness_corrected[satay_fitness_corrected.loc[:,"fitness_domains_corrected"]!= "Not enough reads"]
satay_fitness_corrected=satay_fitness_corrected[satay_fitness_corrected.loc[:,"fitness_domains_corrected"]!= "Not enough insertions"]



# +
array2heatmap=np.zeros((len(chosen_polarity_genes),len(backgrounds_heatmap)))

for i in np.arange(0,len(backgrounds_heatmap)):
    for j in np.arange(0,len(chosen_polarity_genes)):
        #array2heatmap[j,i]=rates_norm_dict_pd.loc[backgrounds_heatmap[i]].tolist()[0][index_polarity_genes[j]]
        tmp=satay_fitness_corrected.loc[backgrounds_heatmap[i]]
        if polarity_genes.index[j] in tmp.index:
            array2heatmap[j,i]=tmp.loc[polarity_genes.index[j]]["fitness_domains_corrected"]
        else:
            array2heatmap[j,i]=np.NaN #np.nan
  

# +
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 15))
# array2heatmap[~np.isfinite(array2heatmap)] = 0
sns.heatmap(array2heatmap,cmap="seismic",vmin=0,vmax=1,xticklabels=backgrounds_heatmap,
yticklabels=chosen_polarity_genes ,cbar=True,annot=True,cbar_kws={'label': 'Fitness reative to neutral genes'})

labels=[]
for i in np.arange(0,len(chosen_polarity_genes)):
    labels.append(chosen_polarity_genes[i]+"-"+str(polarity_genes.loc[chosen_polarity_genes[i],:].unique()))
ax.set_yticklabels(labels);
#fig.savefig("../figures/fig_heatmap_fitness_normalized_wt_with_replicates-and-function.png",dpi=300)
# -

array2heatmap[~np.isfinite(array2heatmap)]=2
g=sns.clustermap(array2heatmap,cmap="seismic",vmin=0,vmax=2,xticklabels=backgrounds_heatmap,
yticklabels=chosen_polarity_genes ,cbar=True,annot=True,cbar_kws={'label': 'Fitness'})
g.fig.set_size_inches((15,15))
#g.savefig("../figures/fig_heatmap_fitness_normalized2WT_coarse_grained_model.png",dpi=400)

# +
import math

def sigmoid(x,k,n):
    a = []
    for item in x:
        a.append(n/(n+math.exp(-item/k)))
    return a

import matplotlib.pyplot as plt
import numpy as np

x = np.arange(-8, 8, 0.2)
for k in np.arange(1,2,0.2):
    plt.plot(x, sigmoid(x,k,1), color='black', linewidth=1.0, linestyle='-')
    plt.xticks([])
    plt.yticks([])

#plt.savefig("../figures/figures_thesis_chapter_2/fig_sigmoid_function.png",dpi=400)


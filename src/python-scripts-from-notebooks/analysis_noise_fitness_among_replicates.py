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
## Importing fitness values 

fitness_all = pd.read_excel('../postprocessed-data/fitness_coarse_grained_all_pd.xlsx',
engine='openpyxl',index_col="Unnamed: 0")

# +
backgrounds=fitness_all["background"].unique()

fitness_wt=fitness_all[fitness_all["background"]=="wt_merged"]
normalized_fitness=[]
for i in backgrounds:
    fitness_all_i=fitness_all[fitness_all["background"]==i]
    index_i=fitness_all_i.index
    normalized_fitness=[]
    for j in np.arange(0,len(fitness_all_i)): 

        if fitness_wt.iloc[j,2]!=0:
            normalized_fitness.append(fitness_all_i.iloc[j,2]/fitness_wt.iloc[j,2])
        else:
            normalized_fitness.append(fitness_all_i.iloc[j,2])
    
    fitness_all.loc[index_i,"normalized_fitness"]=normalized_fitness
    


# -

backgrounds

# +
## Viz of differences between replicates, compared to HO

replicate_1="bem1-aid_a"
replicate_2="bem1-aid_b"

value_a=fitness_all[fitness_all["background"]==replicate_1]
value_b=fitness_all[fitness_all["background"]==replicate_2]

ho=np.where(fitness_all.loc[:,"Gene name"]=="HO")[0][0]

value_a_ho=value_a.iloc[ho,2] # fitness value of HO in replicate_a
value_b_ho=value_b.iloc[ho,2] # fitness value of HO in replicate_b

value_a_norm=[]
for i in value_a.iloc[:,2]:
    if value_a_ho!=np.inf and value_a_ho!=0 and value_a_ho!=np.nan and value_a_ho!=-np.inf:
        value_a_norm.append(i/value_a_ho)
    else:
        value_a_norm.append(i)

value_b_norm=[]
for i in value_b.iloc[:,2]:
    if value_b_ho!=np.inf and value_b_ho!=0 and value_b_ho!=np.nan and value_b_ho!=-np.inf:
        value_b_norm.append(i/value_b_ho)
    else:
        value_b_norm.append(i)



value_a_norm=np.array(value_a_norm)
value_b_norm=np.array(value_b_norm)
values_diff=value_a_norm-value_b_norm

values_diff=values_diff[np.where(values_diff!=-np.inf)]
values_diff=values_diff[np.where(values_diff!=np.inf)]
values_diff=values_diff[np.where(values_diff!=np.nan)]





# +
fig, axes = plt.subplots(3, 1,  gridspec_kw={"height_ratios":(.10, .40,.40)}, figsize = (8, 8))


sns.violinplot(np.abs(values_diff), ax=axes[0],color="gray",orient="h",inner="quartile")

g=sns.histplot(np.abs(values_diff),bins=200,color="gray",ax=axes[1],stat="percent",label="wt_a-wt_b",kde=True,element="step")

axes[1].set_xlabel("Fitness values of " + replicate_1 +"-" + replicate_2+ " compared to HO locus",fontsize=16)
axes[1].set_ylabel("Percent",fontsize=16)


sns.regplot(value_a_norm,value_b_norm,ax=axes[2],scatter_kws={"s":30,"color":"gray"})
#axes[1].vlines(x=1,ymin=0,ymax=2.5,color="red",linestyle="--",linewidth=2,label="HO")
#axes[1].vlines(x=np.mean(values),ymin=0,ymax=2.5,color="black",linestyle="--",linewidth=1)

# axes[2].set_ylim(np.min(value_b_norm),np.max(value_b_norm))
# axes[2].set_xlim(np.min(value_a_norm),np.max(value_a_norm))

axes[2].set_xlabel("Fitness values of "+replicate_1,fontsize=16)
axes[2].set_ylabel("Fitness values of "+ replicate_2,fontsize=16)

for axes in axes:
    axes.tick_params(axis="both",labelsize=16)

plt.tight_layout()
fig.savefig("../figures/figures_thesis_chapter_2/supp_fig_fitness_differences_"+replicate_1+"_"+replicate_2+".png",dpi=400)
    


# +
half_maximun=np.max(g.get_lines()[0].get_data()[1])/2

index_half_maximun=[0,41] # index of the half maximmun value in the histogram

x0,x1=g.get_lines()[0].get_data()[0][0],g.get_lines()[0].get_data()[0][41]
xf=g.get_lines()[0].get_data()[0][-1]

print("The full width at half maximum is:",x1-x0, 
"which represents the", (x1-x0)/xf*100,"% of the total width")
# -

half_maximun, g.get_lines()[0].get_data()[1]



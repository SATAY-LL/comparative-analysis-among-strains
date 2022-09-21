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

# +
# Plot a histogram of the fitness differences between both methods 

## Normalize the SGA values to the value of HO in the SGA dataset (HO is not in the SGA dataset)
# -

import pandas as pd
sga_fitness=pd.read_excel('../data/fitness-SGD.xlsx',index_col=0)
satay_fitness=pd.read_excel('../postprocessed-data/fitness_coarse_grained_all_pd.xlsx',index_col=0)


# +
satay_wt=satay_fitness[satay_fitness.loc[:,"background"]=="wt_merged"]
satay_wt.index=satay_wt.loc[:,"Gene name"]

#satay_wt_ho=satay_wt[satay_wt.loc[:,"Gene name"]=="HO"]
satay_wt_ho=satay_wt.loc["HO"]

#values2ho=satay_wt["fitness"]/satay_wt_ho["fitness"].values[0]
# -

satay_wt["fitness2HO"]=satay_wt["fitness"]/satay_wt_ho["fitness"]

satay_wt.loc["VPS8","fitness2HO"],len(sga_fitness.loc["vps8","fitness"]),type(sga_fitness.loc[i,"fitness"])

# +
from collections import defaultdict
import numpy as np

diff_sga_satay=defaultdict(dict)
for i in sga_fitness.index:
    if i.upper() in satay_wt.index:
        if type(sga_fitness.loc[i,"fitness"])!=np.float64:
            diff_sga_satay[i]["difference"]=np.abs(sga_fitness.loc[i,"fitness"].unique()[0]-satay_wt.loc[i.upper(),"fitness2HO"])
            diff_sga_satay[i]["sga"]=sga_fitness.loc[i,"fitness"].unique()[0]
            diff_sga_satay[i]["satay"]=satay_wt.loc[i.upper(),"fitness2HO"]
        else:
            diff_sga_satay[i]["difference"]=np.abs(sga_fitness.loc[i,"fitness"]-satay_wt.loc[i.upper(),"fitness2HO"])
            diff_sga_satay[i]["sga"]=sga_fitness.loc[i,"fitness"]
            diff_sga_satay[i]["satay"]=satay_wt.loc[i.upper(),"fitness2HO"]
    else:
        diff_sga_satay[i]["difference"]=np.nan
# -

diff_sga_satay_pd=pd.DataFrame(diff_sga_satay).T

diff_sga_satay_pd.head()

# +
from locale import normalize
from matplotlib import pyplot as plt

# remove Nan and inf from dataset

diff_sga_satay_pd=diff_sga_satay_pd[diff_sga_satay_pd.loc[:,"difference"]!=np.inf]
diff_sga_satay_pd=diff_sga_satay_pd[diff_sga_satay_pd.loc[:,"difference"]!=np.nan]
diff_sga_satay_pd=diff_sga_satay_pd[diff_sga_satay_pd.loc[:,"difference"]!=-np.inf]

fig,ax=plt.subplots(1,2,figsize=(10,5))
ax[0].hist(diff_sga_satay_pd.loc[:,"difference"],bins=100,cumulative=True,alpha=0.6);
ax[0].set_xlabel("Difference in fitness between SGA and SATAY per gene")
ax[0].set_ylabel("Number of genes")
ax[0].set_title("Cumulative distribution of fitness differences")
ax[1].scatter(diff_sga_satay_pd.loc[:,"sga"],diff_sga_satay_pd.loc[:,"satay"],s=5,alpha=0.8);
ax[1].set_xlabel("SGA fitness")
ax[1].set_ylabel("SATAY fitness")
ax[1].plot([0,1.4],[0,1.4],color="red")
ax[1].set_xlim(0.2,1.4)
ax[1].set_ylim(0.2,1.4)
ax[1].set_title("How much are the fitness values different?")
# -



# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
# ---

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys
from collections import defaultdict

data_all=pd.read_excel('../postprocessed-data/postprocessed_data_all_backgrounds.xlsx',engine="openpyxl")

# +
## Refine data 
data_all.index=np.arange(0,len(data_all))
data_all.drop(columns=['Unnamed: 1'],inplace=True)
data_all.fillna(0,inplace=True)
data_all.rename(columns={'Unnamed: 0':'background'},inplace=True)

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
indexes_back=[] # indexes for the start of each background

for i in keys:
    indexes_back.append(np.where(data_all.loc[:,"background"]==i)[0])

# Filling the name of the background according to the index
for k in np.arange(0,len(indexes_back)-1):
    
    data_all.loc[np.arange(indexes_back[k][0],indexes_back[k+1][0]),"background"]=keys[k]

data_all.loc[np.arange(indexes_back[-1][0],len(data_all)),"background"]=keys[-1] # for the last key
# -

# Compute average transposon density per background 
data_all_processed=data_all.copy()
for key in keys:
    sum_tr=data_all_processed[data_all_processed.loc[:,"background"]==key]["Ninsertions"].sum()
    genome=data_all_processed[data_all_processed.loc[:,"background"]==key]["Nbasepairs"].sum()
    data_all_processed.loc[data_all_processed.loc[:,"background"]==key,"transposon_density"]=sum_tr/genome

data_all_processed.index=data_all_processed.loc[:,"Standard_name"]
data_subset_processed=data_all_processed.copy()
data_subset_processed=data_subset_processed[data_subset_processed.loc[:,"Feature_type"]!=1]# remove unknown feature type


# +
# Convert all positions to float
from positions_float_dataframe import positions_float_dataframe
positions_float_pd=positions_float_dataframe(data_subset_processed,keys=keys)

positions_float_pd.to_csv("../postprocessed-data/genetic_positions_float_pd_all_backgrounds.csv",index=True)        
        

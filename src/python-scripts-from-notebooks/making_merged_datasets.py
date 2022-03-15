# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
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

# +
import os, sys
import warnings
import timeit
import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd

from collections import defaultdict

# +
data_a=pd.read_excel('../postprocessed-data/dbem3_a_pergene_insertions.xlsx',engine='openpyxl')
data_a=data_a.drop(['Unnamed: 0'],axis=1)

data_b=pd.read_excel('../postprocessed-data/dbem3_b_pergene_insertions.xlsx',engine='openpyxl')
data_b=data_b.drop(['Unnamed: 0'],axis=1)
# data.fillna(1,inplace=True)
data_merged=pd.read_excel('../postprocessed-data/dbem3_merged_pergene_insertions.xlsx',engine='openpyxl')
data_merged=data_merged.drop(['Unnamed: 0'],axis=1)
# -

data_a.head(2)

data_b.head(2)

data_merged.head(2)

from from_excel_to_list import from_excel_to_list

#functional version of recursive flatten which handles both tuples and lists, and lets you throw in any mix of positional arguments. 
# Returns a generator which produces the entire sequence in order,
flatten = lambda *n: (e for a in n
    for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))


def from_replicates_to_merged(replicates):
    """
    Takes a list of replicates and returns a single merged dataset.
  
    """
    #functional version of recursive flatten which handles both tuples and lists, and lets you throw in any mix of positional arguments. 
    # Returns a generator which produces the entire sequence in order,
    flatten = lambda *n: (e for a in n
        for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))
        
    # First, flatten the replicates
    replicates = list(flatten(replicates))

    data_merged_calculation=replicates[0].copy() # make a copy of the first replicate 
    data_merged_calculation.loc[:,"merged_reads"]=np.zeros(len(data_merged_calculation))
    data_merged_calculation.loc[:,"merged_insertions"]=np.zeros(len(data_merged_calculation))
    a=[]
    a_reads=[]
    for i in np.arange(0,len(replicates[0])):
        
        a.append(from_excel_to_list(replicates[0]['Insertion locations'][i]))
        b=from_excel_to_list(replicates[1]['Insertion locations'][i])
        a_reads.append(from_excel_to_list(replicates[0]['Reads per insertion location'][i]))
        b_reads=from_excel_to_list(replicates[1]['Reads per insertion location'][i])

        union=list(flatten(a,b)) # concatenate all the insertions per location for each replicate
        unique_union,unique_index,unique_counts=np.unique(union,return_index=True,return_counts=True)

        union_reads=list(flatten(a_reads,b_reads)) # concatenate all the reads per insertion location of all replicates
            
        union_dict=dict(zip(unique_union,unique_counts)) # make a dictionary of the union

        union_dict_pd=pd.DataFrame.from_dict(union_dict,orient='index',columns=['count']) # make a dataframe of the union

        union_dict_pd['insertions']=union_dict_pd.index
        union_dict_pd.index=np.arange(0,len(union_dict_pd))

        union_dict_pd['reads']=np.zeros(len(union_dict_pd))
        union_dict_pd['merged_reads']=np.zeros(len(union_dict_pd))

        for j in np.arange(0,len(union_dict_pd)):
            # if type(union_reads)!=int:
            union_dict_pd.loc[:,'reads'][j]=union_reads[j]
            union_dict_pd.loc[:,'merged_reads'][j]=union_dict_pd.loc[:,'count'][j]*union_dict_pd.loc[:,'reads'][j]

        total_reads=np.sum(union_dict_pd['merged_reads'])-np.max(union_dict_pd['merged_reads'])

        data_merged_calculation.loc[:,'Insertion locations'][i]=union_dict_pd['insertions'].values
        data_merged_calculation.loc[:,'merged_reads'][i]=total_reads
        data_merged_calculation.loc[:,'Reads per insertion location'][i]=union_dict_pd['reads'].values
        data_merged_calculation.loc[:,'merged_insertions'][i]=len(union_dict_pd['insertions'].values)

    return data_merged_calculation


data_merged=from_replicates_to_merged([data_a,data_b])

# +
data_merged_calculation=data_merged.copy()
data_merged_calculation.loc[:,"merged_reads"]=np.zeros(len(data_merged_calculation))
data_merged_calculation.loc[:,"merged_insertions"]=np.zeros(len(data_merged_calculation))

for i in np.arange(0,len(data_a)):
    a=from_excel_to_list(data_a['Insertion locations'][i])
    b=from_excel_to_list(data_b['Insertion locations'][i])
    a_reads=from_excel_to_list(data_a['Reads per insertion location'][i])
    b_reads=from_excel_to_list(data_b['Reads per insertion location'][i])

    if type(a)!=np.int64 and type(b)!=np.int64 and type(a)!=int and type(b)!=int:
        union=a+b
        unique_union,unique_index,unique_counts=np.unique(union,return_index=True,return_counts=True)
        
    else :
        union=list(flatten(a,b))
        unique_union,unique_index,unique_counts=np.unique(union,return_index=True,return_counts=True)

    if type(a_reads)!=np.int64 and type(b_reads)!=np.int64 and type(a_reads)!=int and type(b_reads)!=int:
        union_reads=a_reads+b_reads
    
    else :
        union_reads=list(flatten(a_reads,b_reads))
        
    union_dict=dict(zip(unique_union,unique_counts))

    union_dict_pd=pd.DataFrame.from_dict(union_dict,orient='index',columns=['count'])

    union_dict_pd['insertions']=union_dict_pd.index
    union_dict_pd.index=np.arange(0,len(union_dict_pd))

    union_dict_pd['reads']=np.zeros(len(union_dict_pd))
    union_dict_pd['merged_reads']=np.zeros(len(union_dict_pd))

    for j in np.arange(0,len(union_dict_pd)):
        if type(union_reads)!=int:
            union_dict_pd.loc[:,'reads'][j]=union_reads[j]
            union_dict_pd.loc[:,'merged_reads'][j]=union_dict_pd.loc[:,'count'][j]*union_dict_pd.loc[:,'reads'][j]

    total_reads=np.sum(union_dict_pd['merged_reads'])-np.max(union_dict_pd['merged_reads'])

    data_merged_calculation.loc[:,'Insertion locations'][i]=union_dict_pd['insertions'].values
    data_merged_calculation.loc[:,'merged_reads'][i]=total_reads
    data_merged_calculation.loc[:,'Reads per insertion location'][i]=union_dict_pd['reads'].values
    data_merged_calculation.loc[:,'merged_insertions'][i]=len(union_dict_pd['insertions'].values)

# -

data_merged_calculation

a=from_excel_to_list(data_a['Insertion locations'][0])
b=from_excel_to_list(data_b['Insertion locations'][0])
a_reads=from_excel_to_list(data_a['Reads per insertion location'][0])
b_reads=from_excel_to_list(data_b['Reads per insertion location'][0])

# +
union=a+b
unique_union,unique_index,unique_counts=np.unique(union,return_index=True,return_counts=True)

union_reads=a_reads+b_reads

union_dict=dict(zip(unique_union,unique_counts))

union_dict_pd=pd.DataFrame.from_dict(union_dict,orient='index',columns=['count'])

union_dict_pd['insertions']=union_dict_pd.index
union_dict_pd.index=np.arange(0,len(union_dict_pd))

union_dict_pd['reads']=np.zeros(len(union_dict_pd))
union_dict_pd['merged_reads']=np.zeros(len(union_dict_pd))

for i in range(len(union_dict_pd)):
    union_dict_pd.loc[:,'reads'][i]=union_reads[i]
    union_dict_pd.loc[:,'merged_reads'][i]=union_dict_pd.loc[:,'count'][i]*union_dict_pd.loc[:,'reads'][i]

np.sum(union_dict_pd['merged_reads'])-np.max(union_dict_pd['merged_reads'])

# +
set_a=set(a)
set_b=set(b)
s = set_a.union(set_b)  # union of both sets

# Difference in two sets
diff_element_a = set_b - set_a
diff_element_b = set_a - set_b
  
# union of difference + first list
output = a + list(diff_element_a) + list(diff_element_b)


# -

d = defaultdict(int)
for k in output:
    d[k] += 1
    

merged_insertions=pd.DataFrame(d.items(), columns=['Insertion locations', 'Count'])
merged_insertions.sort_values(by='Count', inplace=True,ascending=False)

# +
a_reads_dict=defaultdict(int)
j=0
for k in a:
    reads=a_reads[j]
    a_reads_dict[k] = reads
    j=j+1


# -

def get_reads_with_insertions(replicate_reads_per_gene,replicate_insertions_per_gene):
    """Get the reads per insertion location for a given replicate

    Parameters
    ----------
    replicate_reads_per_gene : list
        output from from_excel_to_list function
    """
    reads_dict=defaultdict(int)
    j=0
    for k in replicate_insertions_per_gene:
        reads=replicate_reads_per_gene[j]
        reads_dict[k] = reads
        j=j+1
    return reads_dict


a_reads_dict=get_reads_with_insertions(a_reads,a)
b_reads_dict=get_reads_with_insertions(b_reads,b)

a_reads_df=pd.DataFrame(a_reads_dict.items(), columns=['Insertion locations', 'Reads per insertion location'])
b_reads_df=pd.DataFrame(b_reads_dict.items(), columns=['Insertion locations', 'Reads per insertion location'])

a_reads_merged=pd.merge(pd.merge(merged_insertions,a_reads_df,on='Insertion locations'),b_reads_df,on='Insertion locations',how='left')
a_reads_merged.fillna(0,inplace=True)

print(np.sum(b_reads_df['Reads per insertion location'])-np.max(b_reads_df["Reads per insertion location"]))
print(np.sum(a_reads_df['Reads per insertion location'])-np.max(a_reads_df["Reads per insertion location"]))

print(np.sum(a_reads_merged['Reads per insertion location_x']))
print(np.sum(a_reads_merged['Reads per insertion location_y']))

group=a_reads_merged.groupby("Count").sum()
group

total=[]
for i in group.index:
    print(i)
    total.append(group.loc[i]['Reads per insertion location_x']*i+group.loc[i]['Reads per insertion location_y']*i)

np.sum(total)



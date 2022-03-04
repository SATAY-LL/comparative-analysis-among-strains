# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: 'Python 3.8.10 64-bit (''satay-dev'': conda)'
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
import seaborn as sns

from functions_comparison_replicates import getting_pergene_data
from functions_comparison_replicates import read_pergene_insertions_file
from functions_comparison_replicates import scatter_replicates

from functions_comparison_replicates import array2frame

from functions_comparison_replicates import create_pergene_insertions_excel_file


# +
folder="../data/"
wt_a_datafile = folder+"wt_a/WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_pergene.txt"
wt_b_datafile = folder+"wt_b/WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_pergene.txt"
bem3_a_datafile=folder+ "dbem3_a/bem3_a_pergene_tab.txt"
bem3_b_datafile=folder+"dbem3_b/bem3_b_pergene_tab.txt"
bem1bem3_a_datafile=folder+"dbem1dbem3_a/dbem1dbem3_a_tab_pergene.txt"
bem1bem3_b_datafile=folder+"dbem1dbem3_b/dbem1dbem3_b_tab_pergene.txt"
bem1_a_datafile=folder+"bem1-aid_a/bem1-aid_a_tab_pergene.txt"
bem1_b_datafile=folder+"bem1-aid_b/bem1-aid_b_tab_pergene.txt"
bem1_aid_bem3_a_datafile=folder+"bem1-aid-dbem3_a/bem1-aid-dbem3_a_tab_pergene.txt"
bem1_aid_bem3_b_datafile=folder+"bem1-aid-dbem3_b/bem1-aid-dbem3_b_tab_pergene.txt"
nrp1_a_datafile=folder+"dnrp1_a/nrp1_a_tab_pergene.txt"


wt_gene,wt_a_tn,wt_a_reads = getting_pergene_data(wt_a_datafile)
wt_gene,wt_b_tn,wt_b_reads = getting_pergene_data(wt_b_datafile)


bem3_gene,bem3_a_tn,bem3_a_reads=getting_pergene_data(bem3_a_datafile)
bem3_gene,bem3_b_tn,bem3_b_reads=getting_pergene_data(bem3_b_datafile)

bem1bem3_gene,bem1bem3_a_tn,bem1bem3_a_reads=getting_pergene_data(bem1bem3_a_datafile)
bem1bem3_gene,bem1bem3_b_tn,bem1bem3_b_reads=getting_pergene_data(bem1bem3_b_datafile)

bem1_gene,bem1_b_tn,bem1_b_reads=getting_pergene_data(bem1_b_datafile)
bem1_gene,bem1_a_tn,bem1_a_reads=getting_pergene_data(bem1_a_datafile)

bem1_aid_bem3_gene,bem1_aid_bem3_a_tn,bem1_aid_bem3_a_reads=getting_pergene_data(bem1_aid_bem3_a_datafile)
bem1_aid_bem3_gene,bem1_aid_bem3_b_tn,bem1_aid_bem3_b_reads=getting_pergene_data(bem1_aid_bem3_b_datafile)

# +

scatter_replicates([wt_a_tn,wt_b_tn,wt_a_reads,wt_b_reads],["WT_a","WT_b"],1000,10000)
scatter_replicates([bem3_a_tn,bem3_b_tn,bem3_a_reads,bem3_b_reads],["dbem3_a","dbem3_b"],1000,10000,save=True)
scatter_replicates([bem1_a_tn,bem1_b_tn,bem1_a_reads,bem1_b_reads],["dbem1_a","dbem1_b"],1000,10000,save=True)
scatter_replicates([bem1bem3_a_tn,bem1bem3_b_tn,bem1bem3_a_reads,bem1bem3_b_reads],["dbem1dbem3_a","dbem1dbem3_b"],1000,10000,save=True)
scatter_replicates([bem1_aid_bem3_a_tn,bem1_aid_bem3_b_tn,bem1_aid_bem3_a_reads,bem1_aid_bem3_b_reads],["bem1-aid-dbem3_a","bem1-aid-dbem3_b"],1000,10000,save=True)

# +
bem3_a=[bem3_a_tn,bem3_a_reads]
bem3_b=[bem3_b_tn,bem3_b_reads]

bem1bem3_a=[bem1bem3_a_tn,bem1bem3_a_reads]
bem1bem3_b=[bem1bem3_b_tn,bem1bem3_b_reads]

bem1_b=[bem1_b_tn,bem1_b_reads]
index=bem3_gene

df_bem3_a=array2frame(bem3_a,bem3_gene)
df_bem3_b=array2frame(bem3_b,bem3_gene)

df_bem1bem3_a=array2frame(bem1bem3_a,bem3_gene)
df_bem1bem3_b=array2frame(bem1bem3_b,bem3_gene)

df_bem1_b=array2frame(bem1_b,bem3_gene)


# -

df_bem3=pd.concat([df_bem3_a,df_bem3_b],keys=["a","b"])
df_bem1bem3=pd.concat([df_bem1bem3_a,df_bem1bem3_b],keys=["a","b"])

plt.hist(df_bem3.loc["a","reads"],bins=4000,alpha=0.5,label="dbem3_a")
plt.hist(df_bem3.loc["b","reads"],bins=4000,alpha=0.3,label="dbem3_b"); 
plt.xlim(0,6000)
plt.xlabel("Reads")
plt.title("dbem3")
plt.legend()

plt.hist(df_bem1bem3.loc["a","reads"],bins=4000,alpha=0.5,label="dbem1dbem3_a")
plt.hist(df_bem1bem3.loc["b","reads"],bins=4000,alpha=0.3,label="dbem1dbem3_b"); 
plt.xlim(0,6000)
plt.xlabel("Reads")
plt.title("dbem1dbem3")
plt.legend()

# +

for j in df_bem3.index:
    
    if df_bem3.loc[j,"transposons"]>5:

        df_bem3.loc[j,"reads-per-tr"]=df_bem3.loc[j,"reads"]/(df_bem3.loc[j,"transposons"]-1)

    else :
        df_bem3.loc[j,"reads-per-tr"]=0
# -

for j in df_bem1bem3.index:
    
    if df_bem1bem3.loc[j,"transposons"]>5:

        df_bem1bem3.loc[j,"reads-per-tr"]=df_bem1bem3.loc[j,"reads"]/(df_bem1bem3.loc[j,"transposons"]-1)

    else :
        df_bem1bem3.loc[j,"reads-per-tr"]=0

plt.hist(df_bem3.loc["a","reads-per-tr"],bins=200,alpha=0.5,label="dbem3_a")
plt.hist(df_bem3.loc["b","reads-per-tr"],bins=200,alpha=0.3,label="dbem3_b"); 
plt.xlim(0,60)
plt.xlabel("Reads per transposons per gene")
plt.title("dbem3")
plt.legend()

plt.hist(df_bem1bem3.loc["a","reads-per-tr"],bins=200,alpha=0.5,label="dbem3_a")
plt.hist(df_bem1bem3.loc["b","reads-per-tr"],bins=200,alpha=0.3,label="dbem3_b"); 
plt.xlim(0,60)
plt.xlabel("Reads per transposons per gene")
plt.title("dbem1dbem3")
plt.legend()

len(df_bem1bem3[(df_bem1bem3.loc[:, "transposons"]<5) & (df_bem1bem3.loc[:, "reads"]<10)])

len(df_bem3[(df_bem3.loc[:, "transposons"]<5)& (df_bem3.loc[:, "reads"]<10)])

len(df_bem1_b[(df_bem1_b.loc[:, "transposons"]<5)&(df_bem1_b.loc[:, "reads"]<10)])

# +
folder="../data/"
wt_a_datafile = folder+"wt_a/WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_pergene_insertions.txt"
wt_b_datafile = folder+"wt_b/WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_pergene_insertions.txt"
bem3_a_datafile=folder+ "dbem3_a/all_cleaned_fw_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem3_b_datafile=folder+"dbem3_b/yLIC137_8_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem1bem3_a_datafile=folder+"dbem1dbem3_a/yTW001_4_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem1bem3_b_datafile=folder+"dbem1dbem3_b/yTW001_6_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem1_a_datafile=folder+"bem1-aid_a/yWT03a_16_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem1_b_datafile=folder+"bem1-aid_b/yWT03a_21_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem1_aid_bem3_a_datafile=folder+"bem1-aid-dbem3_a/yWT04a_14_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem1_aid_bem3_b_datafile=folder+"bem1-aid-dbem3_b/yWT04a_23_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
nrp1_a_datafile=folder+"dnrp1_a/dnrp1-1_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene_insertions.txt"
nrp1_b_datafile=folder+"dnrp1_b/dnrp1-2_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene_insertions.txt"
bem3_merged=folder+"dbem3_merged/merged_ylic137_trimmed.sorted.bam_pergene_insertions.txt"
wt_merged=folder+"wt_merged/WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene_insertions.txt"
dnrp1_merged=folder+"dnrp1_merged/dnrp1_merged_dnrp1-1_dnrp1-2_trimmed.sorted.bam_pergene_insertions.txt"

bem3_a_trimmed=folder+"dbem3_a_trimmed_restriction_sites/ylic137_7_trimmed_out_restriction_sites_all_cleaned_fw_reads_trimmed.sorted.bam_pergene_insertions.txt"

datafile=[wt_a_datafile,wt_b_datafile,bem3_a_datafile,bem3_b_datafile,bem1bem3_a_datafile,
bem1bem3_b_datafile,bem1_a_datafile,bem1_b_datafile,bem1_aid_bem3_a_datafile,bem1_aid_bem3_b_datafile,bem1_aid_bem3_a_datafile,
nrp1_a_datafile,nrp1_b_datafile,bem3_merged,wt_merged,dnrp1_merged,bem3_a_trimmed]
# -

datafile=[bem3_a_datafile,bem3_a_trimmed]

result=create_pergene_insertions_excel_file(datafile)

# +
benoit_wt_a=folder+"WT_1_Benoit/ERR1533147_trimmed.sorted.bam_pergene_insertions.txt"
benoit_wt_b=folder+"WT_2_Benoit/ERR1533148_trimmed.sorted.bam_pergene_insertions.txt"

datafile=[benoit_wt_a,benoit_wt_b]

create_pergene_insertions_excel_file(datafile)

# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
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
wt_a_datafile = folder+"wt_a/WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_pergene_insertions.txt"
wt_b_datafile = folder+"wt_b/WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_pergene_insertions.txt"
bem3_a_datafile=folder+ "dbem3_a/all_cleaned_fw_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem3_b_datafile=folder+"dbem3_b/yLIC137_8_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem1bem3_a_datafile=folder+"dbem1dbem3_a/yTW001_4_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem1bem3_b_datafile=folder+"dbem1dbem3_b/yTW001_6_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem1_a_datafile=folder+"bem1-aid_a/yWT03a_16_trimmed_out_restriction_sites_yWT03a_16_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem1_b_datafile=folder+"bem1-aid_b/yWT0321_a_trimmed_out_restriction_sites_yWT0321_a_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem1_aid_bem3_a_datafile=folder+"bem1-aid-dbem3_a/yWT04a_14_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
bem1_aid_bem3_b_datafile=folder+"bem1-aid-dbem3_b/yWT04a_23_merged_cleaned_forward_reads_trimmed.sorted.bam_pergene_insertions.txt"
nrp1_a_datafile=folder+"dnrp1_a/dnrp1-1_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene_insertions.txt"
nrp1_b_datafile=folder+"dnrp1_b/dnrp1-2_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene_insertions.txt"
bem3_merged=folder+"dbem3_merged/merged_ylic137_trimmed.sorted.bam_pergene_insertions.txt"
wt_merged=folder+"wt_merged/WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene_insertions.txt"
dnrp1_merged=folder+"dnrp1_merged/dnrp1_merged_dnrp1-1_dnrp1-2_trimmed.sorted.bam_pergene_insertions.txt"
bem1_aid_merged=folder+"bem1-aid_merged/merged_datasets_trimmed_out_restriction_sites_merged_yTW003_new_trimmed.sorted.bam_pergene_insertions.txt"

datafile=[wt_a_datafile,wt_b_datafile,bem3_a_datafile,bem3_b_datafile,bem1bem3_a_datafile,
bem1bem3_b_datafile,bem1_a_datafile,bem1_b_datafile,bem1_aid_bem3_a_datafile,bem1_aid_bem3_b_datafile,bem1_aid_bem3_a_datafile,
nrp1_a_datafile,nrp1_b_datafile,bem3_merged,wt_merged,dnrp1_merged,bem1_aid_merged]
# -

datafile=[bem1_aid_merged]

result=create_pergene_insertions_excel_file(datafile)

# +
benoit_wt_a=folder+"WT_1_Benoit/ERR1533147_trimmed.sorted.bam_pergene_insertions.txt"
benoit_wt_b=folder+"WT_2_Benoit/ERR1533148_trimmed.sorted.bam_pergene_insertions.txt"

datafile=[benoit_wt_a,benoit_wt_b]

create_pergene_insertions_excel_file(datafile)

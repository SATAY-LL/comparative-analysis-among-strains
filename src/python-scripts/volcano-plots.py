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

# # Transposonmapper output data postprocessing 

## Importing the required python libraries 
import os, sys
import warnings
import timeit
import numpy as np
import pandas as pd 
import pkg_resources
import matplotlib.pyplot as plt
from transposonmapper.statistics import volcano

# # Volcano plots
#
# Do you want to compare two differente libraries to discover which genes stood out from their comparison? 
#
# Then do volcano plots!!
#

# ## Getting the volcano plot
#
# Look at the help of this function , [HERE](https://github.com/SATAY-LL/Transposonmapper/blob/main/transposonmapper/statistics/volcanoplot.py)

# +


path_a = r"../data/"
filelist_a = ["wt_a/WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_pergene.txt",
"wt_b/WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_pergene.txt"]
path_b = r"../data/"
filelist_b = ["dnrp1_a/dnrp1-1_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt",
"dnrp1_b/dnrp1-2_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt"]


variable = 'tn_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.01 #set threshold above which p-values are regarded significant
normalize=True

trackgene_list = ['nrp1','bem3','bem1','bem2'] # ["cdc42"]


figure_title = "WT vs dnrp1"

volcano_df_nrp1_wt = volcano(path_a=path_a, filelist_a=filelist_a,
            path_b=path_b, filelist_b=filelist_b,
            variable=variable,
            significance_threshold=significance_threshold,
            normalize=normalize,
            trackgene_list=trackgene_list,
            figure_title=figure_title)

# +
path_a = r"../data/"
filelist_a = ["wt_a/WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_pergene.txt",
"wt_b/WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_pergene.txt"]

path_b = r"../data/"

filelist_b = ["dbem3_a/bem3_a_pergene_tab.txt",
"dbem3_b/bem3_b_pergene_tab.txt"]

variable = 'tn_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.01 #set threshold above which p-values are regarded significant
normalize=True

trackgene_list = ['bem3','bem1','nrp1','bem2'] # ["cdc42"]


figure_title = "WT vs dbem3 "

volcano_df_bem3_wt = volcano(path_a=path_a, filelist_a=filelist_a,
            path_b=path_b, filelist_b=filelist_b,
            variable=variable,
            significance_threshold=significance_threshold,
            normalize=normalize,
            trackgene_list=trackgene_list,
            figure_title=figure_title)

# +
path_a = r"../data/"
filelist_a = ["dnrp1_a/dnrp1-1_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt",
"dnrp1_b/dnrp1-2_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt"]

path_b = r"../data/"

filelist_b = ["dbem3_a/bem3_a_pergene_tab.txt",
"dbem3_b/bem3_b_pergene_tab.txt"]

variable = 'tn_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.01 #set threshold above which p-values are regarded significant
normalize=True

trackgene_list = ['bem3','bem1','nrp1','bem2'] # ["cdc42"]


figure_title = "dnrp1 vs dbem3"

volcano_df_bem3_nrp1 = volcano(path_a=path_a, filelist_a=filelist_a,
            path_b=path_b, filelist_b=filelist_b,
            variable=variable,
            significance_threshold=significance_threshold,
            normalize=normalize,
            trackgene_list=trackgene_list,
            figure_title=figure_title)

# +
path_a = r"../data/"
filelist_a = ["dbem1dbem3_a/dbem1dbem3_a_tab_pergene.txt",
"dbem1dbem3_b/dbem1dbem3_b_tab_pergene.txt"]

path_b = r"../data/"

filelist_b = ["dbem3_a/bem3_a_pergene_tab.txt",
"dbem3_b/bem3_b_pergene_tab.txt"]

variable = 'tn_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.01 #set threshold above which p-values are regarded significant
normalize=True

trackgene_list = ['bem3','bem1','nrp1','bem2'] # ["cdc42"]


figure_title = "dbem1dbem3 vs dbem3"

volcano_df = volcano(path_a=path_a, filelist_a=filelist_a,
            path_b=path_b, filelist_b=filelist_b,
            variable=variable,
            significance_threshold=significance_threshold,
            normalize=normalize,
            trackgene_list=trackgene_list,
            figure_title=figure_title)

# +
path_a = r"../data/"
filelist_a = ["wt_a/WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_pergene.txt",
"wt_b/WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_pergene.txt"]

path_b = r"../data/"

filelist_b = ["bem1-aid_a/bem1-aid_a_tab_pergene.txt",
"bem1-aid_b/bem1-aid_b_tab_pergene.txt"]

variable = 'tn_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.01 #set threshold above which p-values are regarded significant
normalize=True

trackgene_list = ['bem3','bem1','nrp1','bem2'] # ["cdc42"]


figure_title = "WT vs bem1-aid"

volcano_df = volcano(path_a=path_a, filelist_a=filelist_a,
            path_b=path_b, filelist_b=filelist_b,
            variable=variable,
            significance_threshold=significance_threshold,
            normalize=normalize,
            trackgene_list=trackgene_list,
            figure_title=figure_title)
# -

from annotate_volcano import annotate_volcano   #import annotate_volcano function
annotate_volcano(volcano_df,[2.1,-1.5],[2,2])

# +
path_a = r"../data/"
filelist_b = ["bem1-aid_a/bem1-aid_a_tab_pergene.txt",
"bem1-aid_b/bem1-aid_b_tab_pergene.txt"]

path_b = r"../data/"

filelist_b = ["bem1-aid-dbem3_a/bem1-aid-dbem3_a_tab_pergene.txt",
"bem1-aid-dbem3_b/bem1-aid-dbem3_b_tab_pergene.txt"]

variable = 'tn_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.01 #set threshold above which p-values are regarded significant
normalize=True

trackgene_list = ['bem3','bem1','nrp1','bem2'] # ["cdc42"]


figure_title = "bem1-aid vs bem1-aid-dbem3"

volcano_df = volcano(path_a=path_a, filelist_a=filelist_a,
            path_b=path_b, filelist_b=filelist_b,
            variable=variable,
            significance_threshold=significance_threshold,
            normalize=normalize,
            trackgene_list=trackgene_list,
            figure_title=figure_title)

# +
path_a = r"../data/"
filelist_a = ["dbem1dbem3_a/dbem1dbem3_a_tab_pergene.txt",
"dbem1dbem3_b/dbem1dbem3_b_tab_pergene.txt"]

path_b = r"../data/"

filelist_b = ["bem1-aid-dbem3_a/bem1-aid-dbem3_a_tab_pergene.txt",
"bem1-aid-dbem3_b/bem1-aid-dbem3_b_tab_pergene.txt"]

variable = 'tn_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.01 #set threshold above which p-values are regarded significant
normalize=True

trackgene_list = ['bem3','bem1','nrp1','bem2'] # ["cdc42"]


figure_title = "bem1dbem3d vs bem1-aid-dbem3"

volcano_df = volcano(path_a=path_a, filelist_a=filelist_a,
            path_b=path_b, filelist_b=filelist_b,
            variable=variable,
            significance_threshold=significance_threshold,
            normalize=normalize,
            trackgene_list=trackgene_list,
            figure_title=figure_title)

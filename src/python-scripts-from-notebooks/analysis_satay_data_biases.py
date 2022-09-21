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

import wiggelen as wg

# +
## Importing the required python libraries and using the transposon_bias2centromeres function

wigfile="../data/wt_merged/WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_clean.wig"
centromeres_file="../data/centromeres_genomic_locations.csv"

wigfile_mutant="../data/dbem3_merged/merged_ylic137_trimmed.sorted.bam_clean.wig"
from functions_satay_biases import transposon_bias2centromeres

fig, distance2cent_all = transposon_bias2centromeres(wigfile,centromeres_file,save=True)

# +
### Independence of this behaviour in a different genetic background (for sumplementary material)

fig,distance2cent_all=transposon_bias2centromeres(wigfile_mutant,centromeres_file,save=True)

# +
## Bias 2: More transposon enrichment in non coding regions than coding regions

# import pergene file for the Coding regions

#file="../postprocessed-data/wt_merged_pergene_insertions.xlsx"
file="../postprocessed-data/dbem3_merged_pergene_insertions.xlsx"
pergene_data=pd.read_excel(file,engine='openpyxl',index_col="Unnamed: 0")

# import wig file for all the regions

wigfile_path="../data/wt_merged/WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_clean.wig"
wigfile_mutant="../data/dbem3_merged/merged_ylic137_trimmed.sorted.bam_clean.wig"


wig_reads=[]

for x in wg.walk(open(wigfile_mutant)):
    wig_reads.append(x)

wig_reads = np.array(wig_reads)


# +
pergene=defaultdict(dict)
for i in pergene_data.index:
    pergene[i]["chrom"]=pergene_data.loc[i,"Chromosome"]
    pergene[i]["start"]=pergene_data.loc[i,"Start location"]
    pergene[i]["end"]=pergene_data.loc[i,"End location"]
    pergene[i]["tr"]=pergene_data.loc[i,"Insertions"]

pergene_pd=pd.DataFrame.from_dict(pergene,orient="index")

# +
pos=[]
chrom=[]
for i in wig_reads:
    chrom.append(i[0]) ## chromosomes
    pos.append(float(i[1])) ## positions

chrom=np.array(chrom)
pos=np.array(pos)
# -

chrom_conversion={"chrI":"I","chrII":"II","chrIII":"III",
"chrIV":"IV","chrV":"V","chrVI":"VI","chrVII":"VII","chrVIII":"VIII",
"chrIX":"IX","chrX":"X","chrXI":"XI","chrXII":"XII","chrXIII":"XIII",
"chrXIV":"XIV","chrXV":"XV","chrXVI":"XVI"}

non_cod=[]
for i in chrom_conversion.items():
    tmp=pergene_pd[pergene_pd["chrom"]==i[1]]
    tmp_pos=pos[chrom==i[0]]
    non_cod_1=np.setdiff1d(tmp_pos,tmp["start"])
    non_cod_2=np.setdiff1d(tmp_pos,tmp["end"])
    non_cod.append(len(np.unique(np.concatenate((non_cod_1,non_cod_2)))))

total_tr_non_coding=np.sum(non_cod)
total_tr_non_coding

## intergenic regions larger than 100bp per chromosome from yeast mine 
# https://yeastmine.yeastgenome.org/yeastmine/template.do?name=Chromosome_IntergenicSequence&scope=all
non_coding_total=100+400+163+751+289+129+548+287+217+366+319+511+467+394+549+461
non_coding_total

# +
average_over_cdcs=pergene_pd["tr"].sum()/len(pergene_data)
average_over_non_coding=total_tr_non_coding/non_coding_total

std_over_cdcs=pergene_pd["tr"].std()/np.sqrt(len(pergene_data))
std_over_non_coding=np.std(non_cod)/np.sqrt(non_coding_total)



std_over_cdcs,std_over_non_coding

# +
figure,axes=plt.subplots(1,1,figsize=(5,5))
plotline1, caplines1, barlinecols1=axes.errorbar(["Coding"],average_over_cdcs,
yerr=std_over_cdcs,fmt="",color="black",lolims=True,capsize = 0, ls='None')

plotline2, caplines2, barlinecols2=axes.errorbar(["Non-coding"],average_over_non_coding,yerr=std_over_non_coding,
fmt="",color="black",lolims=True,capsize = 0, ls='None')

caplines1[0].set_marker('_')
caplines1[0].set_markersize(20)

caplines2[0].set_marker('_')
caplines2[0].set_markersize(20)

axes.bar(["Coding","Non-coding"],[average_over_cdcs,average_over_non_coding],color=["#1f77b4","#ff7f0e"])

axes.set_ylabel("Average of transposon insertions")
plt.tight_layout()
figure.savefig("../figures/fig_satay_bias_transposon_insertions_coding_noncoding_dbem3.png",dpi=300,transparent=True)
# -



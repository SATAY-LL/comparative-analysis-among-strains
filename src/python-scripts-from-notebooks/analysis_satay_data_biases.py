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


# +
from functions_satay_biases import transposon_bias2centromeres

fig, distance2cent_all = transposon_bias2centromeres(wigfile,centromeres_file,save=True)

# +
### Independence of this behaviour in a different genetic background (for sumplementary material)

fig,distance2cent_all=transposon_bias2centromeres(wigfile_mutant,centromeres_file,save=False)

# +
## Bias 2: More transposon enrichment in non coding regions than coding regions

# import pergene file for the Coding regions

file="../postprocessed-data/wt_merged_pergene_insertions.xlsx"
file_mutant="../postprocessed-data/dbem3_merged_pergene_insertions.xlsx"
pergene_data=pd.read_excel(file,engine='openpyxl',index_col="Unnamed: 0")

# import wig file for all the regions

wigfile_path="../data/wt_merged/WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_clean.wig"
wigfile_mutant="../data/dbem3_merged/merged_ylic137_trimmed.sorted.bam_clean.wig"


wig_reads=[]

for x in wg.walk(open(wigfile_path)):
    wig_reads.append(x)

wig_reads = np.array(wig_reads)


# +
pergene=defaultdict(dict)
for i in pergene_data.index:
    pergene[i]["chrom"]=pergene_data.loc[i,"Chromosome"]
    pergene[i]["start"]=pergene_data.loc[i,"Start location"]
    pergene[i]["end"]=pergene_data.loc[i,"End location"]
    pergene[i]["tr"]=pergene_data.loc[i,"Insertions"]
    pergene[i]["length"]=pergene_data.loc[i,"End location"]-pergene_data.loc[i,"Start location"]
    pergene[i]["tr_density"]=pergene_data.loc[i,"Insertions"]/pergene[i]["length"]

pergene_pd=pd.DataFrame.from_dict(pergene,orient="index")
# -

pergene_pd.head()

# +
## Distangling the wig reads array into two arrays: one with the chromosome and the other with the position
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

# Computing number of transposons in non coding regions
non_cod=[]
tr_non_coding_chrom={}
for i in chrom_conversion.items():
    tmp=pergene_pd[pergene_pd["chrom"]==i[1]]
    tmp_pos=pos[chrom==i[0]]
    non_cod_1=np.setdiff1d(tmp_pos,tmp["start"])
    non_cod_2=np.setdiff1d(tmp_pos,tmp["end"])
    non_cod.append(len(np.unique(np.concatenate((non_cod_1,non_cod_2)))))
    tr_non_coding_chrom[i[1]]=len(np.unique(np.concatenate((non_cod_1,non_cod_2))))

## intergenic regions larger than 100bp per chromosome from yeast mine 
# https://yeastmine.yeastgenome.org/yeastmine/template.do?name=Chromosome_IntergenicSequence&scope=all
non_coding_elements={"I":100,"II":400,"III":163,
"IV":751,"V":289,"VI":129,"VII":548,"VIII":287,
"IX":217,"X":366,"XI":319,"XII":511,"XIII":467,"XIV":394,
"XV":549,"XVI":461}


# +
# size of every chromosome in bp from https://www-ncbi-nlm-nih-gov.tudelft.idm.oclc.org/genome/?term=Saccharomyces%20cerevisiae%5BOrganism%5D&cmd=DetailsSearch
chromosome_size={"I":230218,"II":813184,"III":316620,
"IV":1531933,"V":576874,"VI":270161,"VII":1090940,
"VIII":562643,"IX":439888,"X":745751,"XI":666816,
"XII":1078177,"XIII":924431,"XIV":784333,"XV":1091291,
"XVI":948066}

# Size of coding regions per chromosome
size_coding_regions={}
for i in pergene_pd["chrom"].unique():
    size_coding_regions[i]=pergene_pd[pergene_pd["chrom"]==i]["length"].sum()
    
size_coding_regions

# +
# density of non coding regions per chromosome
non_coding_density={k: v / (chromosome_size[k]-size_coding_regions[k]) for k, v in tr_non_coding_chrom.items()}


# density of coding regions per chromosome
coding_density_per_chrom={}
for chrom in pergene_pd["chrom"].unique():
    coding_density_per_chrom[chrom]=pergene_pd[pergene_pd["chrom"]==chrom]["tr_density"].mean()


# +
average_density_non_coding=np.mean(list(non_coding_density.values()))
average_density_non_coding_std=np.std(list(non_coding_density.values()))


average_density_cds=np.mean(list(coding_density_per_chrom.values()))
average_density_cds_std=np.std(list(coding_density_per_chrom.values()))

average_density_cds,average_density_non_coding

# +
figure,axes=plt.subplots(1,1,figsize=(6,5))
plotline1, caplines1, barlinecols1=axes.errorbar(["Coding"],average_density_cds,
yerr=average_density_cds_std,fmt="",color="black",lolims=True,capsize = 0, ls='None')

plotline2, caplines2, barlinecols2=axes.errorbar(["Non-coding"],average_density_non_coding,
yerr=average_density_non_coding_std,fmt="",color="black",lolims=True,capsize = 0, ls='None')

caplines1[0].set_marker('_')
caplines1[0].set_markersize(20)

caplines2[0].set_marker('_')
caplines2[0].set_markersize(20)

axes.bar(["Coding","Non-coding"],[average_density_cds,average_density_non_coding],color=["#1f77b4","#ff7f0e"])

axes.set_ylabel("Average of transposon insertions per bp")
plt.tight_layout()
figure.savefig("../figures/fig_satay_bias_transposon_insertions_coding_noncoding.png",dpi=300,transparent=True)
# -



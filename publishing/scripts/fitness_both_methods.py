
"""
Created on Thu Sep 14 14:38:28 2023

@author: ekingma1/linigodelacruz
"""

# general basic python libraries
import numpy as np
import pandas as pd

# special home made helper python functions
from data_import_tools import per_gene, wiggle,fitness_average,fitness_domains,compute_fitness_diff



## Computation of the fitness average and fitness per domain for all genetic backgrounds used in the research

# +
## Importing pergene files 

pergene_files=["WT_merged_pergene_insertions.txt"]

wig_files=["WT_merged.wig"]


# +
# Import data

gene_inserts=[]
wig_inserts=[]
for i in range(len(pergene_files)):
    gene_inserts.append(per_gene(pergene_files[i]))
    wig_inserts.append(wiggle(wig_files[i]))

# +
data_domains=pd.read_excel("genomic-domains-wt.xlsx",index_col="Unnamed: 0")

## data from yeastmine
domains_names=pd.read_csv('Domains_all_genes_protein_coordinates_yeastmine.tsv',sep="\t")
domains_names.index=domains_names["Gene Name"]


if __name__=="__main__":
    

    ## Fitness average 
    fitness_all_pd=fitness_average(pergene_files=pergene_files,
                                gene_inserts=gene_inserts,
                                wig_inserts=wig_inserts)

    # ### Fitness for domains 

    fitness_domain_all_pd=fitness_domains(pergene_files=pergene_files,
                                            gene_inserts=gene_inserts,
                                            wig_inserts=wig_inserts,
                                            data_domains=data_domains)


    # Normalize each fitness value by the median fitness value of all genes 
        
    x=fitness_domain_all_pd.loc[:,"fitness_domains"].tolist()
    flattened_list = [item for sublist in x if isinstance(sublist, list) for item in sublist] + [item for item in x if isinstance(item, (int, float))]
    ref=np.nanmedian(flattened_list)
    # divide each fitness value by the reference value
    for j in fitness_domain_all_pd.index:
        fitness_domain_all_pd.loc[j,"fitness_domains"]=fitness_domain_all_pd.loc[j,"fitness_domains"]/ref

    # +
    ## Compute the differences in all backgrounds
    fitness_corrected_pd=compute_fitness_diff(fitness_average=fitness_all_pd,
                                            fitness_domains=fitness_domain_all_pd)





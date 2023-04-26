
import pandas as pd
import numpy as np

def filter_fitness(data_fitness,backgrounds,goi=["BEM1","BEM3","NRP1"],discard=["Not enough flanking regions"],set2zero=["Not enough reads","Not enough insertions"],cols=["fitness_gene","fitness_domains_corrected"]):

    data=[]
    for k in backgrounds:
        f=data_fitness.loc[k]
        if len(discard)==1:
            f=f[f.loc[:,cols[0]]!=discard[0]]
        else:
            for i in discard:
                f=f[f.loc[:,cols[0]]!=i]

        for i in f.index:
            if i!=goi[0] and i!=goi[1] and i!=goi[2]:
                if f.loc[i,cols[1]]==set2zero[0] or f.loc[i,cols[1]]==set2zero[1]:
                    f.drop(i,inplace=True)
                else:
                    pass
            else:
                if f.loc[i,cols[1]]==set2zero[0] or f.loc[i,cols[1]]==set2zero[1]:
                   f.loc[i,cols[1]]=0
                else:
                    pass
        
        f[f.loc[:,cols[0]]==set2zero[0]]=0
        f[f.loc[:,cols[0]]==set2zero[1]]=0
        
        data.append(f)

    return pd.concat(data,axis=0,keys=backgrounds)
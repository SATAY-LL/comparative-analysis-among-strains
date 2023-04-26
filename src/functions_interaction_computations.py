
import pandas as pd
import numpy as np

from collections import defaultdict
from scipy import stats

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



def digenic_GI(data,goi="BEM1",col_fitness="fitness_gene",backg=["wt_a","wt_b","bem1-aid_a","bem1-aid_b"],significance=0.1):
   
## Removing genes with zero fitness from all backgrounds for the interaction analysis

    
    data_backg=data.loc[backg]
    data_fitness2interact=data_backg[data_backg.loc[:,col_fitness]!=0]
     

    if goi=="BEM1":
        goi_f_a=(data_backg.loc[backg[0],col_fitness].loc[goi]+data_backg.loc[backg[0]].loc[goi,"fitness_domains_corrected"])/2
        goi_f_b=(data_backg.loc[backg[1],col_fitness].loc[goi]+data_backg.loc[backg[1]].loc[goi,"fitness_domains_corrected"])/2
    else:
        goi_f_a=data_backg.loc[backg[0],col_fitness].loc[goi]
        goi_f_b=data_backg.loc[backg[1],col_fitness].loc[goi]

    data_0=data_fitness2interact.loc[backg[0],col_fitness]
    data_1=data_fitness2interact.loc[backg[1],col_fitness]
    data_2=data_fitness2interact.loc[backg[2],col_fitness]*goi_f_a/(2**np.median(data_backg.loc[backg[2],col_fitness]))## scale to the fitness of the gene of interest
    data_3=data_fitness2interact.loc[backg[3],col_fitness]*goi_f_b/(2**np.median(data_backg.loc[backg[2],col_fitness]))## scale to the fitness of the gene of interest


   ## Computing the interactors of BEM1 using the whole fitness of the gene 
    
    intersection_genes=list((set(data_0.index)&set(data_1.index)&
    set(data_2.index)&set(data_3.index)))

    # index_num=[]
    # for i in intersection_genes:
    #     index_num.append(data_a.index.get_loc(i))


    significance_threshold = significance #set significance threshold
    gi=defaultdict(dict)
    ttest_tval_list = [np.nan]*2 #initialize list for storing t statistics
    ttest_pval_list = [np.nan]*2 #initialize list for storing p-values
    signif_thres_list = False #initialize boolean list for indicating datapoints with p-value above threshold
    fc_list = [np.nan]*2
    for gene in intersection_genes:
        geneX=gene
        
        #geneX_f_a=data_a.loc[geneX,"fitness_gene"]
        geneX_f_a=data_0[geneX]
        
        geneX_f_b=data_1[geneX]
            
        geneXbem1_f_a=data_2[geneX]
        geneXbem1_f_b=data_3[geneX]
                
        variable_a_array=[geneXbem1_f_a,geneXbem1_f_b]# in satay the double mutant fitness is only the fitness of the genex on dbem1 but we dont have the other way around data
        variable_b_array=[(geneX_f_a*goi_f_a),(geneX_f_b*goi_f_b)]
        
        ttest_val = stats.ttest_ind(variable_a_array, variable_b_array) #T-test
        gi[gene]["p_statistic"]=ttest_val[1]
        ttest_tval_list = ttest_val[0]
        gi[gene]["t_statistic"]=ttest_tval_list
        if not ttest_val[1] == 0: #prevent p=0 to be inputted in log
            ttest_pval_list = -1*np.log10(ttest_val[1])
            
        else:
            ttest_pval_list = 0
        gi[gene]["p_value"]=ttest_pval_list
        if ttest_pval_list >= -1*np.log10(significance_threshold):
            gi[gene]["significance"]=True
        else:
            gi[gene]["significance"]=False
        #DETERMINE FOLD CHANGE PER GENE
        
        fc_list=(np.mean(variable_a_array)-np.mean(variable_b_array))
        gi[gene]["fold_change"]=fc_list
        # gi[gene]["multiple-delete-fitness"]=np.mean(variable_a_array)
        # gi[gene]["single-delete-fitness-product"]=np.mean(variable_b_array)
        
        gi[gene]["e_a"]=variable_a_array[0]-variable_b_array[0]
        gi[gene]["e_b"]=variable_a_array[1]-variable_b_array[1]
        gi[gene]["fold_change_std"]= np.std([gi[gene]["e_a"],gi[gene]["e_b"]])
            
    return gi
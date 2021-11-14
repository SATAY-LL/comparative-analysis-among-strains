# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 11:31:37 2021

@author: linigodelacruz
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

def getting_r(datasets): 
    """
    This function computes the maximum rates per gene,effectively the fitness for 
    a strain with a knockout in gene X, given an intergenic model . The intergenic model assumes a 
    finite population growth until the carrying capacity (K). This population growth describes a 
    logistic growth. 
    
    Parameters
    ----------
    datasets : TYPE - list of dataframes
        These are dataframes containing the information on insertions per gene and reads per gene 
        given by the transposon data sequencing.
        The columns should have the following names:
            number_of_read_per_gene
            number_of_transposon_per_gene
            

    Returns
    -------
    r : array of length of the list given as input 
        array containing the maximum rates per gene per dataset in the population according an intergenic model
   

    """
    ## add a for loop over times and compute the rates over those times and averaged them out and std 
 
    T=90
    r=[]  
    r_non_filter=[]
    volume=3000 #mL
    C=1# aprox density of cells per ml at the start of the reseed 
    #C_end=C*8 # aprox density of cells per ml at the end of the reseed (counting the number of generations )
 
    #K=C_end
    # check if datasets is a llist of dataframes or just one 
    if type(datasets)=="list": 
        for i in np.arange(0,len(datasets)):
        
        
            
            
            K=np.sum(datasets[i]['reads-per-tr']) # carrying capacity
            N=datasets[i]['reads-per-tr'] # contribution of each mutant to the population density 
            # it will compute a carrying capacity per dataset 
        
                        
            r.append(np.log(N*K/(K-N))/T)
            
            # K_n=np.sum(datasets[i]['Nreadsperinsrt']) # it will compute a carrying capacity per dataset 
            # N_n=datasets[i]['Nreadsperinsrt']
        

                        
            # r_non_filter.append(np.log(N_n*K_n/(K_n-N_n))/T)
    else: 
            datasets.replace([np.inf, -np.inf], np.nan, inplace=True) # replace inf by the max number of the dataset
            datasets.fillna(max(datasets["reads-per-tr"]),inplace=True)
            K=np.sum(datasets['reads-per-tr']) # carrying capacity
            N=datasets['reads-per-tr'] # contribution of each mutant to the population density 
            # it will compute a carrying capacity per dataset 
        
                        
            r.append(np.log(N*K/(K-N))/T)
            
            # K_n=np.sum(datasets['Nreadsperinsrt']) # it will compute a carrying capacity per dataset 
            # N_n=datasets['Nreadsperinsrt']
        

                        
            # r_non_filter.append(np.log(N_n*K_n/(K_n-N_n))/T)

            datasets["rates-intergenic"]=r[0]
              
        
        
    return r

#### adding features for analysis . File required: pergene_insertions files 
def adding_features2dataframe(data):
    #data["N_reads"]=np.nan
    data["tr-density"]=np.nan
    #data["std-reads"]=np.nan
    #data["mean-reads"]=np.nan
    data["reads-per-tr"]=np.nan
    #data["tranposons"]=np.nan
    for j in data.index:
        #data["N_reads"][j]=np.sum(data["Reads per insertion location"][j])
        data["tr-density"][j]=data["Insertions"][j]/(data["End location"][j]-data["Start location"][j])
        #data["std-reads"][j]=np.std(data["Reads per insertion location"][j])
        #data["mean-reads"][j]=np.mean(data["Reads per insertion location"][j])
        if data["Insertions"][j]>5:
            data["reads-per-tr"][j]=data["Reads"][j]/((data["Insertions"][j])-1)
        else:
            data["reads-per-tr"][j]=0
        #data["tranposons"][j]=len(data["Reads per insertion location"][j])
    return data 


### Configuring the dataframes for analyses

def pd_to_compare_libraries(data_list,filesnames,datasets,norm=True):
    #pd_list=pd.DataFrame(data_list,index=np.arange(0,len(datasets_for_index)))
    pd_list=pd.DataFrame(data_list)
    pd_list=pd_list.T
       
    columns=[]
    
    for i in np.arange(0,len(filesnames)):
        if filesnames[i].startswith('d'): #dnrp1
            columns.append(filesnames[i][0:7]+filesnames[i][27:29])
        elif filesnames[i].startswith('W'):
            columns.append(filesnames[i][0:2]+filesnames[i][22:24])
            
        
    pd_list.columns=columns
    
    #pd_list.index=datasets_for_index[0]['gene_name']
    pd_list.index=datasets[0].loc[pd_list.index,'gene_name']
    pd_list.fillna(0,inplace=True)
    pd_list.replace(float('-inf'),0,inplace=True)
    
    pd_dnrp1=pd_list.filter(regex='nrp')
    
    pd_list['average_dnrp1']=np.mean(pd_dnrp1,axis=1)
    pd_list['std_dnrp1']=np.std(pd_dnrp1,axis=1)
    
    pd_list['merged_dnrp1']=np.sum(pd_dnrp1,axis=1)
    
    pd_wt=pd_list.filter(regex='WT')
    pd_list['average_wt']=np.mean(pd_wt,axis=1)
    pd_list['std_wt']=np.std(pd_wt,axis=1)
    
    pd_list['merged_wt']=np.sum(pd_wt,axis=1)
    
    
    pd_list['average_ratio_mutant_vs_wt']=pd_list['average_dnrp1']/pd_list['average_wt']
    if norm==True:
        
        # pd_list['norm2max_wt']=pd_list.loc[:,'average_wt']/pd_list.loc[:,'average_wt'].max()
        # pd_list['norm2max_dnrp1']=pd_list.loc[:,'average_dnrp1']/pd_list.loc[:,'average_dnrp1'].max()
        
        pd_list['norm2HO_wt']=pd_list.loc[:,'average_wt']/pd_list.loc['HO','average_wt']
        pd_list['norm2HOdnrp1']=pd_list.loc[:,'average_dnrp1']/pd_list.loc['HO','average_dnrp1']

        
    return pd_list
### function to filter the data 
def filtering_significant_genes(dataset,p_value_th,fc_th):
    dataset_highfc=dataset[dataset.loc[:,'fc_log2']>fc_th] 
    dataset_highp=dataset[dataset.loc[:,'-log_p']>p_value_th]
    dataset_significant_genes_wt= pd.merge(dataset_highfc, dataset_highp, how='inner', on=['-log_p'],left_index=True, right_index=True)
    dataset_significant_genes_wt.drop(columns=['fc_log2_y'],inplace=True)
    dataset_significant_genes_wt.rename(columns={'fc_log2_x': 'fc_log2'},inplace=True)
    
    
    dataset_highfc=dataset[dataset.loc[:,'fc_log2']<-fc_th] 
    dataset_highp=dataset[dataset.loc[:,'-log_p']>p_value_th] 
    dataset_significant_genes_mutant= pd.merge(dataset_highfc, dataset_highp, how='inner', on=['-log_p'],left_index=True, right_index=True)
    dataset_significant_genes_mutant.drop(columns=['fc_log2_y'],inplace=True)
    dataset_significant_genes_mutant.rename(columns={'fc_log2_x': 'fc_log2'},inplace=True)
    
    return dataset_significant_genes_wt,dataset_significant_genes_mutant


def viz_data(dataset,p_value_th,fc_th,filenames,saving=True,plotting=True):
    
    results=[]
    if type(dataset)==list:
        
        for i in np.arange(0,len(dataset)):
            
            ref,exp=filtering_significant_genes(dataset[i],p_value_th,fc_th)
            results.append([ref,exp])
    else : 
        
        ref,exp=filtering_significant_genes(dataset,p_value_th,fc_th)
        results.append([ref,exp])        
### Plotting all datasets 
    for data in np.arange(0,len(results)):
        fig = plt.figure(figsize=(10,5))
        ax = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        
        data_ref=results[data][0]
        data_exp=results[data][1]
        
        for i in range(len(data_ref)):
            x = data_ref.loc[data_ref.index[i],'fc_log2']
            y = data_ref.loc[data_ref.index[i],'-log_p']
            ax.plot(x, y, 'bo')
        
            ax.text(x*(1-0.02) , y*(1+0.02)  , data_ref.index[i], fontsize=4)
         
        for i in range(len(data_exp)):
            x = data_exp.loc[data_exp.index[i],'fc_log2']
            y = data_exp.loc[data_exp.index[i],'-log_p']
            ax2.plot(x, y, 'bo')
           
            ax2.text(x*(1-0.02) , y*(1+0.02)  , data_exp.index[i], fontsize=4)
        
        for axes in [ax,ax2]:
            axes.set_ylim((0, 10))
            axes.set_xlabel('fc_log2')
            axes.set_ylabel('-log_p')
            if axes==ax:
                axes.set_xlim((1, 5*fc_th))
                axes.set_title('Signif for WT-'+ filenames,fontsize=8)
            else:
                axes.set_xlim((-5*fc_th, -1))
                axes.set_title('Signif for mutant-'+ filenames,fontsize=8)
                
        if plotting==True:
            plt.show()
        if saving==True:
            fig.savefig(filenames + '.png',dpi=300,format='png',transparent=False)

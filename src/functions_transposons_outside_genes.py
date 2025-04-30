## Convert to a function the assembly of the positions_float_pd 
import pandas as pd
from collections import defaultdict
import numpy as np

def positions_float_dataframe(data,keys):
    """[summary]

    Parameters
    ----------
    data : dataframe
        A dataframe whose index are the standard names of the postprocessed dataset, excluding the classified as 1.

    Returns
    -------
    dataframe
        A dataframe with the locations of each gene in the genome as a list of floats, for further processing. 
    """

    from from_excel_to_list import from_excel_to_list


    positions_float=defaultdict(dict)
    for k in keys:
        data_new=data[data.loc[:,"background"]==k]

        for i in data_new.index:

            if type(data_new.loc[i,"Position"])==pd.core.series.Series:
                tmp=[]
                tmp_insertions=[]
                for j in np.arange(0,len(data_new.loc[i,"Position"])):
                    tmp.append(from_excel_to_list(data_new.loc[i,"Position"][j]))
                    positions_float[i,k]["Positions_float"]=tmp
                    tmp_insertions.append(data_new.loc[i,"Ninsertions"][j])
                    positions_float[i,k]["Ninsertions"]=tmp_insertions
            else:
                positions_float[i,k]["Positions_float"]=(from_excel_to_list(data_new.loc[i,"Position"]))
                positions_float[i,k]["Ninsertions"]=data_new.loc[i,"Ninsertions"]
            
            
            positions_float[i,k]["Feature_type"]=np.unique(data_new.loc[i,"Feature_type"])[0]
    
    positions_float_pd=pd.DataFrame.from_dict(positions_float,orient='index')

    return positions_float_pd


def defining_threshold_given_tr_density(data,windows_size,background):
    """[summary]

    Parameters
    ----------
    data : [type]
        [description]
    windows_size : [type]
        [description]
    key : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """    """"""

    threshold=np.round(data[data.loc[:,"background"]==background]["transposon_density"].unique()[0]*windows_size,decimals=0)
    density=data[data.loc[:,"background"]==background]["transposon_density"].unique()[0]
    return threshold,density


def get_discarded_genes_by_duplication(data_positions_pd,genes_names):
    """[summary]

    Parameters
    ----------
    data_positions_pd : [type]
        [description]
    genes_names : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """

    ## To discard genes that have multiple locations in the genome in the wt background

    discard_genes_duplicated=[]
    

    for i in genes_names:
        
        tmp=(data_positions_pd.loc[(data_positions_pd.background=="wt_merged") & (data_positions_pd["Gene name"]==i), "Positions_float"]).tolist()
        # tmp=positions_float_pd.loc[i,"Positions_float"][keys[0]]
        if tmp!=[]:
            if type(tmp[0][1])==list:
                discard_genes_duplicated.append(i)

    return discard_genes_duplicated


def get_amenable_genes_coverage_neighborhood(data_positions_pd,genes_names,discarded_genes_by_duplication,windows_size):
    """To get the genes that are inside the allow region of the genome based on the
    windows size chosen to look at the flanks for each gene. 
    For example if the windows size is 3kb then the genes whose initial genomic
    location is less than 3000 are discarded for the analsysis because there is not enough 
    data to make the calculation. Same if the end location is less than 3KB to the end of the genome. 

    Parameters
    ----------
    data_positions_pd : dataframe
        dataframe with the positions of the transposons along the genome.
    genes_names : list
        all considered genes
    discarded_genes_by_duplication : list
        duplictaed genes in the genome that are discarded

    Returns
    -------
    list
        target genes that are inside the allow region of the genome
    list 
        discarded genes that are outside the allow region of the genome
    """

    key="wt_merged"
    # extract the positions from the wt background 
    positions_float_pd_wt=data_positions_pd[data_positions_pd.loc[:,"background"]==key]

    # Renumber the indexes of the genes to start from 0
    positions_float_pd_wt.index=np.arange(0,len(positions_float_pd_wt))

    targets=[]
    valid_gene_names=np.setdiff1d(genes_names,discarded_genes_by_duplication)
    max_location=positions_float_pd_wt.loc[positions_float_pd_wt.index[-1],"Positions_float"][1]

    # keeping the genes whose initial location is above the windows size and the last one bellow
    # the end location - windows size 
    for names in valid_gene_names:
        if positions_float_pd_wt[positions_float_pd_wt.loc[:,"Gene name"]==names]["Positions_float"].tolist()!=[]:
            if positions_float_pd_wt[positions_float_pd_wt.loc[:,"Gene name"]==names]["Positions_float"].tolist()[0][1]<max_location-windows_size:
                targets.append(names)
            elif positions_float_pd_wt[positions_float_pd_wt.loc[:,"Gene name"]==names]["Positions_float"].tolist()[0][0]>windows_size:
                targets.append(names)

    genes_not_discarded_by_location=np.setdiff1d(valid_gene_names,targets)

    return targets,genes_not_discarded_by_location



def local_discrimination_genes_by_neighbors_coverage(data_positions_pd,background,gene_of_interest,windows_size,threshold):
    
    # extract the positions from the desired background
    data_background=data_positions_pd[data_positions_pd.loc[:,"background"]==background]

    # Renumber the indexes of the genes to start from 0
    data_background.index=np.arange(0,len(data_background))

    # Extract the index of the gene of interest
    index_target=np.where(data_background.loc[:,"Gene name"]==gene_of_interest)[0][0] 

    # Extract the genomic positions of the gene of interest
    tmp=(data_positions_pd.loc[(data_positions_pd.background==background) & (data_positions_pd["Gene name"]==gene_of_interest), "Positions_float"]).tolist()

    # Vectors that define the neighborhood of the target gene
    max_location=data_background.loc[data_background.index[-1],"Positions_float"][1]

    vector1=[tmp[0][0]-windows_size,tmp[0][0]]
    vector2=[tmp[0][1],tmp[0][1]+windows_size]
    if vector1[0]<0:
        vector1[0]=0
    if vector2[1]>max_location:
        vector2[1]=max_location

    # defining maximum upstream and downstream neighborhoods around the gene of interest

    neighborhood_upstream_insertions=[]
    neighborhood_upstream_insertions.append(data_background.loc[(data_background.index<index_target),"Ninsertions"])

    neighborhood_downstream_insertions=[]
    neighborhood_downstream_insertions.append(data_background.loc[(data_background.index>index_target),"Ninsertions"])

    neighborhood_upstream_locations=[]
    neighborhood_upstream_locations.append(data_background.loc[(data_background.index<index_target),"Positions_float"])

    neighborhood_downstream_locations=[]
    neighborhood_downstream_locations.append(data_background.loc[(data_background.index>index_target),"Positions_float"])
    
    ## Computing the downstream and upstream search by locations , collecting the insertions and compare with
    # the threshold
    upstream_insertions=[]
    discard_genes_upstream=[]
    for i in np.arange(0,index_target): # upstream search 
        if type(neighborhood_upstream_locations[0].tolist()[i][1])== float: # if the gene is not duplicated
            if neighborhood_upstream_locations[0].tolist()[i][1]<=vector1[0]:
                upstream_insertions.append(neighborhood_upstream_insertions[0].tolist()[i])
            
    sum_upstream_insertions=np.sum(upstream_insertions)

    if sum_upstream_insertions<threshold:
        discard_genes_upstream.append(gene_of_interest)

    discard_genes_downstream=[]
    downstream_insertions=[]
    downstream_search=data_background.index[-1]-index_target
    for i in np.arange(0,downstream_search,dtype=int): # downstream search 
        if type(neighborhood_downstream_locations[0].tolist()[i][1])== float: # if the gene is not duplicated
            if neighborhood_downstream_locations[0].tolist()[i][1]<=vector2[1]: 
                downstream_insertions.append(neighborhood_downstream_insertions[0].tolist()[i])

    sum_downstream_insertions=np.sum(downstream_insertions)

    if sum_downstream_insertions<threshold:
        discard_genes_downstream.append(gene_of_interest)

    

    return discard_genes_upstream,discard_genes_downstream,sum_upstream_insertions,sum_downstream_insertions


def get_genes_names(data_raw,keys):
    """[summary]

    Parameters
    ----------
    data_raw : [type]
        [description]
    keys : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """    """"""


    ## Refine data 
    data_raw.index=np.arange(0,len(data_raw))
    data_raw.drop(columns=['Unnamed: 1'],inplace=True)
    data_raw.fillna(0,inplace=True)
    data_raw.rename(columns={'Unnamed: 0':'background'},inplace=True)

    indexes_back=[] # indexes for the start of each background

    for i in keys:
        indexes_back.append(np.where(data_raw.loc[:,"background"]==i)[0])

    # Filling the name of the background according to the index
    for k in np.arange(0,len(indexes_back)-1):
        
        data_raw.loc[np.arange(indexes_back[k][0],indexes_back[k+1][0]),"background"]=keys[k]

    data_raw.loc[np.arange(indexes_back[-1][0],len(data_raw)),"background"]=keys[-1] # for the last key
        
    for key in keys:
        sum_tr=data_raw[data_raw.loc[:,"background"]==key]["Ninsertions"].sum()
        genome=data_raw[data_raw.loc[:,"background"]==key]["Nbasepairs"].sum()
        data_raw.loc[data_raw.loc[:,"background"]==key,"transposon_density"]=sum_tr/genome


    genes_names=data_raw[data_raw.loc[:,"Feature_type"]=="Gene; Verified"]["Standard_name"].unique()
    np.append(genes_names,data_raw[data_raw.loc[:,"Feature_type"]=="Gene; Dubious"]["Standard_name"].unique())
    np.append(genes_names,data_raw[data_raw.loc[:,"Feature_type"]=="Gene; Uncharacterized"]["Standard_name"].unique())
    np.append(genes_names,data_raw[data_raw.loc[:,"Feature_type"]=="Gene; Verified|silenced_gene"]["Standard_name"].unique())

    return data_raw,genes_names
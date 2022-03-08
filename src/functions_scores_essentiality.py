## Step 1: Selecting target genes for the essentiality analysis
import numpy as np
import pandas as pd
from collections import defaultdict
from from_excel_to_list import from_excel_to_list

def get_genes_names_for_essentiality(data_genes_names,data_discarded_genes_by_duplication,discarded_genes_by_low_coverage,background):
    """[summary]

    Parameters
    ----------
    data_genes_names : [type]
        [description]
    data_discarded_genes_names : [type]
        [description]
    background : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """

    a=np.array(data_genes_names.loc[background,"Gene name"])
    b=data_discarded_genes_by_duplication
    c=np.array(discarded_genes_by_low_coverage.loc[background,"discarded_genes_neighborhood"])

    genes_names_without_duplicates=np.setdiff1d(a,b)
    genes_names_without_low_coverage=np.setdiff1d(genes_names_without_duplicates,c)

    return genes_names_without_low_coverage


def get_no_duplicates_gene_names(data_genes_names,data_discarded_genes_by_duplication):
    """[summary]

    Parameters
    ----------
    data_genes_names : [type]
        [description]
    data_discarded_genes_names : [type]
        [description]
    background : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """

    a=np.array(data_genes_names.loc["wt_merged","Gene name"])
    b=data_discarded_genes_by_duplication

    genes_names_without_duplicates=np.setdiff1d(a,b)

    return genes_names_without_duplicates


def get_essentiality_score_per_gene_per_background(useful_genes,background,data_insertions):
    """[summary]

    Parameters
    ----------
    useful_genes : [type]
        [description]
    background : [type]
        [description]
    data_insertions : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """    """"""
    data_background=data_insertions.loc[background]
    data_background.index=data_background["Gene name"]
    valid_index=[x for x in data_background.index if x in useful_genes] # useful indexes

    score=defaultdict(dict)
    for gene in valid_index:
        insertion_float=from_excel_to_list(data_background.loc[gene,"Insertion locations"])

        if type(insertion_float)!=int : # for genes with more than one insertion location
            n = 5 # to ensure the number of locations is big enough, more than 5 
            if len(insertion_float)>n : # for genes which number of insertion locations bigger than 5

                ## Step 2: Compute largest interval free of transposons between the transposon n and the n+5
                
                x=np.array(insertion_float)
                #L=np.max(np.diff(insertion_float))
                
                L=np.max(x[n:] - x[:-n])

                ## Step 2: Compute length gene

                l=data_background.loc[gene,"End location"]-data_background.loc[gene,"Start location"]

                ## Step 3: Compute the number of insertions per gene 

                N=data_background.loc[gene,"Insertions"]
                #N=data_background.loc[gene,"tr_normalized_windows"]
                N_median=np.median(data_background.loc[:,"Insertions"])

                ## Compute the score per gene per background

                if all((N>N_median,L>200,L>0.1*l,L<0.9*l)):
                    #print(L,N,l)
                    score["value"][gene]=L*N/pow(l,1.5)
                    
                else:
                    score["value"][gene]=0
                    
                
                
            else:
                score["value"][gene]=0
        else:
            score["value"][gene]=0
    return score


def write_ones_if_essential(data_scores,background,essential_genes):
    """Writing 1 if the gene is essential , zero otherwise in the desired background

    Parameters
    ----------
    data_scores : [type]
        [description]
    background : [type]
        [description]
    essential_genes : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """

    scores_data=data_scores.loc[background,:]
    for true_genes in essential_genes:
        
        index_essential=np.argwhere(scores_data.index==true_genes)
        
        if len(index_essential)>0:
            #scores_wt.loc[index_essential[0],"true essential"]=1
            tmp=scores_data.index[index_essential[0][0]]
            scores_data.loc[tmp,"true essential"]=1
        
        
    scores_data.fillna(0,inplace=True)

    return scores_data

## Genes that only have insertions with one or two reads
def exclude_genes_with_one_or_two_reads(data,background):
    """[summary]

    Parameters
    ----------
    data : [type]
        [description]
    background : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """

    tmp=data.loc[background]
    tmp=tmp[tmp.loc[:,"Reads"]>2]["Gene name"] # genes with more than 2 reads
    
    return tmp
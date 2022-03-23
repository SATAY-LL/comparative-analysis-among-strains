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
    """Get the domain likelihood scores pointed out by Benoit in his
    paper https://elifesciences.org/articles/23570#fig2s1

    Parameters
    ----------
    useful_genes : array
        Array of genes names as str that have multiple genomic locations 
    background : str
        The name of the genetic background yhe user wish to analyze 
    data_insertions : pd.DataFrame
        The data where you can find the insertions locations per genomic locations and the 
        total number of insertions. 

    Returns
    -------
    dict
        The scores per genomic location
    """    

    ## Get the data from the background inserted by the user 
    data_background=data_insertions.loc[background]
    data_background.index=data_background["Gene name"]

    ## Get the name of the genes that are valid , that is , the one inserted in the useful_genes variable. 
    valid_index=[x for x in data_background.index if x in useful_genes] # useful indexes

    ## Compute the scores for each gene
    score=defaultdict(dict)
    longest_interval=defaultdict(dict)
    for gene in valid_index:
        # Get the insertion locations for each gene 
        insertion_float=from_excel_to_list(data_background.loc[gene,"Insertion locations"])
        # for genes with more than one insertion location

        if type(insertion_float)!=int : 
            # Setting a the minimum number of insertions a gene should have to compute a score
            n = 5 
            # for genes which number of insertion locations bigger than 5
            if len(insertion_float)>n : 

                ##Compute largest interval free of transposons between the transposon n and the n+5
                
                x=np.array(insertion_float)
                ### Shifting the array by n to the right and then to left and substract them and take the maximum
                L=np.max(x[n:] - x[:-n])
                longest_interval["value"][gene]=L


                ## Compute length gene

                l=data_background.loc[gene,"End location"]-data_background.loc[gene,"Start location"]

                ## Compute the number of insertions per gene 

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
    return score,longest_interval


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
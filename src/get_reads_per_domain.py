import pandas as pd
import numpy as np
from from_excel_to_list import from_excel_to_list
import matplotlib.pyplot as plt

## Function to get the reads per domain per gene 

def get_reads_per_domain(data,gene):
    """Get the reads per domain per gene and the genomic coordinates of the protein domains, assuming
    the gene does not contain introns. Thereby I can do a linear transformation from protein 
    domains into genomic locations.

    Parameters
    ----------
    data : dataframe
        Dataframe that contains the insertion locations and the reads per insertion location.
    gene : str
        The gene of interest to know the reads per domain . 

    Returns
    -------
    np.array
        Arrays with the reads per domain and the genomic coordinates of the protein domains
        of the gene of interest.
    """    

    

    genomic_length=data.loc[gene,"End location"]-data.loc[gene,"Start location"]
    protein_length=data.loc[gene,"protein length"]

    insertions_vector=from_excel_to_list(data.loc[gene,"Insertion locations"])
    reads_vector=from_excel_to_list(data.loc[gene,"Reads per insertion location"])

    
    domain_genomic=[]
    tmp=data.loc[gene,"domain locations"]
    
    if type(tmp)!=float : 

        if type(tmp.tolist()[0])!=int:
            
            for i in np.arange(0,int(len(np.concatenate(tmp).tolist())/2)):
            
                domain_protein=[data.loc[gene,"domain locations"][1][i],data.loc[gene,"domain locations"][0][i]]
                 
                if type(protein_length)==list: 
                    domain_protein_rel=np.array(domain_protein)/protein_length[0]
                else:
                    domain_protein_rel=domain_protein/protein_length

                domain_genomic_rel=np.round(domain_protein_rel*genomic_length)

                domain_genomic.append([domain_genomic_rel[1]+data.loc[gene,"Start location"],
        domain_genomic_rel[0]+data.loc[gene,"Start location"]]) 

        else:

                domain_protein=[data.loc[gene,"domain locations"][1],data.loc[gene,"domain locations"][0]]

                if type(protein_length)==list: 
                    domain_protein_rel=np.array(domain_protein)/protein_length[0]
                else:
                    domain_protein_rel=domain_protein/protein_length
               

                domain_genomic_rel=np.round(domain_protein_rel*genomic_length)

                domain_genomic.append([domain_genomic_rel[1]+data.loc[gene,"Start location"],
        domain_genomic_rel[0]+data.loc[gene,"Start location"]]) # Revert the array from start to end


    

    
    reads_domain=[]
    insertions_domain=[]
    for j in np.arange(0,len(domain_genomic)):
        pos_domain=[]
        pos_domain_hits=[]
        if type(insertions_vector)!=int:
            for i in np.arange(0,len(insertions_vector)):     

                if domain_genomic[j][0]<=insertions_vector[i] and insertions_vector[i]<=domain_genomic[j][1]:
                    pos_domain.append(reads_vector[i]) 
                    #pos_domain_hits.append(insertions_vector[i]) 
                    pos_domain_hits.append(1)  # append 1 if the insertion is in the domain

        reads_domain.append(np.sum(pos_domain))
        insertions_domain.append(np.sum(pos_domain_hits))

    return reads_domain,domain_genomic,insertions_domain



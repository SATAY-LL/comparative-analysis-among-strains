import numpy as np

from from_excel_to_list import from_excel_to_list
import wiggelen as wg
from collections import defaultdict

import pandas as pd


def reads_per_insertion_along_gene_length(data_pergene,background,number_of_parts=10):
    """Compute the reads per insertion along the gene length for a given background.
    It returns an array with the number of genes as rows and the number of columns the number
    of parts defined along the gene length.

    Parameters
    ----------
    data_pergene : pandas.DataFrame
        The pergene dataframe of all the backgrounds
    background : str
        The genetic background to be analyzed
    number_of_parts : int, optional
        The number of parts the user wants to divide the gene, by default 10

    Returns
    -------
    ndarray
        Array with the number of genes as rows and the number of columns the number
        of parts defined along the gene length.
    """

    data=data_pergene.loc[background]

    ## Extract as lists the gene coordinates, and the locations of every insertion and reads per gene 
    reads_locations=[]
    insertion_locations=[]
    gene_coordinates=[]


    for j in data.index:
            coi=from_excel_to_list(data.loc[j]["Reads per insertion location"])
            coi_insertions=from_excel_to_list(data.loc[j]["Insertion locations"])
            gene_coordinates.append([data.loc[j,"Start location"],data.loc[j,"End location"]])
            reads_locations.append(coi)
            insertion_locations.append(coi_insertions)

    # compute the reads per insertions for each gene along the gene length
    gene_parts=np.linspace(0,1,number_of_parts+1)
    r=np.zeros(shape=(len(insertion_locations),len(gene_parts))) # reads_per_insertion_parts array , every gene in the rows 

    for i in np.arange(0,len(insertion_locations)):
        if (insertion_locations[i])!=0:
            g=np.array(insertion_locations[i]) # insertion locations
            f=np.linspace(gene_coordinates[i][0],gene_coordinates[i][1],len(gene_parts)) # gene coordinates
            #f=np.array(gene_coordinates[i][0]+gene_coordinates[i][1]*gene_parts)
            binedges = g.searchsorted(f)

            rngs = [list(range(binedges[k], binedges[k+1])) for k in range(len(binedges)-1)]

            for k in np.arange(0,len(rngs)):
                readsperinsert=[]
                for j in np.arange(0,len(rngs[k])):
                    readsperinsert.append(reads_locations[i][rngs[k][j]])
                    
                if len(readsperinsert)>1:
                    r[i,k]=np.sum(readsperinsert)/(len(reads_locations[i])-1)#discarding the insertion with the highest read count
                else:
                    r[i,k]=0

    return r,gene_coordinates,reads_locations,insertion_locations




def genes_discarded4fitness_from_flanking_regions(data_pergene,background,wigfile_path,windows_size=20000):
    """_summary_

    Parameters
    ----------
    data_pergene : _type_
        _description_
    background : _type_
        _description_
    wigfile_path : _type_
        _description_
    windows_size : int, optional
        _description_, by default 20000

    Returns
    -------
    _type_
        _description_
    """    

 

    wig_reads=[]

    for x in wg.walk(open(wigfile_path)):
        wig_reads.append(x)

    wig_reads = np.array(wig_reads)

    ## reading wig file
    pos=[]
    chrom=[]
    reads=[]
    for i in wig_reads:
        chrom.append(i[0]) ## chromosomes
        pos.append(float(i[1])) ## positions
        reads.append(float(i[2])) ## reads

    chrom=np.array(chrom)
    pos=np.array(pos)
    reads=np.array(reads)

    # size of every chromosome in bp from https://www-ncbi-nlm-nih-gov.tudelft.idm.oclc.org/genome/?term=Saccharomyces%20cerevisiae%5BOrganism%5D&cmd=DetailsSearch
    chromosome_size={"I":230218,"II":813184,"III":316620,
    "IV":1531933,"V":576874,"VI":270161,"VII":1090940,
    "VIII":562643,"IX":439888,"X":745751,"XI":666816,
    "XII":1078177,"XIII":924431,"XIV":784333,"XV":1091291,
    "XVI":948066}

    chrom2latin={"I":1,"II":2,"III":3,"IV":4,"V":5,"VI":6,"VII":7,"VIII":8,"IX":9,
    "X":10,"XI":11,"XII":12,"XIII":13,"XIV":14,"XV":15,"XVI":16}

    chrom_conversion={"chrI":"I","chrII":"II","chrIII":"III",
    "chrIV":"IV","chrV":"V","chrVI":"VI","chrVII":"VII","chrVIII":"VIII",
    "chrIX":"IX","chrX":"X","chrXI":"XI","chrXII":"XII","chrXIII":"XIII",
    "chrXIV":"XIV","chrXV":"XV","chrXVI":"XVI"}

    chromosome_density=[]

    for i in chrom_conversion.keys():
        tmp=len(pos[chrom==i])/chromosome_size[chrom_conversion[i]]
        chromosome_density.append(tmp)

    expected_flanking_enrichment=np.round(np.multiply(chromosome_density,windows_size),0).astype(int)# this means 5kb flanking regions each side

        
    data=data_pergene.loc[background]


    flanking_regions_data=defaultdict(dict)
    for gene in data["Gene name"]:

        data_loc_gene=data[data.loc[:,"Gene name"]==gene]
        data_loc=[data_loc_gene["End location"],data_loc_gene["Start location"]]

        gene_chromosome=data_loc_gene["Chromosome"].tolist()[0]

        if gene_chromosome!="Mito":

            
            # flanking regions
            flanking_regions=[data_loc[1].tolist()[0]-windows_size/2,data_loc[0].tolist()[0]+windows_size/2]
            pos_chrom=pos[chrom=="chr"+gene_chromosome]


            ## find the closest index of the flanking regions in the chromosome insertions from the wig file
            pos_up=min(range(len(pos_chrom)), key=lambda i: abs(pos_chrom[i]-flanking_regions[0]))
            pos_down=min(range(len(pos_chrom)), key=lambda i: abs(pos_chrom[i]-flanking_regions[1]))

            insertions_flanking=len(pos_chrom[pos_up:pos_down])

            flanking_regions_data[gene]["gene location"]=data_loc
            flanking_regions_data[gene]["chromosome"]=gene_chromosome
            flanking_regions_data[gene]["flanking_regions"]=flanking_regions
            flanking_regions_data[gene]["insertions_flanking"]=insertions_flanking
            flanking_regions_data[gene]["expected_flanking_enrichment"]=expected_flanking_enrichment[chrom2latin[gene_chromosome]-1] # it should be 2* expected value but actually many insertions goes to centromeres and then we are overestimating this number for the flanking regions of any gene
            if insertions_flanking<0.4*flanking_regions_data[gene]["expected_flanking_enrichment"]: # astringent cut off (0.5 is an intermediate cut off)
                flanking_regions_data[gene]["classification"]="Not enough flanking regions"
            else:
                flanking_regions_data[gene]["classification"]="OK"

    flanking_regions_data_pd=pd.DataFrame.from_dict(flanking_regions_data,orient="index")

    discarded_genes2fitness=flanking_regions_data_pd[flanking_regions_data_pd["classification"]=="Not enough flanking regions"].index

    return flanking_regions_data_pd,discarded_genes2fitness




def protein_domains_info(data_pergene,background,gene_coordinates,reads_locations,insertion_locations,discarded_genes):
    """_summary_

    Parameters
    ----------
    data_pergene : _type_
        _description_
    background : _type_
        _description_
    gene_coordinates : _type_
        _description_
    reads_locations : _type_
        _description_
    insertion_locations : _type_
        _description_
    discarded_genes : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    ## import protein domains
    ## data generated in the notebook "analysis_proteins_domains.ipynb"
    data_domains=pd.read_excel("../postprocessed-data/genomic-domains-wt.xlsx",index_col="Unnamed: 0")


    ## data from yeastmine
    domains_names=pd.read_csv('../data/Domains_all_genes_protein_coordinates_yeastmine.tsv',sep="\t")
    domains_names.index=domains_names["Gene Name"]

    protein_domains_data=defaultdict(dict)

    data=data_pergene.loc[background]

    

    for i in data.index:
        gene=data.loc[i,"Gene name"]
        if gene in data_domains.index:
            
            
            protein_domains_data[gene]["gene coordinates"]= np.round(np.linspace(gene_coordinates[i][0],
            gene_coordinates[i][1],10),0)
            central_gene_part=protein_domains_data[gene]["gene coordinates"][1:9]
            protein_domains_data[gene]["gene coordinates central"]=central_gene_part
            
            if data_domains.loc[gene][0]!="[]":
                tmp=from_excel_to_list(data_domains.loc[gene][0])
                protein_domains_data[gene]["domains coordinates"]= tmp
                tmp_central=[]
                for i in np.arange(0,len(tmp),2):
                    if np.min(central_gene_part)<tmp[i] and np.max(central_gene_part)>tmp[i]: # if the first protein coordinate from the domain is inside the central part
                        tmp_central.append(tmp[i:i+2])
                    elif np.min(central_gene_part)<tmp[i+1] and np.max(central_gene_part)>tmp[i+1]: # if the first one is not but the second one is
                        tmp_central.append(tmp[i:i+2])
                protein_domains_data[gene]["domains coordinates central"]=tmp_central
                        
                protein_domains_data[gene]["domains name"]= np.array(domains_names.loc[gene]["Protein Domain"])
                protein_domains_data[gene]["domain description"]= np.array(domains_names.loc[gene]["Protein domain description"])
                protein_domains_data[gene]["reads_domain"]= np.zeros_like(protein_domains_data[gene]["domains name"])
                protein_domains_data[gene]["insertions_domain"]= np.zeros_like(protein_domains_data[gene]["domains name"])
            else:
                protein_domains_data[gene]["domains coordinates"]= np.nan
                protein_domains_data[gene]["domains name"]= np.nan
                protein_domains_data[gene]["domain description"]= np.nan
                protein_domains_data[gene]["reads_domain"]= np.sum(reads_locations[i])
                protein_domains_data[gene]["insertions_domain"]= insertion_locations[i]

            if gene in discarded_genes:
                protein_domains_data[gene]["classification"]="Not enough flanking regions"
            else:
                protein_domains_data[gene]["classification"]="OK"
        else:
            protein_domains_data[gene]["gene coordinates"]=np.nan
            protein_domains_data[gene]["gene coordinates central"]=np.nan
            protein_domains_data[gene]["domains coordinates"]= np.nan
            protein_domains_data[gene]["domains name"]= np.nan
            protein_domains_data[gene]["domain description"]= np.nan
            protein_domains_data[gene]["reads_domain"]= np.nan
            protein_domains_data[gene]["insertions_domain"]= np.nan
            

            
    protein_domains_data_pd=pd.DataFrame.from_dict(protein_domains_data,orient="index")

    return protein_domains_data_pd



def reads_and_insertions_per_domain(data_pergene,background,data_domains,reads_locations,insertion_locations):
    """_summary_

    Parameters
    ----------
    data_pergene : _type_
        _description_
    background : _type_
        _description_
    data_domains : _type_
        _description_
    reads_locations : _type_
        _description_
    """

    data=data_pergene.loc[background]
    domain_coordinates=data_domains["domains coordinates"]
    domain_coordinates=domain_coordinates[domain_coordinates.isna()==False]

    for i in domain_coordinates.index: # over all genes 
        i_index=data_domains.index.get_loc(i)
       
        b=np.array(insertion_locations[i_index])
        if b.size>1:
            f=domain_coordinates.loc[i] # protein domains
            binedges = b.searchsorted(f) # find the index of the insertion in the domains

            rngs = [list(range(binedges[k], binedges[k+1])) for k in range(0,len(binedges)-1,2)]

            totalreads=[]
            totalinsertions=[]
            for k in np.arange(0,len(rngs)):
                readsperdomain=[]
                insertionsperdomain=[]
                for j in np.arange(0,len(rngs[k])):
                    tmp=reads_locations[i_index][rngs[k][j]]
                    readsperdomain.append(tmp)
                    if type(tmp)!=float:
                        insertionsperdomain.append(len(tmp))
                    else:
                        insertionsperdomain.append(1)

                totalreads.append(np.sum(readsperdomain))
                totalinsertions.append(np.sum(insertionsperdomain))
        else:
            totalreads=1
            totalinsertions=1
                    
        data_domains.loc[i,"reads_domain"]=totalreads
        data_domains.loc[i,"insertions_domain"]=totalinsertions


    return data_domains


def fitness_models(data_pergene,background,data_domains_extended,reads_per_insertion_array,discarded_genes):
    """_summary_

    Parameters
    ----------
    data_pergene : _type_
        _description_
    background : _type_
        _description_
    data_domains_extended : _type_
        _description_
    reads_per_insertion_array : _type_
        _description_
    discarded_genes : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """    

    ref=np.log2(np.median(np.sum(reads_per_insertion_array[1:9],axis=1))) # reference fitness, assumption: most genes are neutral in the wild type
    
    data=data_pergene.loc[background]
    fitness_models=defaultdict(dict)



    for i in np.arange(0,len(data)):
        gene=data.loc[i,"Gene name"]
        if gene in discarded_genes:
            fitness_models[gene]["fitness_gene"]="Not enough flanking regions"
            fitness_models[gene]["fitness_gene_std"]="Not enough flanking regions"
            fitness_models[gene]["fitness_domains_vector"]="Not enough flanking regions"
            fitness_models[gene]["fitness_domains_average"]="Not enough flanking regions"
            fitness_models[gene]["fitness_domains_corrected"]="Not enough flanking regions"
            fitness_models[gene]["domains"]="Not enough flanking regions"
        
        else:
            if np.sum(reads_per_insertion_array[i,1:9])!=0: # if there are no reads in the 80% central part of the gene, the fitness is not calculated
                fitness_models[gene]["fitness_gene"]=np.log2(np.sum(reads_per_insertion_array[i,1:9]))/ref # getting the 80% central part of the reads per insertions
                fitness_models[gene]["fitness_gene_std"]=(np.std(reads_per_insertion_array[i,1:9]))
                if type(data_domains_extended.loc[gene,"reads_domain"])==list:
                    nume=np.array(data_domains_extended.loc[gene,"reads_domain"])
                    deno=np.array(data_domains_extended.loc[gene,"insertions_domain"])
                    H=data_domains_extended.loc[gene,"exclude domains"]
                    
                    if len(H)==1:
                        if H[0]==True:
                            fitness_models[gene]["fitness_domains_vector"]="Not enough insertions"
                        elif H[0]==False and nume[0]==0:
                            fitness_models[gene]["fitness_domains_vector"]="Not enough insertions"
                        elif H[0]==False and deno[0]==0:
                            fitness_models[gene]["fitness_domains_vector"]="Not enough insertions"
                        elif H[0]==False and deno[0]!=0:
                            fitness_models[gene]["fitness_domains_vector"]=np.log2(nume/deno)/ref
                    else:
                        f=[]
                        for j in np.arange(0,len(H)):
                            
                            if H[j]==False and deno[j]!=0:
                                y=np.log2(nume[j]/deno[j])/ref
                                
                            elif H[j]==True : # the domain do not have enough insertions to compute fitness 
                                y="Not enough insertions"
                            elif H[j]==False and deno[j]==0:
                                y=0
                            
                            f.append(y)
                            fitness_models[gene]["fitness_domains_vector"]=f
                    if "Not enough insertions"==fitness_models[gene]["fitness_domains_vector"]:
                        fitness_models[gene]["fitness_domains_average"]="Not enough insertions"
                        fitness_models[gene]["fitness_domains_std"]="Not enough insertions"
                        fitness_models[gene]["domains"]="Not enough insertions"
                        fitness_models[gene]["fitness_domains_corrected"]="Not enough insertions"
                    elif "Not enough insertions"  in fitness_models[gene]["fitness_domains_vector"]: 
                        J=np.where(np.array(fitness_models[gene]["fitness_domains_vector"])=="Not enough insertions")[0]
                        # exclude J from the fitness vector
                        f_0=np.delete(fitness_models[gene]["fitness_domains_vector"],J)
                        
                        f_0=np.array(f_0,dtype=float)
                        if len(f_0)>0:
                            for i in J:
                                f[i]=0 # assign zero fitness if the domain that does not have enough insertions is inside other that have enough insertions
                            fitness_models[gene]["fitness_domains_average"]=np.mean(f)
                            fitness_models[gene]["fitness_domains_std"]=np.std(f)
                            fitness_models[gene]["domains"]="annotated"
                            ## computing the corrected fitness as the fitness domain that has the maximum difference with the average fitnes
                            tmp=np.absolute(f-fitness_models[gene]["fitness_gene"])
                            index_max=np.where(tmp==np.max(tmp))[0]
                            if len(index_max)==1:
                                fitness_models[gene]["fitness_domains_corrected"]=f[int(index_max)]
                            else: # just take the first domain as the one with the maximum difference
                                fitness_models[gene]["fitness_domains_corrected"]=f[int(index_max[0])]
                        else:
                            fitness_models[gene]["fitness_domains_average"]="Not enough insertions"
                            fitness_models[gene]["fitness_domains_std"]="Not enough insertions"
                            fitness_models[gene]["domains"]="Not enough insertions"
                            fitness_models[gene]["fitness_domains_corrected"]="Not enough insertions"
                    else:
                        fitness_models[gene]["fitness_domains_average"]=np.mean(fitness_models[gene]["fitness_domains_vector"])
                        fitness_models[gene]["fitness_domains_std"]=np.std(fitness_models[gene]["fitness_domains_vector"])
                        fitness_models[gene]["domains"]="annotated"
                        ## computing the corrected fitness as the fitness domain that has the maximum difference with the average fitnes
                        tmp=np.absolute(fitness_models[gene]["fitness_domains_vector"]-fitness_models[gene]["fitness_gene"])
                        index_max=np.where(tmp==np.max(tmp))[0]
                        if len(index_max)==1:
                            fitness_models[gene]["fitness_domains_corrected"]=fitness_models[gene]["fitness_domains_vector"][int(index_max)]
                        else: # just take the first domain as the one with the maximum difference
                            fitness_models[gene]["fitness_domains_corrected"]=fitness_models[gene]["fitness_domains_vector"][int(index_max[0])]
                else:
                    fitness_models[gene]["fitness_domains_vector"]=fitness_models[gene]["fitness_gene"]
                    fitness_models[gene]["fitness_domains_std"]=fitness_models[gene]["fitness_gene_std"]
                    fitness_models[gene]["fitness_domains_average"]=fitness_models[gene]["fitness_gene"]
                    fitness_models[gene]["fitness_domains_corrected"]=fitness_models[gene]["fitness_gene"]
                    fitness_models[gene]["domains"]="non annotated domains"
            else:
                fitness_models[gene]["fitness_gene"]=0
                fitness_models[gene]["fitness_gene_std"]=0
                fitness_models[gene]["fitness_domains_vector"]=0
                fitness_models[gene]["fitness_domains_average"]=0
                fitness_models[gene]["fitness_domains_corrected"]=0
                fitness_models[gene]["domains"]="zero reads"


    fitness_models_pd=pd.DataFrame.from_dict(fitness_models,orient="index")
    # fitness_models_pd["fitness_domains_average"].replace(-np.inf,0,inplace=True)
    # fitness_models_pd["fitness_domains_corrected"].replace(-np.inf,0,inplace=True)
    # fitness_models_pd["fitness_domains_average"].replace(np.inf,1.5,inplace=True)
    # fitness_models_pd["fitness_domains_corrected"].replace(np.inf,1.5,inplace=True)
    # fitness_models_pd["fitness_gene"].replace(-np.inf,0,inplace=True)
    # fitness_models_pd["fitness_gene"].replace(np.inf,1.5,inplace=True)



    return fitness_models_pd


def excluding_domains(data_pergene,background,data_domains_extended):
    
    data=data_pergene.loc[background]
    data_domains_corrected=defaultdict(dict)
    k=0
    for gene in data_domains_extended.index:
        data_domains_corrected[gene]["transposon density"]=data.loc[k,"Insertions"]/(data.loc[k,"End location"]-data.loc[k,"Start location"])

        x=data_domains_extended.loc[gene,"domains coordinates"]
        if type(x)!=float:
            
            # substract the  even elements by the odd elements
            x=np.array(x)

            x1=x[1::2]
            x2=x[::2]

            x3=x1-x2 ## getting the length of every domain
            exclude_dom_all=[]
            for i in np.arange(0,len(x3)):
                
                if data_domains_corrected[gene]["transposon density"]*x3[i]<5:
                    exclude_dom=True
                else:
                    exclude_dom=False
                exclude_dom_all.append(exclude_dom)
            data_domains_corrected[gene]["exclude domains"]=exclude_dom_all
            data_domains_corrected[gene]["domain length"] = x3
    
        k=k+1
    ## compute average density
    density=[]
    for i in np.arange(0,len(data)):
        length=data.loc[i,"End location"]-data.loc[i,"Start location"]
        density.append(data.loc[i,"Insertions"]/length)

    average_density=np.median(density)

    
    data_domains_corrected_pd=pd.DataFrame.from_dict(data_domains_corrected,orient="index")

    ## adding that dataframe to the data_domains_extended dataframe

    data_domains_extended["exclude domains"]=data_domains_corrected_pd["exclude domains"]
    data_domains_extended["transposon density"]=data_domains_corrected_pd["transposon density"]
    data_domains_extended["domain length"]=data_domains_corrected_pd["domain length"]
    data_domains_extended["average density library"]=average_density

    return data_domains_extended
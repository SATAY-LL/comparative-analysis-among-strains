# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 15:16:28 2021

@author: ekingma1/leilaicruz
"""
import pandas as pd
import numpy as np
import operator
from  from_excel_to_list import from_excel_to_list
from scipy import stats
from collections import defaultdict
import ast


# from insertion_bias.corr_func import power_law

def profile(tn_file,step_size=10000,
            chr_file = 'info_files/genetic/chromosome_sizes.txt',
            cen_file = 'info_files/genetic/centromere_locations.txt'):
    """
    

    Parameters
    ----------
    tn_file : TYPE
        DESCRIPTION.
    step_size : TYPE, optional
        DESCRIPTION. The default is 10000.
    chr_file : TYPE, optional
        DESCRIPTION. The default is 'info_files/genetic/chromosome_sizes.txt'.
    cen_file : TYPE, optional
        DESCRIPTION. The default is 'info_files/genetic/centromere_locations.txt'.

    Returns
    -------
    d_cen : TYPE
        DESCRIPTION.
    tn_csum : TYPE
        DESCRIPTION.

    """
    
    # Importing from a file if the file is in wiggle format
    if type(tn_file) == str and tn_file.split(".")[-1] == "wig":
        df_tn = wiggle(tn_file)
    
    # Use tn_file variable directly if data is presented as DataFrame    
    elif type(tn_file) == pd.DataFrame: 
        df_tn = tn_file
        
    # Exit if data is not presented in one of the correct formats
    else:
        exit(print("The variable tn_file is not properly defined"))
       
    tn_pos = np.array( df_tn["Position"] )
    tn_chrom = np.array( df_tn["Chrom"] )
     
    # Import information about the length of the chromosomes
    chr_sizes = pd.read_csv(chr_file,sep="\t")
    longest_chr = max(chr_sizes['Length (bp)'])
    
    # Import the location of the centromeres
    cen = pd.read_csv(cen_file,sep="\t")
    
    # Pre-initialize the array for the cumulatitive sum of insertions
    n_cols = int( np.ceil(longest_chr/step_size) )
    tn_csum = np.empty( ((len(chr_sizes)-1)*2,n_cols) )
    tn_csum[:] = np.nan
    
    # Determine the distances over which the cumulatative sum is calculated
    d_cen = np.arange( 0,n_cols*step_size,step_size )
    
    # Iterate over all chromosomes
    count = 0
    for i in cen.index:    
        # Get the size of the chromosome and the start and stop of the centromere
        chrom_size = chr_sizes.loc[i,'Length (bp)']
        cen_start = cen.loc[i,'start']
        cen_stop = cen.loc[i,'stop']
        
        # Create an array of the positions to bin over
        pos_upstr = np.arange( 0,cen_start,step_size )
        pos_down = np.arange( cen_stop,chrom_size,step_size ) - cen_stop
        
        # Get the insertions that map to the current chromosome
        ins_chrom = tn_pos[ tn_chrom == cen.loc[i,'chrom'] ]
        
        # Calculate the distance of each insertion site to the centromere
        d_upstr = cen_start - ins_chrom[ ins_chrom<cen_start ]
        d_downstr = ins_chrom[ ins_chrom>cen_stop ] - cen_stop
        
        # Sort the distance arrays in ascending order
        d_upstr.sort(); d_downstr.sort()
        
        # Calculate the cumulative number of insertions with the distance to the centromere
        csum_up = [ len( d_upstr[ 0:np.argmin(d_upstr<=n) ] ) for n in pos_upstr 
                       if len(d_upstr)!=0 and np.argmin(d_upstr<=n)!=0 ]
        
        csum_down = [ len( d_downstr[ 0:np.argmin(d_downstr<=n) ] ) for n in pos_down
                        if len(d_downstr)!=0 and np.argmin(d_downstr<=n)!=0 ] 
        
        # Store the sums of different chromosomes in a single array
        tn_csum[ count,0:len(csum_up)+1 ] = [0]+csum_up
        tn_csum[ count+1,0:len(csum_down) +1] = [0]+csum_down
        
        count+=2
        
    return d_cen, tn_csum, cen


def correct(d,tn_csum,maxdist=200000,method="poly",deg=4): 
    """
    

    Parameters
    ----------
    d : TYPE
        DESCRIPTION.
    tn_csum : TYPE
        DESCRIPTION.
    maxdist : TYPE, optional
        DESCRIPTION. The default is 300000.
    method : TYPE, optional
        DESCRIPTION. The default is "plaw".
    deg : TYPE, optional
        DESCRIPTION. The default is 4.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    # Remove columns with only nan values from the tn cumulative matrix
    last_nan = max( np.where( ~np.isnan(tn_csum) )[1] )+1
    tn_csum_nonNaN = tn_csum[:,0:last_nan]
    d = d[0:last_nan]
    
    # Initiallize the array for the average cumulative tn insertions    
    tn_csum_av = np.nanmean(tn_csum_nonNaN,axis=0)
    
    # Correct the centromere bias by fitting a power-law function
    if method == "plaw": 
        d_log = np.log( d[tn_csum_av>0] )
        tn_csum_log = np.log( tn_csum_av[tn_csum_av>0] ) 
        
        p,rss = np.polyfit(d_log, tn_csum_log, 1,full = True)[0:2] 
        tn_rate = power_law( p[0],np.exp(p[1]) )
    
        return tn_rate
    
    if method == "poly":
        p = np.polyfit(d[d<maxdist],tn_csum_av[d<maxdist],deg)
        s_fun = np.poly1d(p)
        
        return s_fun
    


#text2array = data_format_conversion.text2array

def per_gene(datafile,gene=[],remove_ade=True):
    
    data = pd.read_csv(datafile, sep='\t', header=0, index_col=0) 
    
    if remove_ade:
        
        data = data.drop("ADE2")
    
    if not gene:
        
        data['Reads per insertion location'] =\
            data['Reads per insertion location'].apply(from_excel_to_list)
            
        data['Insertion locations'] =\
            data['Insertion locations'].apply(from_excel_to_list)
    
    else:
        
        data = data.loc[gene]
        data['Reads per insertion location'] =\
            data['Reads per insertion location'].apply(from_excel_to_list)
            
        data['Insertion locations'] =\
            data['Insertion locations'].apply(from_excel_to_list)
            
    return data


def wiggle(wig_file,remove_ade=True,as_dict=False):
    
    # Predefine the lists
    N_reads=[]
    Chrom=[]
    pos=[]
    
    # Walk through the wig file to extract the data   
    w = open(wig_file)
    for line in w:
        
        if "chrom" in line:      
            chrm = line.split("=")[-1].strip()
            
        elif "".join( line.split() ).isdigit():
            pos.append( int(line.split()[0]) )
            N_reads.append( int(line.split()[1]) )
            Chrom.append( chrm )
    
    # Store the data in a dataframe
    out_data=pd.DataFrame({'Reads':N_reads,'Position':pos,'Chrom':Chrom})
    
    # (OPTIONAL) Remove all reads mapping to the ADE2 locus. Default is True
    if remove_ade: 
        ade2_reads= np.logical_and(out_data['Chrom']=='chrXV',
                    np.logical_and(out_data['Position']>564476,
                                   out_data['Position']<566191))
        out_data = out_data.loc[map(operator.not_, ade2_reads)]
    
    # (OPTIONAL) Convert the data format into a dictionary for each chromosome
    if as_dict:      
        out_data_dict = {}
        chrom_names = out_data['Chrom'].unique()
        
        for i in chrom_names:
            chrom_indx = out_data['Chrom'] == i
            chrom_insertions = out_data.loc[chrom_indx,['Reads','Position']]
            out_data_dict[i] = chrom_insertions
            
        out_data = out_data_dict
        
    return out_data


## output the index of a vector whose value shows the maximum difference with a number
def max_difference_index(vector, number):
    # Calculate the element-wise absolute differences
    differences = np.abs(vector - number)
    
    # Find the index where the difference is maximum
    max_index = np.argmax(differences)
    
    return max_index


# Ensure the entire structure is converted to floats
def convert_to_float(item):
    if isinstance(item, list):
        return [convert_to_float(elem) for elem in item]
    else:
        return float(item)
    

# +
# Correct for the centromere bias in the insertion density
def fitness_average(pergene_files,wig_inserts,gene_inserts):

    maxdist = 200000
    chr_file="chromosome_sizes.txt"
    cen_file="centromere_locations.txt"
    fitness_all=[]

    for file in np.arange(len(pergene_files)):
        c_dis, tn_csum, cen = profile(wig_inserts[file],chr_file=chr_file,cen_file=cen_file)
        poly_fit = correct(c_dis, tn_csum,method="poly",maxdist=maxdist)
        poly_der = np.polyder(poly_fit, 1)

        # Estimate the expected number of insertions for each gene
        exp_tn = []
        reads_av = []
        gene_name = []
        sem=[]

        for g in gene_inserts[file].index:
            
            # Get the info of the current gene
            gene_info = gene_inserts[file].loc[g]
            gene_size = gene_info['End location'] - gene_info['Start location']
            
            reads_b =  ( wig_inserts[file]["Chrom"] == "chr"+gene_info["Chromosome"] ) & \
                    ( wig_inserts[file]["Position"] > gene_info['Start location']+gene_size*0.1 ) & \
                    ( wig_inserts[file]["Position"] < gene_info['Start location']+gene_size*0.9 )
            
            reads = wig_inserts[file]["Reads"][reads_b]
            
            # Calculate smallest distance from the gene to the centromere
            curr_cen = cen[ "chr"+gene_info["Chromosome"] == cen["chrom"] ].reset_index()
            
            if len(curr_cen):
                
                gene_center = gene_info['Start location'] + gene_size/2
                d_cen = min( abs( gene_center - curr_cen["start"][0] ), 
                            abs( gene_center - curr_cen["stop"][0] ) )
                
                # Determine the expected number of transposon insertions
                if d_cen <= maxdist:           
                    E_tn = np.floor( poly_der(d_cen) *gene_size *0.8 ) 
                    
                else:
                    E_tn = np.floor( poly_der(maxdist) *gene_size *0.8 )
                
                # Add zeros to the reads if the expected tn count is higher than observed
                if E_tn > len(reads):
                    
                    tn_add = np.zeros( int( E_tn-len(reads) ) )
                    reads = np.append(reads,tn_add)
                        
                # Remove outliers using the interquantile range
                iqr = stats.iqr(reads, axis=None, rng=(5, 95), scale=1.0)
                
                if any(reads):
                    perc = np.percentile(reads,95)
                
                if iqr>0:
                    reads_filt = reads[ reads < perc+1.5*iqr ]
                    
                else:
                    reads_filt = reads
            
                # Store gene name when the number of tn is below the threshold
                if len(reads_filt) < 5 or sum(reads_filt)==0:
                    
                    reads_av.append(np.nan)
                    gene_name.append( gene_info.name )
                    sem.append(np.nan)
                    continue
                
                # Store the average read count for the mean fitness
                reads_av.append( np.mean( reads_filt ) )
                gene_name.append( gene_info.name )
                sem.append(np.sqrt(np.std(reads_filt)/len(reads_filt)))
                    
        # Construct fitness array 
        data = np.empty( ( len(reads_av),2 ) )
        data[:,0] = [np.log2(x) for x in reads_av]
        data[:,1]=sem # standardr error of the mean
        fitness_average = pd.DataFrame(data,columns=['mean','sem'],index=gene_name)

        fitness_norm=fitness_average["mean"]/np.nanmedian(fitness_average["mean"])
        fitness_average["mean_norm"]=fitness_norm
            
        fitness_all.append(fitness_average)
    # -

        fitness_all_pd=pd.concat([fitness_all[i] for i in np.arange(0,len(fitness_all))],axis=0)

        return fitness_all_pd
    

def fitness_domains(pergene_files,wig_inserts,gene_inserts,data_domains):

    fitness_domain=defaultdict(dict)
    fitness_domain_all=[]
    cen_file="centromere_locations.txt"
    chr_file="chromosome_sizes.txt"

    maxdist = 200000
    for file in np.arange(0,len(pergene_files)): 
        c_dis, tn_csum, cen = profile(wig_inserts[file],chr_file=chr_file,cen_file=cen_file)
        poly_fit = correct(c_dis, tn_csum,method="poly",maxdist=1000000)
        poly_der = np.polyder(poly_fit, 1)

        # Estimate the expected number of insertions for each gene
        exp_tn = []
        reads_av = []
        gene_name = []

        # common indexes between domains_names and gene_inserts
        common_indexes=list(set(data_domains.index).intersection(set(gene_inserts[file].index)))
        for g in common_indexes:

            # Get the domains from Pfam of the current gene 
        
            x=convert_to_float(ast.literal_eval(data_domains.loc[g].tolist()[0])) 
            
                

            # Get the info of the current gene
            gene_info = gene_inserts[file].loc[g]
            gene_size = gene_info['End location'] - gene_info['Start location']
            

            # Loop over the domains of x
            fitness_perdomain=[]
            sem=[]
            for i in np.arange(0,len(x)):

                reads_b =  ( wig_inserts[file]["Chrom"] == "chr"+gene_info["Chromosome"] ) & \
                    ( wig_inserts[file]["Position"] > int(x[i][0]) ) & \
                    ( wig_inserts[file]["Position"] < int(x[i][1]) )
                reads = wig_inserts[file]["Reads"][reads_b]
                domain_size=int(x[i][1])-int(x[i][0])

                # Calculate smallest distance from the gene to the centromere
                curr_cen = cen[ "chr"+gene_info["Chromosome"] == cen["chrom"] ].reset_index()
                
                if len(curr_cen):
                    
                    gene_center = gene_info['Start location'] + gene_size/2
                    d_cen = min( abs( gene_center - curr_cen["start"][0] ), 
                                abs( gene_center - curr_cen["stop"][0] ) )
                    
                    # Determine the expected number of transposon insertions
                    if d_cen <= maxdist:           
                        E_tn_domain = np.floor( poly_der(d_cen) *domain_size ) 
                        E_tn_gene = np.floor( poly_der(d_cen) *gene_size *0.8 )
                        
                    else:
                        E_tn_domain = np.floor( poly_der(maxdist) *domain_size )
                        E_tn_gene = np.floor( poly_der(maxdist) *gene_size *0.8 )
                    
                    # Add zeros to the reads if the expected tn count is higher than observed
                    if E_tn_domain > len(reads):
                        
                        tn_add = np.zeros( int( E_tn_domain-len(reads) ) )
                        reads = np.append(reads,tn_add)
                    
                    # Store gene name when the number of tn is below the threshold
                    if (len(reads) < 2 and E_tn_gene>=5) or (sum(reads)==0 and E_tn_gene>=5): # domains with less than 2 insertions in genes with more than 5 insertions have zero fitness                 
                        fitness_perdomain.append(0)
                        gene_name.append( gene_info.name )
                        sem.append(0)
                    elif (len(reads) < 2 and E_tn_gene<5) or (sum(reads)==0 and E_tn_gene<5): # domains with less than 2 insertions in genes with less than 5 insertions are discarded                 
                        fitness_perdomain.append(np.nan)
                        gene_name.append( gene_info.name )
                        sem.append(np.nan)
                
                    else:
            
                        fitness_perdomain.append(np.log2(np.nanmean(reads)))
                        gene_name.append( gene_info.name )
                        sem.append(np.sqrt(np.std(reads)/len(reads)))

            fitness_domain[g]["fitness_domains"]=fitness_perdomain
            fitness_domain[g]["sem"]=sem
        fitness_domain_all.append(fitness_domain)    
    # -

    fitness_domain_all_pd=pd.concat([pd.DataFrame.from_dict(fitness_domain_all[i],orient="index") for i in np.arange(0,len(fitness_domain_all))],axis=0)
    return fitness_domain_all_pd


def compute_fitness_diff(fitness_average,fitness_domains):
    
    ## Compute the differences in all backgrounds
    fitness_corrected_pd=pd.DataFrame(index=fitness_average.index,columns=["fitness_domain_corrected","sem"])


    x_d=fitness_domains.index
    x_m=fitness_average.index
    common_indexes=list(set(x_d).intersection(set(x_m)))
    for j in common_indexes:
        vector = fitness_domains.loc[j,"fitness_domains"]
        number = fitness_average.loc[j,"mean_norm"]
        if vector.size==1 or vector.size==0:
            if vector<number:
                fitness_corrected_pd.loc[j,"fitness_domain_corrected"]=vector
                fitness_corrected_pd.loc[j,"sem"]=fitness_domains.loc[j,"sem"]
            else:
                fitness_corrected_pd.loc[j,"fitness_domain_corrected"]=number
                fitness_corrected_pd.loc[j,"sem"]=fitness_average.loc[j,"sem"]
        else:
            index = max_difference_index(vector, number)
            vector_max = vector[index]
            fitness_corrected_pd.loc[j,"fitness_domain_corrected"]=vector_max
            fitness_corrected_pd.loc[j,"sem"]=fitness_domains.loc[j,"sem"][index]

    return fitness_corrected_pd 
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 14:38:28 2023

@author: ekingma1
"""
import scipy
from scipy import stats
import numpy as np
import pandas as pd
import data_import_tools as read_data
from cen_bias import profile,correct 


maxdist = 200000
# Import data
per_gene_file="WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene_insertions.txt"
wiggle_file="WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_clean.wig_roman.wig"
gene_inserts = read_data.per_gene(per_gene_file)
wig_inserts = read_data.wiggle(wiggle_file)

# Correct for the centromere bias in the insertion density
chr_file="chromosome_sizes.txt"
cen_file="centromere_locations.txt"
c_dis, tn_csum, cen = profile(wig_inserts,chr_file=chr_file,cen_file=cen_file)
poly_fit = correct(c_dis, tn_csum,method="poly",maxdist=maxdist)
poly_der = np.polyder(poly_fit, 1)

# Estimate the expected number of insertions for each gene
exp_tn = []
reads_av = []
gene_name = []

for g in gene_inserts.index:
    
    # Get the info of the current gene
    gene_info = gene_inserts.loc[g]
    gene_size = gene_info['End location'] - gene_info['Start location']
    
    # Exclude insertions in the first and last 10% of the coding region of a gene
    ##
    # tn_pos = gene_info['Insertion locations']
    # tn_sel = ( tn_pos - gene_info['Start location'] > gene_size*0.1 ) \
    #         & ( tn_pos - gene_info['Start location'] < gene_size*0.9 )
     
    # reads = gene_info['Reads per insertion location'][tn_sel]
    ##
    reads_b =  ( wig_inserts["Chrom"] == "chr"+gene_info["Chromosome"] ) & \
               ( wig_inserts["Position"] > gene_info['Start location']+gene_size*0.1 ) & \
               ( wig_inserts["Position"] < gene_info['Start location']+gene_size*0.9 )
    
    reads = wig_inserts["Reads"][reads_b]
    
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
            continue
        
        # Store the average read count for the mean fitness
        reads_av.append( np.mean( reads_filt ) )
        gene_name.append( gene_info.name )
            
# Construct fitness array including confidence intervals
data = np.empty( ( len(reads_av),1 ) )
data[:,0] = [np.log2(x) for x in reads_av]
fitness = pd.DataFrame(data,columns=['mean'],index=gene_name)





# # Calculate shared dispersion
# disp = confidence_level.estimate_shared_disp(reads=y, mu=x)

# # Determine confidence intervals
# eps = confidence_level.conf_int(read_list=reads_per_tn_non_ess,dispersion=disp.params[0],CL=0.05)

# # Construct fitness array including confidence intervals
# data = np.empty((len(reads_av),4))
# data[:,0] = [np.log2(x) for x in reads_av]
# data[:,1] = [np.log2(x-e) for x, e in zip(reads_av,eps)]
# data[:,2] = [np.log2(x+e) for x, e in zip(reads_av,eps)]
# data[:,3] = [ x for x in  tn_count ]
# fitness = pd.DataFrame(data,columns=['mean','lower_CI','upper_CI','n_tn'],index=gene_names)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:03:05 2023

@author: linigodelacruz
"""

import scipy
from scipy import stats
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import ast
from collections import defaultdict
from from_excel_to_list import from_excel_to_list
from enzo_functions.data_import_tools import per_gene, wiggle
from enzo_functions.cen_bias import profile,correct 

# Correct for the centromere bias in the insertion density

maxdist = 200000
chr_file="../data/chromosome_sizes.txt"
cen_file="../data/centromere_locations.txt"

fitness_all=[]

## Importing pergene files 

pergene_files=[]
#data_dir= "../satay/data_files/data_unmerged/"
#data_dir="../transposonmapper/data_files/files4test/"
data_dir="../data/"
#data_dir="../transposonmapper/data_files/"
for root, dirs, files in os.walk(data_dir):
    for file in files:
        if file.endswith("_pergene_insertions.txt"):
            pergene_files.append(os.path.join(root, file))

wig_files=[]
for root, dirs, files in os.walk(data_dir):
    for file in files:
        if file.endswith("_clean.wig_roman.wig"):
            wig_files.append(os.path.join(root, file))
            
# Extract infor from per gene nd wig files 

gene_inserts=[]
wig_inserts=[]
for i in range(len(pergene_files)):
    gene_inserts.append(per_gene(pergene_files[i]))
    wig_inserts.append(wiggle(wig_files[i]))
    
keys=[]
for i in np.arange(len(pergene_files)):
    keys.append(pergene_files[i].split("/")[2])
    
# Compute the profile and correct it 

c_dis, tn_csum, cen = profile(wig_inserts[8],chr_file=chr_file,cen_file=cen_file)
poly_fit = correct(c_dis, tn_csum,method="poly",maxdist=maxdist)
poly_der = np.polyder(poly_fit, 1)

# Estimate the expected number of insertions for each gene
exp_tn = []
reads_av = []
gene_name = []
sem=[]
file=8
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
    
   
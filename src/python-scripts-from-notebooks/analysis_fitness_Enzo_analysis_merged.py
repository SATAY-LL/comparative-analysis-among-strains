# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: transposonmapper
#     language: python
#     name: python3
# ---

# +
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 14:38:28 2023

@author: ekingma1/linigodelacruz
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
from functions_analysis_frompergene2fitness import genes_discarded4fitness_from_flanking_regions



# -

# Ensure the entire structure is converted to floats
def convert_to_float(item):
    if isinstance(item, list):
        return [convert_to_float(elem) for elem in item]
    else:
        return float(item)


# +
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

import pickle
with open("../postprocessed-data/discarded_genes_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

# +
# Import data

gene_inserts=[]
wig_inserts=[]
for i in range(len(pergene_files)):
    gene_inserts.append(per_gene(pergene_files[i]))
    wig_inserts.append(wiggle(wig_files[i]))


# -

keys=[]
for i in np.arange(len(pergene_files)):
    keys.append(pergene_files[i].split("/")[2].split("\\")[0])


# +
## Big foor loop to compute the genes in with poor flanking regions in each background

# gene_inserts_all=pd.concat(gene_inserts,axis=0,keys=keys)

# for i in keys:
#     # add an extra column that is the index of the dataframe 
#     gene_inserts_all.loc[i,"Gene name"]=gene_inserts_all.loc[i].index

# # getting all discarded genes per background
# i=0
# discarded_genes_all=[]
# for background in keys:
    

#     flanking_regions,discarded_genes=genes_discarded4fitness_from_flanking_regions(gene_inserts_all,background,wig_files[i],windows_size=20000)
    
#     i=i+1
#     discarded_genes_all.append(discarded_genes)

# import pickle

# with open("../postprocessed-data/discarded_genes_all_backgrounds_enzo_keys", "wb") as fp:   #Pickling
#     pickle.dump(discarded_genes_all, fp)


# +
## BIG for loop across all backgrounds 

# +
# Correct for the centromere bias in the insertion density

maxdist = 200000
chr_file="../data/chromosome_sizes.txt"
cen_file="../data/centromere_locations.txt"
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

fitness_all_pd=pd.concat(fitness_all,axis=0,keys=keys)

# +
## Assign fitness of genes with poor flanking regions to NaN 

# for i in np.arange(0,len(keys)):
#     for j in fitness_all_pd.loc[keys[i]].index:
#         if j in b[i]:
#             fitness_all_pd.loc[keys[i],"mean"][j]=np.nan
#             fitness_all_pd.loc[keys[i],"sem"][j]=np.nan




##Saving
# fitness_discarding_genes_pd.to_csv("../postprocessed-data/fitness_all_Enzo_analysis.csv")

# with open("../postprocessed-data/fitness_models_all_backgrounds_enzo_analysis_merged", "wb") as fp:   #Pickling
#     pickle.dump(fitness_discarding_genes_pd, fp)
# -

with open("../postprocessed-data/fitness_models_all_backgrounds_enzo_analysis_merged", "wb") as fp:   #Pickling
     pickle.dump(fitness_all_pd, fp)

# +
# with open("../postprocessed-data/fitness_models_all_backgrounds_enzo_analysis_merged", "rb") as fp:   # Unpickling
#     fitness_average_all = pickle.load(fp)
# -

# ### Fitness for domains 

# +
data_domains=pd.read_excel("../postprocessed-data/genomic-domains-wt.xlsx",index_col="Unnamed: 0")

## data from yeastmine
domains_names=pd.read_csv('../data/Domains_all_genes_protein_coordinates_yeastmine.tsv',sep="\t")
domains_names.index=domains_names["Gene Name"]

# +
# From the loop per gene select the domains per gene from the data domains file and convert them to floats. 
# ref=np.nanmedian(np.log2(wig_inserts[0]["Reads"]))

fitness_domain=defaultdict(dict)
fitness_domain_all=[]
cen_file="../data/centromere_locations.txt"
chr_file="../data/chromosome_sizes.txt"

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

fitness_domain_all_pd=pd.concat([pd.DataFrame.from_dict(fitness_domain_all[i],orient="index") for i in np.arange(0,len(fitness_domain_all))],axis=0,keys=keys)

# +
## Assign fitness of genes with poor flanking regions to NaN 

# for i in np.arange(0,len(keys)):
#     for j in fitness_domain_all_pd.loc[keys[i]].index:
#         if j in b[i]:
#             fitness_domain_all_pd.loc[keys[i],"fitness_domains"][j]=np.nan
#             fitness_domain_all_pd.loc[keys[i],"sem"][j]=np.nan

# +
# replace the nan values with zero of the fitness per domains (optional but distorts the distribution)

# for i in np.arange(0,len(keys)):
#     for j in fitness_domain_all_pd.loc[keys[i]].index:
#         for k in np.arange(0,len(fitness_domain_all_pd.loc[keys[i]].loc[j,"fitness_domains"])):
#             if np.isnan(fitness_domain_all_pd.loc[keys[i]].loc[j,"fitness_domains"][k]):
#                 fitness_domain_all_pd.loc[keys[i]].loc[j,"fitness_domains"][k]=0
# -

# to compute the median of the fitness of each domain across all backgrounds
for i in np.arange(0,len(keys)):
    
    x=fitness_domain_all_pd.loc[keys[i],"fitness_domains"].tolist()
    flattened_list = [item for sublist in x if isinstance(sublist, list) for item in sublist] + [item for item in x if isinstance(item, (int, float))]
    ref=np.nanmedian(flattened_list)
    # divide each fitness value by the reference value
    for j in fitness_domain_all_pd.loc[keys[i]].index:
        fitness_domain_all_pd.loc[keys[i]].loc[j,"fitness_domains"]=fitness_domain_all_pd.loc[keys[i]].loc[j,"fitness_domains"]/ref

with open("../postprocessed-data/fitness_models_all_backgrounds_domain_analysis_merged", "wb") as fp:   #Pickling
     pickle.dump(fitness_domain_all_pd, fp)


# +

def max_difference_index(vector, number):
    # Calculate the element-wise absolute differences
    differences = np.abs(vector - number)
    
    # Find the index where the difference is maximum
    max_index = np.argmax(differences)
    
    return max_index



# -

with open("../postprocessed-data/fitness_models_all_backgrounds_enzo_analysis_merged", "rb") as fp:   # Unpickling
    fitness_all_pd = pickle.load(fp)

# +
## Compute the differences in all backgrounds
fitness_corrected_pd=pd.DataFrame(index=fitness_all_pd.index,columns=["fitness_domain_corrected","sem"])

for i in keys:
    x_d=fitness_domain_all_pd.loc[i].index
    x_m=fitness_all_pd.loc[i].index
    common_indexes=list(set(x_d).intersection(set(x_m)))
    for j in common_indexes:
        vector = fitness_domain_all_pd.loc[i,"fitness_domains"].loc[j]
        number = fitness_all_pd.loc[i,"mean_norm"].loc[j]
        if vector.size==1 or vector.size==0:
            if vector<number:
                fitness_corrected_pd.loc[i,"fitness_domain_corrected"].loc[j]=vector
                fitness_corrected_pd.loc[i,"sem"].loc[j]=fitness_domain_all_pd.loc[i,"sem"].loc[j]
            else:
                fitness_corrected_pd.loc[i,"fitness_domain_corrected"].loc[j]=number
                fitness_corrected_pd.loc[i,"sem"].loc[j]=fitness_all_pd.loc[i,"sem"].loc[j]
        else:
            index = max_difference_index(vector, number)
            vector_max = vector[index]
            fitness_corrected_pd.loc[i,"fitness_domain_corrected"].loc[j]=vector_max
            fitness_corrected_pd.loc[i,"sem"].loc[j]=fitness_domain_all_pd.loc[i,"sem"].loc[j][index]
# -

# replace -inf values of te fitness with zero 
for i in keys: 
    fitness_corrected_pd.loc[i].replace(-np.inf,0,inplace=True)

with open("../postprocessed-data/fitness_models_all_backgrounds_domain_analysis_corrected", "wb") as fp:   #Pickling
     pickle.dump(fitness_corrected_pd, fp)

fitness_corrected_pd

fitness_corrected_pd.loc["wt_merged"].fitness_domain_corrected.hist(bins=100)
fitness_all_pd.loc["wt_merged"].mean_norm.hist(bins=100)

# +
## Plots for the figure of the average transposon density based on the missing sites 

# +
maxdist = 200000
chr_file="../data/chromosome_sizes.txt"
cen_file="../data/centromere_locations.txt"
# Compute the profile and correct it 
file=8 # wt_merged
c_dis, tn_csum, cen = profile(wig_inserts[file],chr_file=chr_file,cen_file=cen_file)
poly_fit = correct(c_dis, tn_csum,method="poly",maxdist=1000000)
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
    
   

# +
 ## Linear fit

for i in np.arange(0,tn_csum.shape[0]):
    plt.plot(c_dis,tn_csum[i],color="grey",alpha=0.3)

mean_tn_csum=np.nanmean(tn_csum,axis=0)
plt.plot(c_dis,mean_tn_csum,color="blue",linewidth=1,label="Average")

x=np.arange(0,400000,10000)
x2fit=np.arange(200000,400000,10000)
model = np.polyfit(x2fit, mean_tn_csum[20:40], 1)
plt.plot(x,model[0]*x+model[1],color="red",label="Linear fit",alpha=0.8)
#plt.plot(model[0]*x+model[1],color="blue",label="Linear fit")
plt.text(200000, 5000, 'y=%.3fx+%.1f' % (model[0], model[1]), fontsize=10)
xticks_values = [100000, 200000, 300000, 400000]  # Specify the x-axis tick values
xticks_labels = np.round(np.divide(xticks_values,1000),0)  # Specify the labels for the tick values

# Set the x-axis tick locations and labels
plt.xticks(xticks_values, xticks_labels)
plt.xlim(0,400000)
plt.ylim(0,25000)
plt.legend()

#plt.savefig("../figures/centromere_bias_linear_fit.png",dpi=300)

# +
## average cumulative insertions vs distance to centromere 

mean_tn_csum=np.nanmean(tn_csum,axis=0)
plt.plot(c_dis,mean_tn_csum,color="blue",linewidth=1,label="Average")
plt.plot(c_dis,poly_fit(c_dis),color="red",linewidth=1,label="Polynomial fit")
plt.xlim(0,400000)
plt.ylim(0,25000)

plt.xlabel("Distance to centromere (bp)")
plt.ylabel("Cumulative insertions")
plt.legend()
#plt.savefig("../figures/centromere_bias_polynomial_fit.png",dpi=300)
# -

poly_fit

# +
## Plot insertion rate vs distance to centromere and polynomial fit

plt.plot(c_dis[0:-1],np.diff(mean_tn_csum)/10000,color="blue",linewidth=1,label="Average insertion rate")
plt.plot(c_dis,poly_der(c_dis),color="red",linewidth=1,label="Polynomial fit")
plt.xlim(0,400000)
plt.ylim(0,0.1)
plt.ylabel("Insertion rate ($bp^{-1}$)")
plt.xlabel("Distance to centromere (bp)")

plt.legend()

#plt.savefig("../figures/centromere_bias_polynomial_fit_insertion_rate.png",dpi=300)
# -

poly_der

# +
## Expected insertion density minus the observed insertion density over essential and not essential genes

standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str) 
diff_essentials=[]
diff=[]

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
            
            # Determine the expected insertion rate
            if d_cen <= maxdist:           
                E_tn = np.floor( poly_der(d_cen)*gene_size*0.8 ) 
                
            else:
                E_tn = np.floor( poly_der(maxdist)*gene_size*0.8  )
            if g in standard_essentials:
            # compute the differences from the expected insertion density and the observed insertion density
                diff_essentials.append(E_tn-len(reads))
            else:
                diff.append(E_tn-len(reads))

# +
plt.hist(diff,bins=600,histtype="stepfilled",color="gray",label="Non essential");
plt.hist(diff_essentials,bins=600,histtype="stepfilled",color="pink",label="Essential");

plt.xlabel("E(x)-O(x)")
plt.ylabel("Count")

plt.legend()
plt.tight_layout()
plt.savefig("../figures/centromere_bias_expected_observed_insertions.png",dpi=300)

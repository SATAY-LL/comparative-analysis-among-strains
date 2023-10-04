# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
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

from data_import_tools import per_gene, wiggle
from cen_bias import profile,correct 




# -

# Ensure the entire structure is converted to floats
def convert_to_float(item):
    if isinstance(item, list):
        return [convert_to_float(elem) for elem in item]
    else:
        return float(item)


# +
## Importing pergene files 

pergene_files=["WT_merged_pergene_insertions.txt"]

wig_files=["WT_merged.wig"]

# +
# Import data

gene_inserts=[]
wig_inserts=[]
for i in range(len(pergene_files)):
    gene_inserts.append(per_gene(pergene_files[i]))
    wig_inserts.append(wiggle(wig_files[i]))



# +
# Correct for the centromere bias in the insertion density

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
# -

# ### Fitness for domains 

# +
data_domains=pd.read_excel("genomic-domains-wt.xlsx",index_col="Unnamed: 0")

## data from yeastmine
domains_names=pd.read_csv('Domains_all_genes_protein_coordinates_yeastmine.tsv',sep="\t")
domains_names.index=domains_names["Gene Name"]

# +
# From the loop per gene select the domains per gene from the data domains file and convert them to floats. 
# ref=np.nanmedian(np.log2(wig_inserts[0]["Reads"]))

fitness_domain=defaultdict(dict)
fitness_domain_all=[]
cen_file="centromere_locations.txt"
chr_file="data/chromosome_sizes.txt"

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


# to compute the median of the fitness of each domain

    
x=fitness_domain_all_pd.loc["fitness_domains"].tolist()
flattened_list = [item for sublist in x if isinstance(sublist, list) for item in sublist] + [item for item in x if isinstance(item, (int, float))]
ref=np.nanmedian(flattened_list)
# divide each fitness value by the reference value
for j in fitness_domain_all_pd.index:
    fitness_domain_all_pd.loc.loc[j,"fitness_domains"]=fitness_domain_all_pd.loc[j,"fitness_domains"]/ref




# +

def max_difference_index(vector, number):
    # Calculate the element-wise absolute differences
    differences = np.abs(vector - number)
    
    # Find the index where the difference is maximum
    max_index = np.argmax(differences)
    
    return max_index





# +
## Compute the differences in all backgrounds
fitness_corrected_pd=pd.DataFrame(index=fitness_all_pd.index,columns=["fitness_domain_corrected","sem"])


x_d=fitness_domain_all_pd.index
x_m=fitness_all_pd.index
common_indexes=list(set(x_d).intersection(set(x_m)))
for j in common_indexes:
    vector = fitness_domain_all_pd.loc["fitness_domains"].loc[j]
    number = fitness_all_pd.loc["mean_norm"].loc[j]
    if vector.size==1 or vector.size==0:
        if vector<number:
            fitness_corrected_pd.loc["fitness_domain_corrected"].loc[j]=vector
            fitness_corrected_pd.loc["sem"].loc[j]=fitness_domain_all_pd.loc["sem"].loc[j]
        else:
            fitness_corrected_pd.loc["fitness_domain_corrected"].loc[j]=number
            fitness_corrected_pd.loc["sem"].loc[j]=fitness_all_pd.loc["sem"].loc[j]
    else:
        index = max_difference_index(vector, number)
        vector_max = vector[index]
        fitness_corrected_pd.loc["fitness_domain_corrected"].loc[j]=vector_max
        fitness_corrected_pd.loc["sem"].loc[j]=fitness_domain_all_pd.loc["sem"].loc[j][index]
# -






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

plt.savefig("../figures/centromere_bias_linear_fit.png",dpi=300)

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
plt.savefig("../figures/centromere_bias_polynomial_fit.png",dpi=300)

# +
## Plot insertion rate vs distance to centromere and polynomial fit

plt.plot(c_dis[0:-1],np.diff(mean_tn_csum)/10000,color="blue",linewidth=1,label="Average insertion rate")
plt.plot(c_dis,poly_der(c_dis),color="red",linewidth=1,label="Polynomial fit")
plt.xlim(0,400000)
plt.ylim(0,0.1)
plt.ylabel("Insertion rate ($bp^{-1}$)")
plt.xlabel("Distance to centromere (bp)")

plt.legend()



# +
## Expected insertion density minus the observed insertion density over essential and not essential genes

standard_essentials=np.loadtxt("standard_essentials.txt",dtype=str) 
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


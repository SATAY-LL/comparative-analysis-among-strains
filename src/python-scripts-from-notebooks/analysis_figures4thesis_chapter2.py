# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3.9.7 ('transposonmapper')
#     language: python
#     name: python3
# ---

# ## Transposon insertion data
# ### Describe our first dataset on WT. 
# - Show the differences in read count and transposon counts in technical replicates( they are replicates splited before the PCR). Cite Enzo thesis to say that PCR amplification is one of the major noise sources in the read counts of satay libraries, to explain the differences across replicates.
# - Show the known biases of this type of transposon system to genes located near to centromeres. 
# - Plot the cumulative plot of transposons against the distance of every gene to the centromere of its chromosome. 
# - Show the tendency of transposons to insert in non coding regions and that this does not change  across different genetic backgrounds (wt,dnrp1,dbem3,dbem1dbem3) 
#

# +
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys
from collections import defaultdict
from ast import literal_eval
from scipy.stats import norm

from from_excel_to_list import from_excel_to_list

# +
## Plot transposons along the whole genome


## Importing pergene files 

pergene_files=[]
#data_dir= "../satay/data_files/data_unmerged/"
#data_dir="../transposonmapper/data_files/files4test/"
data_dir="../postprocessed-data/"
#data_dir="../transposonmapper/data_files/"
for root, dirs, files in os.walk(data_dir):
    for file in files:
        if file.endswith("pergene_insertions.xlsx"):
            pergene_files.append(os.path.join(root, file))

list_data=[]
for i in pergene_files:
    list_data.append(pd.read_excel(i,engine='openpyxl',index_col="Unnamed: 0"))

keys=[]
for i in np.arange(0,len(pergene_files)):
    keys.append(pergene_files[i].split("/")[-1].split("_")[0]+"_"+pergene_files[i].split("/")[-1].split("_")[1])

list_data_pd=pd.concat(list_data,axis=0,keys=keys)

# -

keys= ['wt_merged','wt_a','wt_b','dnrp1_1','dnrp1_2']
data=list_data_pd.loc[keys] # Take only data from targeted genotypes

# +
from sklearn.linear_model import LinearRegression

#initiate linear regression model
model = LinearRegression()

#define predictor and response variables



figure,ax=plt.subplots(1,2,figsize=(10,5))
plt.subplots_adjust(wspace=0.4,hspace=0.5)

magnitudes=["Reads","Insertions"]
strains=["wt_a","wt_b","dnrp1_1","dnrp1_2"]

for i in np.arange(0,len(magnitudes)):
    tmp_0=data.loc[strains[0],magnitudes[i]]/data.loc[strains[0],magnitudes[i]].sum()
    tmp_1=data.loc[strains[1],magnitudes[i]]/data.loc[strains[1],magnitudes[i]].sum()

    X, y = tmp_0.values.reshape(-1,1), tmp_1.values.reshape(-1,1)

    #fit regression model
    model.fit(X, y)

    #calculate R-squared of regression model
    r_squared = model.score(X, y)


    ax[i].scatter(tmp_0,tmp_1,s=20,alpha=0.5)
    ax[i].set_title("Normalized "+magnitudes[i])
    ax[i].set_xlabel("tech replicate 1")
    ax[i].set_ylabel("tech replicate 2")
    ax[i].set_xlim(0,0.001)
    ax[i].set_ylim(0,0.001)
    ax[i].plot([0,0.001],[0,0.001],color="black",linestyle="--")
    ax[i].text(0, 0.0005, '$R^2=%.3f$' % (r_squared),fontsize=12)

    # ax[i,1].set_title( magnitudes[i] + " difference")
    # ax[i,1].hist(tmp_1-tmp_0,bins=100)
    # data2fit = tmp_1-tmp_0
    # mu, std = norm.fit(data2fit)
    
    # xmin, xmax = ax[i,1].get_xlim()
    # x = np.linspace(xmin, xmax, 100)
    # p = norm.pdf(x, mu, std)

    # ax[i,1].plot(x, p, 'k', linewidth=2)

   
figure,ax=plt.subplots(1,2,figsize=(10,5))
plt.subplots_adjust(wspace=0.4,hspace=0.5)

for i in np.arange(0,len(magnitudes)):
    tmp_0=data.loc[strains[2],magnitudes[i]]/data.loc[strains[2],magnitudes[i]].sum()
    tmp_1=data.loc[strains[3],magnitudes[i]]/data.loc[strains[3],magnitudes[i]].sum()

    X, y = tmp_0.values.reshape(-1,1), tmp_1.values.reshape(-1,1)

    #fit regression model
    model.fit(X, y)

    #calculate R-squared of regression model
    r_squared = model.score(X, y)


    ax[i].scatter(tmp_0,tmp_1,s=20,alpha=0.5)
    ax[i].set_title("Normalized "+magnitudes[i])
    ax[i].set_xlabel("biolog replicate 1")
    ax[i].set_ylabel("biolog replicate 2")
    ax[i].set_xlim(0,0.001)
    ax[i].set_ylim(0,0.001)
    ax[i].plot([0,0.001],[0,0.001],color="black",linestyle="--")
    ax[i].text(0, 0.0005, '$R^2=%.3f$' % (r_squared),fontsize=12)


#figure.savefig("../figures/fig_differences_WT_technical_replicates.png",dpi=300)


# -

# ## Using read counts for fitness 
# - Simplified fitness model (simple malthusian model where the differences in mutant abundance is proportional to its growth rate, and there is no limiting factor)
#
# For the fitness calculation:
#
# - Take only the 80% central part of the gene (show the overall enrichment towards the first and last part in average across all genes.)
# - Annotate which genes "we can not analyse" because they dont have enough data. So, the challenge here is to define what is a region where we dont have sufficient data to say something about it. 
#     - use the density of the flanking regions as the reference for the gene transposon density.
#
# - For the "analysable genes" compute the fitness using the maltusian model of the:
#     - whole gene by averaging the read counts over the number of transposons
#     - subdivide the gene into known domains and compute the fitness of individual domains
#     - assign a new corrected fitness to the every gene by taking the fitness of the domain with strongest effect with respect to the average , so the max|f_mean-f_domain_i|
#  
#  Then :
#  - Plot the distributions of every fitness calculation referenced to the median fitness value of the population.
#  - Plot the average fitness vs the corrected fitness 
#
#
#

data=list_data_pd.loc["wt_merged"]
data.head()

# ## Computing the number of reads per insertion along the gene length for all genes 

# +
reads_locations_center=[]
reads_locations=[]
insertion_locations=[]
insertion_locations_center=[]

gene_coordinates=[]
data=list_data_pd.loc["wt_merged"]
gene_len_data=[]

read_insert_pair=defaultdict(dict)
for j in data.index:
        coi=from_excel_to_list(data.loc[j]["Reads per insertion location"])
        coi_insertions=from_excel_to_list(data.loc[j]["Insertion locations"])
        gene_coordinates.append([data.loc[j,"Start location"],data.loc[j,"End location"]])
        reads_locations.append(coi)
        insertion_locations.append(coi_insertions)
        read_insert_pair[j]=(coi,coi_insertions)

        gene_len=data.loc[j,"End location"]-data.loc[j,"Start location"]
        gene_len_data.append(gene_len)
        if type(coi)!=int: # for genes that have more than one read
            start=gene_len*0.1
            end=gene_len*0.9
            data_center=coi[int(start):int(end)]
            data_center_insertions=coi_insertions[int(start):int(end)]
            reads_locations_center.append(data_center)
            insertion_locations_center.append(data_center_insertions)
        else:
            data_center=coi
            reads_locations_center.append(coi)
            insertion_locations_center.append(coi_insertions)

# +
# compute the reads per insertions for each gene along the gene length

gene_parts=np.linspace(0,1,21)
r=np.zeros(shape=(len(insertion_locations),len(gene_parts))) # reads_per_insertion_parts array , every gene in the rows 

for i in np.arange(0,len(insertion_locations)):
    if (insertion_locations[i])!=0:
        g=np.array(insertion_locations[i])
        f=np.linspace(gene_coordinates[i][0],gene_coordinates[i][1],len(gene_parts))
        #f=np.array(gene_coordinates[i][0]+gene_coordinates[i][1]*gene_parts)
        binedges = g.searchsorted(f)

        rngs = [list(range(binedges[k], binedges[k+1])) for k in range(len(binedges)-1)]

        for k in np.arange(0,len(rngs)):
            readsperinsert=[]
            for j in np.arange(0,len(rngs[k])):
                readsperinsert.append(reads_locations[i][rngs[k][j]])
            r[i,k]=np.sum(readsperinsert)/len(readsperinsert)

    


# -

gene_parts,np.linspace(0,1,11)

# +
# replace nan values with zeros

r[np.isnan(r)]=0

# Sum all values of r along the rows

r_sum=np.sum(r,axis=0)

# Plot r values along gene parts

ax=plt.figure(figsize=(8,5))

# for i in np.arange(1,len(insertion_locations)):
#     plt.plot(gene_parts,r[i,:],color="black",alpha=0.3,marker="o",markersize=1)

# plt.plot(gene_parts,r[0,:],color="black",alpha=0.3,marker="o",markersize=1,
#     label="Data per gene")
# plt.yscale("log")
plt.xlabel("Gene length")
plt.ylabel("Reads per insertion per gene")

plt.bar(gene_parts,r_sum,color="black",alpha=0.3,width=gene_parts[-1]/len(gene_parts),label="Sum of all genes")
plt.xticks(np.linspace(0,1,11),labels=["0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"]);

plt.legend()
plt.tight_layout()

#plt.savefig("../figures/fig_reads_per_insertion_all_genes_along_gene_length.png",dpi=300,format="png")
## from the figure it follows that I should remove only the first 10% of the reads to compute the fitness of the whole gene


# +

keys= ['bem1-aid_a','dbem1dbem3_b','wt_merged','dbem1dbem3_a', 
'dnrp1_merged','bem1-aid_b','dbem3_merged']


# -

# ## Compute poor enriched flanking regions for each gene, to classify genes as with poor data or not , for the fitness 
# ### Procedure:
# 1. Read the wig file and annotate per genomic location how many reads are found.
# 2. Go gene by gene start and end location and count how many insertions are 10kb upstream and downstream
# 3. Compute the average insertion density per chromosome being the total number of insertions per chromosome over its length. 
# 4. For every gene xompute the number of insertions 10kb upstream and downstream using the wig file. Discard the genes that the sum of insertions in either the flanking regions is less than the expected number assuming it is 10kb*chromosome insertion density. 
#  
# note: the centromere bias is neglected for now, so we dont expect to find low enrich areas close to centromeres (due to bias) .
#

# +
import wiggelen as wg

## Importing the required python libraries and using the transposon_bias2centromeres function

wigfile="../data/wt_merged/WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_clean.wig"
centromeres_file="../data/centromeres_genomic_locations.csv"

## Import wig files 
wig_reads=[]

for x in wg.walk(open(wigfile)):
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

# +
# size of every chromosome in bp from https://www-ncbi-nlm-nih-gov.tudelft.idm.oclc.org/genome/?term=Saccharomyces%20cerevisiae%5BOrganism%5D&cmd=DetailsSearch
chromosome_size={"I":230218,"II":813184,"III":316620,
"IV":1531933,"V":576874,"VI":270161,"VII":1090940,
"VIII":562643,"IX":439888,"X":745751,"XI":666816,
"XII":1078177,"XIII":924431,"XIV":784333,"XV":1091291,
"XVI":948066}

chrom2latin={"I":1,"II":2,"III":3,"IV":4,"V":5,"VI":6,"VII":7,"VIII":8,"IX":9,
"X":10,"XI":11,"XII":12,"XIII":13,"XIV":14,"XV":15,"XVI":16}

cumulative_bp_along_chrom=np.cumsum(list(chromosome_size.values()))

# +
## densities over each chromosome
chrom_conversion={"chrI":"I","chrII":"II","chrIII":"III",
"chrIV":"IV","chrV":"V","chrVI":"VI","chrVII":"VII","chrVIII":"VIII",
"chrIX":"IX","chrX":"X","chrXI":"XI","chrXII":"XII","chrXIII":"XIII",
"chrXIV":"XIV","chrXV":"XV","chrXVI":"XVI"}

chromosome_density=[]

for i in chrom_conversion.keys():
    tmp=len(pos[chrom==i])/chromosome_size[chrom_conversion[i]]
    chromosome_density.append(tmp)

expected_flanking_enrichment=np.round(np.multiply(chromosome_density,10000),0).astype(int)

expected_flanking_enrichment

# +

data=list_data_pd.loc["wt_merged"]
windows_size=10000

flanking_regions_data=defaultdict(dict)
for gene in data["Gene name"]:

    data_loc_gene=data[data.loc[:,"Gene name"]==gene]
    data_loc=[data_loc_gene["End location"],data_loc_gene["Start location"]]

    gene_chromosome=data_loc_gene["Chromosome"].tolist()[0]

    if gene_chromosome!="Mito":

        
        # flanking regions
        flanking_regions=[data_loc[0].tolist()[0]-windows_size,data_loc[1].tolist()[0]+windows_size]
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
        if insertions_flanking<expected_flanking_enrichment[chrom2latin[gene_chromosome]-1]:
            flanking_regions_data[gene]["classification"]="Not enough flanking regions"
        else:
            flanking_regions_data[gene]["classification"]="OK"



# +
flanking_regions_data_pd=pd.DataFrame.from_dict(flanking_regions_data,orient="index")

discarded_genes2fitness=flanking_regions_data_pd[flanking_regions_data_pd["classification"]=="Not enough flanking regions"].index

# +

from functions_satay_biases import transposon_bias2centromeres

fig, distance2cent_all = transposon_bias2centromeres(wigfile,centromeres_file,save=False)

# -

distance2cent_all_chrom=(list(distance2cent_all.values()))
figure=plt.subplots(4,4,figsize=(16,12))
plt.subplots_adjust(hspace=0.5,wspace=0.4)
# using colormap plasma 
color = plt.get_cmap('jet') 
j=1
for i in np.arange(0,len(distance2cent_all_chrom)):
    plt.subplot(4,4,j)
    #plt.hist(distance2cent_all_chrom[i],bins=50,alpha=0.4,label=list(distance2cent_all.keys())[i],density=True,color=color(i/len(distance2cent_all_chrom)));
    #sns.histplot(distance2cent_all_chrom[i],bins=50,alpha=0.4,kde=True,label=list(distance2cent_all.keys())[i],color=color(i/len(distance2cent_all_chrom)));
    sns.kdeplot(distance2cent_all_chrom[i],alpha=0.4,
    label=list(distance2cent_all.keys())[i],color=color(i/len(distance2cent_all_chrom)));

    plt.xlabel("Distance to centromere (bp)")
    plt.ylabel("Transposons insertions")
    #plt.ylim(0,3000)
    plt.xlim(-1000000,2500000)
    plt.legend()

    j=j+1
  

# +
## import protein domains

## data generated in the notebook "analysis_proteins_domains.ipynb"
data_domains=pd.read_excel("../postprocessed-data/genomic-domains-wt.xlsx",index_col="Unnamed: 0")
data_reads_domains=pd.read_excel("../postprocessed-data/reads-per-domain-all-backgrounds-new.xlsx",index_col="Unnamed: 0")
data_insertions_domains=pd.read_excel("../postprocessed-data/insertions-per-domain-all-backgrounds-new.xlsx",index_col="Unnamed: 0")


## data from yeastmine
domains_names=pd.read_csv('../data/Domains_all_genes_protein_coordinates_yeastmine.tsv',sep="\t")
domains_names.index=domains_names["Gene Name"]

# -

domains_names

# +
protein_domains_data=defaultdict(dict)

for i in data.index:

    gene=data.loc[i,"Gene name"]
    
    protein_domains_data[gene]["gene coordinates"]= np.round(np.linspace(gene_coordinates[i][0],
    gene_coordinates[i][1],10),0)
    central_gene_part=protein_domains_data[gene]["gene coordinates"][2:8]
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
        protein_domains_data[gene]["reads_domain"]= from_excel_to_list(data_reads_domains.loc[gene][0])
        protein_domains_data[gene]["insertions_domain"]= from_excel_to_list(data_insertions_domains.loc[gene][0])
    else:
        protein_domains_data[gene]["domains coordinates"]= np.nan
        protein_domains_data[gene]["domains name"]= np.nan
        protein_domains_data[gene]["domain description"]= np.nan
        protein_domains_data[gene]["reads_domain"]= np.sum(reads_locations[i])
        protein_domains_data[gene]["insertions_domain"]= insertion_locations[i]

    if gene in discarded_genes2fitness:
        protein_domains_data[gene]["classification"]="Not enough flanking regions"
    else:
        protein_domains_data[gene]["classification"]="OK"
        


# -

protein_domains_data_pd=pd.DataFrame.from_dict(protein_domains_data,orient="index")
protein_domains_data_pd.head(5)

# +
# how many genes have domains?

protein_domains_data_pd[protein_domains_data_pd["domains coordinates"].isna()==False].shape

# +
## fitness - malthusian model




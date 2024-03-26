# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
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

## standard for plots
plt.rc('font', family='serif',size=14)
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

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
    tmp=pd.read_excel(i,engine='openpyxl',index_col="Unnamed: 0")
    ## remove ADE2 genes
    tmp=tmp[tmp.loc[:,"Gene name"]!="ADE2"]
    tmp.index=np.arange(0,len(tmp))
    list_data.append(tmp)


keys=[]
for i in np.arange(0,len(pergene_files)):
    keys.append(pergene_files[i].split("/")[-1].split("_")[0]+"_"+pergene_files[i].split("/")[-1].split("_")[1])

list_data_pd=pd.concat(list_data,axis=0,keys=keys)

# -

standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str)

# keys= ['wt_merged','wt_a','wt_b','dnrp1_1','dnrp1_2']
# data=list_data_pd.loc[keys] # Take only data from targeted genotypes


for i in keys:
    print(i,"The number of reads are",list_data_pd.loc[i,"Reads"].sum())


for i in keys:
    print(i,"The number of insertions are",list_data_pd.loc[i,"Insertions"].sum())


# +
data=list_data_pd
data_wta=data.loc["wt_a"]
data_wtb=data.loc["wt_b"]

data_wta.index=data_wta["Gene name"]
data_wtb.index=data_wtb["Gene name"]

# are the index of data_wta and data_wtb the same?
print(data_wta.index.equals(data_wtb.index))




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


    ax[i].scatter(tmp_0,tmp_1,s=1,alpha=0.5,color="black")
    ax[i].set_title("Normalized "+magnitudes[i])
    ax[i].set_xlabel("tech replicate 1")
    ax[i].set_ylabel("tech replicate 2")
    ax[i].set_xlim(1e-7,1e-1)
    ax[i].set_ylim(1e-7,1e-1)
    ax[i].plot([1e-7,1e-1],[1e-7,1e-1],color="black",linestyle="--")
    ax[i].text(1e-4, 5e-3, '$R^2=%.3f$' % (r_squared),fontsize=12)

    ax[i].set_aspect('equal', 'box')
    ax[i].set_yscale('log')
    ax[i].set_xscale('log')
    # ax[i,1].set_title( magnitudes[i] + " difference")
    # ax[i,1].hist(tmp_1-tmp_0,bins=100)
    # data2fit = tmp_1-tmp_0
    # mu, std = norm.fit(data2fit)
    
    # xmin, xmax = ax[i,1].get_xlim()
    # x = np.linspace(xmin, xmax, 100)
    # p = norm.pdf(x, mu, std)

    # ax[i,1].plot(x, p, 'k', linewidth=2)

plt.tight_layout()
#figure.savefig("../figures/fig_differences_technical_replicates.png",dpi=300)
   
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


    ax[i].scatter(tmp_0,tmp_1,s=1,alpha=0.5,color="black")
    ax[i].set_title("Normalized "+magnitudes[i])
    ax[i].set_xlabel("biolog replicate 1")
    ax[i].set_ylabel("biolog replicate 2")
    ax[i].set_xlim(1e-7,1e-1)
    ax[i].set_ylim(1e-7,1e-1)
    
    ax[i].plot([1e-7,1e-1],[1e-7,1e-1],color="black",linestyle="--")
    ax[i].text(1e-4, 5e-3, '$R^2=%.3f$' % (r_squared),fontsize=12)
    ax[i].set_aspect('equal', 'box')
    ax[i].set_yscale('log')
    ax[i].set_xscale('log')


plt.tight_layout()
#figure.savefig("../figures/fig_differences_biological_replicates.png",dpi=300)



# +
from sklearn.linear_model import LinearRegression

#initiate linear regression model
model = LinearRegression()

#define predictor and response variables



figure,ax=plt.subplots(1,1,figsize=(10,5))
plt.subplots_adjust(wspace=0.4,hspace=0.5)

magnitudes=["Reads","Insertions"]
strains=["wt_a","wt_b","dnrp1_1","dnrp1_2"]

# compute the normalized reads per insertions per gene for each strain

tmp_0=(data.loc[strains[0],magnitudes[0]]*data.loc[strains[0],magnitudes[1]].sum())/(data.loc[strains[0],magnitudes[1]]*data.loc[strains[0],magnitudes[0]].sum())
tmp_1=(data.loc[strains[1],magnitudes[0]]*data.loc[strains[1],magnitudes[1]].sum())/(data.loc[strains[1],magnitudes[1]]*data.loc[strains[1],magnitudes[0]].sum())

## Remove Nan and infinite
tmp_0=tmp_0[np.isfinite(tmp_0)]
tmp_1=tmp_1[np.isfinite(tmp_1)]

## Make both variables same length
tmp_0=tmp_0[tmp_0.index.isin(tmp_1.index)]
tmp_1=tmp_1[tmp_1.index.isin(tmp_0.index)]

X, y = tmp_0.values.reshape(-1,1), tmp_1.values.reshape(-1,1)

#fit regression model
model.fit(X, y)

#calculate R-squared of regression model
r_squared = model.score(X, y)


ax.scatter(tmp_0,tmp_1,s=1,alpha=0.5,color="black")
ax.set_title("Normalized "+"Reads per insertions per gene")
ax.set_xlabel("tech replicate 1")
ax.set_ylabel("tech replicate 2")
ax.set_xlim(0,4)
ax.set_ylim(0,4)
ax.plot([0,4],[0,4],color="black",linestyle="--")
ax.text(0, 2, '$R^2=%.3f$' % (r_squared),fontsize=12)

ax.set_aspect('equal', 'box')
    
    # ax[i,1].set_title( magnitudes[i] + " difference")
    # ax[i,1].hist(tmp_1-tmp_0,bins=100)
    # data2fit = tmp_1-tmp_0
    # mu, std = norm.fit(data2fit)
    
    # xmin, xmax = ax[i,1].get_xlim()
    # x = np.linspace(xmin, xmax, 100)
    # p = norm.pdf(x, mu, std)

    # ax[i,1].plot(x, p, 'k', linewidth=2)

plt.tight_layout()
#figure.savefig("../figures/fig_differences_technical_replicates.png",dpi=300)



# -

# ### Comparing libraries with the wig files , therefore the insertions along all genome 

# +
import wiggelen as wg

## Importing the required python libraries and using the transposon_bias2centromeres function

wigfiles={"wt_a":"../data/wt_a/WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_clean.wig",
"wt_b":"../data/wt_b/WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_clean.wig",
"dnrp1_1":"../data/dnrp1_1/dnrp1-1_merged-DpnII-NlaIII-a_trimmed.sorted.bam_clean.wig",
"dnrp1_2": "../data/dnrp1_2/dnrp1-2_merged-DpnII-NlaIII-a_trimmed.sorted.bam_clean.wig",
"bem1-aid_a":"../data/bem1-aid_a/yWT03a_16_trimmed_out_restriction_sites_yWT03a_16_merged_cleaned_forward_reads_trimmed.sorted.bam_clean.wig",
"bem1-aid_b":"../data/bem1-aid_b/yWT0321_a_trimmed_out_restriction_sites_yWT0321_a_merged_cleaned_forward_reads_trimmed.sorted.bam_clean.wig",
"dbem1dbem3_a": "../data/dbem1dbem3_a/yTW001_4_merged_cleaned_forward_reads_trimmed.sorted.bam_clean.wig",
"dbem1dbem3_b": "../data/dbem1dbem3_b/yTW001_6_merged_cleaned_forward_reads_trimmed.sorted.bam_clean.wig",
"dbem3_a": "../data/dbem3_a/all_cleaned_fw_reads_trimmed.sorted.bam_clean.wig",
"dbem3_b": "../data/dbem3_b/yLIC137_8_merged_cleaned_forward_reads_trimmed.sorted.bam_clean.wig",}


## Import wig files 
wig_reads=defaultdict(dict)


for wigfile in wigfiles.keys():
    tmp=[]
    for x in wg.walk(open(wigfiles[wigfile])):
        tmp.append(x)
    wig_reads[wigfile]=tmp
    

## reading the wig files

data_wigfiles=defaultdict(dict)
for wigfile in wig_reads.keys():
    reads=[]
    pos=[]
    chrom=[]
    for i in wig_reads[wigfile]:
        reads.append(i[2])
        pos.append(i[1])
        chrom.append(i[0])
    data_wigfiles[wigfile]["Reads"]=reads
    data_wigfiles[wigfile]["Positions"]=pos
    data_wigfiles[wigfile]["Chromosomes"]=chrom



# +
data_wigfiles_pd=pd.DataFrame(data_wigfiles)

data_wigfiles_pd=data_wigfiles_pd.T
# -

data_wigfiles_pd

# +
## Compute the sum of the reads and length of Positions for each key

data_wigfiles_pd["Reads_sum"]=data_wigfiles_pd["Reads"].apply(lambda x: sum(x))
data_wigfiles_pd["Insertions_sum"]=data_wigfiles_pd["Positions"].apply(lambda x: len(x))
# -

data_wigfiles_pd

np.sum(data_wigfiles_pd["Reads"]["wt_a"]),np.sum(data_wigfiles_pd["Reads"]["wt_b"]),np.sum(data_wigfiles_pd["Reads"]["dnrp1_1"]),np.sum(data_wigfiles_pd["Reads"]["dnrp1_2"])

len(data_wigfiles_pd["Positions"]["wt_a"]),len(data_wigfiles_pd["Positions"]["wt_b"]),len(data_wigfiles_pd["Positions"]["dnrp1_1"]),len(data_wigfiles_pd["Positions"]["dnrp1_2"])

list_data_pd.loc["wt_a","Reads"].sum(),list_data_pd.loc["wt_b","Reads"].sum(),list_data_pd.loc["dnrp1_1","Reads"].sum(),list_data_pd.loc["dnrp1_2","Reads"].sum()

list_data_pd.loc["wt_a","Insertions"].sum(),list_data_pd.loc["wt_b","Insertions"].sum(),list_data_pd.loc["dnrp1_1","Insertions"].sum(),list_data_pd.loc["dnrp1_2","Insertions"].sum()

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
#     - subdivide the gene into known domains and compute the fitness of individual domains, and take the average as fitness . For genes with no domains then the fitness is from the whole gene. 
#     - assign a new corrected fitness to the every gene by taking the fitness of the domain with strongest effect with respect to the average , so the max|f_mean-f_domain_i|
#  
#  Then :
#  - Plot the distributions of every fitness calculation referenced to the median fitness value of the population.
#  - Plot the average fitness vs the corrected fitness 
#
#
#

list_data_pd.loc["wt_merged"].sort_values(by="Reads",ascending=False).head(10)

# ## Computing the number of reads per insertion along the gene length for all genes 

# +
# data=list_data_pd.loc["wt_merged"]
data=list_data_pd.loc["wt_a"]
## remove ade2 data from the dataset
# data.drop(labels=5805,axis=0,inplace=True)
# ## remove the ura3 data from the dataset
# data.drop(labels=1927,axis=0,inplace=True)


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
        



# +
# compute the reads per insertions for each gene along the gene length

gene_parts=np.linspace(0,1,11)
r=np.zeros(shape=(len(insertion_locations),len(gene_parts))) # reads_per_insertion_parts array , every gene in the rows 
insertions_array=np.zeros(shape=(len(insertion_locations),len(gene_parts))) # insertions_parts array , every gene in the rows
reads_array=np.zeros(shape=(len(insertion_locations),len(gene_parts))) # reads_parts array , every gene in the rows
for i in np.arange(0,len(insertion_locations)):
    if (insertion_locations[i])!=0:
        g=np.array(insertion_locations[i]) # insertion locations
        f=np.linspace(gene_coordinates[i][0],gene_coordinates[i][1],len(gene_parts)) # gene coordinates
        #f=np.array(gene_coordinates[i][0]+gene_coordinates[i][1]*gene_parts)
        binedges = g.searchsorted(f)

        rngs = [list(range(binedges[k], binedges[k+1])) for k in range(len(binedges)-1)]
        

        for k in np.arange(0,len(rngs)):
            readsperinsert=[]
            insertions_array[i,k]=len(rngs[k])
            for j in np.arange(0,len(rngs[k])):
                readsperinsert.append(reads_locations[i][rngs[k][j]])
            reads_array[i,k]=np.sum(readsperinsert) # number of reads binned in 10 equally sized bins
            
                
            if len(readsperinsert)>1:
                r[i,k]=np.sum(readsperinsert)/(len(readsperinsert)-1)#discarding the insertion with the highest read count
            else:
                r[i,k]=0

    



# +
data_new=data.copy()
data_new.index=data_new["Gene name"]
r_essentials=[]
for i in standard_essentials:
    if i in data_new.index:
        #ind=np.where(np.array(data.loc[:,"Gene name"])==i)[0]
        ind=data_new.index.get_loc(i)
        r_essentials.append(r[ind,:])
        # if ind.size>0:
        #     r_essentials.append(r[ind[0],:])

r_essentials=np.array(r_essentials)

# +
# replace nan values with zeros

r[np.isnan(r)]=0

## remove data from ADe2 and URA3 genes

# r_new=np.delete(r,5805,axis=0)
# r_new=np.delete(r_new,1927,axis=0)

# r[5805,:]=1
# r[1927,:]=1
# Take the median of  all values of r along the rows

r_sum=np.mean(r,axis=0)

# Plot r values along gene parts

ax=plt.figure(figsize=(8,5))

# for i in np.arange(1,len(insertion_locations)):
#     plt.plot(gene_parts,r[i,:],color="black",alpha=0.3,marker="o",markersize=1)

# plt.plot(gene_parts,r[0,:],color="black",alpha=0.3,marker="o",markersize=1,
#     label="Data per gene")
# plt.yscale("log")
plt.xlabel("Gene length",fontsize=20)
plt.ylabel("Reads per insertion per gene",fontsize=20)

plt.bar(gene_parts,r_sum,color="gray",alpha=0.7,width=gene_parts[-1]/len(gene_parts),
label="Mean of all genes")
plt.xticks(np.linspace(-0.05,0.95,11),
labels=["0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"],
fontsize=16);
plt.yticks(fontsize=20);

# plt.legend(fontsize=16)
# plt.grid(linewidth=0.5)
plt.tight_layout()

#plt.savefig("../figures/fig_reads_per_insertion_all_genes_along_gene_length.png",dpi=300,format="png")



# +
# replace nan values with zeros

r_essentials[np.isnan(r_essentials)]=0

# Sum all values of r along the rows

r_sum=np.mean(r_essentials,axis=0)

# Plot r values along gene parts

ax=plt.figure(figsize=(8,5))

# for i in np.arange(1,len(insertion_locations)):
#     plt.plot(gene_parts,r[i,:],color="black",alpha=0.3,marker="o",markersize=1)

# plt.plot(gene_parts,r[0,:],color="black",alpha=0.3,marker="o",markersize=1,
#     label="Data per gene")
# plt.yscale("log")
plt.xlabel("Gene length",fontsize=20)
plt.ylabel("Reads per insertion per gene",fontsize=20)

plt.bar(gene_parts,r_sum,color="pink",alpha=0.7,width=gene_parts[-1]/len(gene_parts),
label="Mean of all essential genes")
plt.xticks(np.linspace(-0.05,0.95,11),
labels=["0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"],
fontsize=16);
plt.yticks(fontsize=20);
# plt.legend(fontsize=16)
# plt.grid(linewidth=0.5)
plt.tight_layout()

#plt.savefig("../figures/fig_reads_per_insertion_essential_genes_along_gene_length.png",dpi=300,format="png")
## from the figure it follows that I should remove only the first 10% of the reads to compute the fitness of the whole gene

# -

r_sum

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
windows_size=20000
chrom_conversion={"chrI":"I","chrII":"II","chrIII":"III",
"chrIV":"IV","chrV":"V","chrVI":"VI","chrVII":"VII","chrVIII":"VIII",
"chrIX":"IX","chrX":"X","chrXI":"XI","chrXII":"XII","chrXIII":"XIII",
"chrXIV":"XIV","chrXV":"XV","chrXVI":"XVI"}

chromosome_density=[]

for i in chrom_conversion.keys():
    tmp=len(pos[chrom==i])/chromosome_size[chrom_conversion[i]]
    chromosome_density.append(tmp)

expected_flanking_enrichment=np.round(np.multiply(chromosome_density,windows_size),0).astype(int)# this means 5kb flanking regions each side

expected_flanking_enrichment

# +

data=list_data_pd.loc["wt_merged"]


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



# +
flanking_regions_data_pd=pd.DataFrame.from_dict(flanking_regions_data,orient="index")

discarded_genes2fitness=flanking_regions_data_pd[flanking_regions_data_pd["classification"]=="Not enough flanking regions"].index
# -

#

flanking_regions_data_pd[flanking_regions_data_pd["classification"]=="Not enough flanking regions"].sort_values(by="insertions_flanking")

flanking_regions_data_pd[flanking_regions_data_pd["classification"]=="OK"].sort_values(by="insertions_flanking",ascending=False)[10:20]

# +

from functions_satay_biases import transposon_bias2centromeres

fig, distance2cent_all = transposon_bias2centromeres(wigfile,centromeres_file,save=False)

# -

## Divide distance to chromosome centromeres by the length of each chromosome 
distance2cent_all_density={}
for i in distance2cent_all.keys():
    distance2cent_all_density[i]=distance2cent_all[i]*2/chromosome_size[chrom_conversion[i]]

# +
## Distance to centromere bias in every chromosome

distance2cent_all_chrom=(list(distance2cent_all.values()))
distance2cent_all_chrom_density=(list(distance2cent_all_density.values()))

figure=plt.subplots(4,4,figsize=(16,12))
plt.subplots_adjust(hspace=0.5,wspace=0.5)
plt.suptitle("Centromere bias in every chromosome",fontsize=20,y=0.95)
# using colormap plasma 
color = plt.get_cmap('jet') 

j=1
for i in np.arange(0,len(distance2cent_all_chrom)):
    plt.subplot(4,4,j)
    #plt.hist(distance2cent_all_chrom[i],bins=50,alpha=0.4,label=list(distance2cent_all.keys())[i],density=True,color=color(i/len(distance2cent_all_chrom)));
    #sns.histplot(distance2cent_all_chrom[i],bins=50,alpha=0.4,kde=True,label=list(distance2cent_all.keys())[i],color=color(i/len(distance2cent_all_chrom)));
    # sns.kdeplot(distance2cent_all_chrom[i],alpha=0.4,
    # label=list(distance2cent_all.keys())[i],color=color(i/len(distance2cent_all_chrom)));
    sns.kdeplot(distance2cent_all_chrom[i]/1000,alpha=0.4,
    label=list(distance2cent_all.keys())[i],color="black");

    plt.xlabel("Distance to centromere (kbp)")
    plt.ylabel("Transposons density")
    #plt.ylim(0,3000)
    plt.xlim(0,2500)
    plt.legend()

    j=j+1
  
# -

# ## fitness - malthusian model
#
#
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
#     - subdivide the gene into known domains and compute the fitness of individual domains, and take the average as fitness . For genes with no domains then the fitness is from the whole gene. 
#     - assign a new corrected fitness to the every gene by taking the fitness of the domain with strongest effect with respect to the average , so the max|f_mean-f_domain_i|
#  
#  Then :
#  - Plot the distributions of every fitness calculation referenced to the median fitness value of the population.
#  - Plot the average fitness vs the corrected fitness 

# +
import pickle
with open("../postprocessed-data/fitness_models_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

fitness_all_pd=pd.concat(b,axis=0,keys=keys)

standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str)

# -

from functions_interaction_computations import filter_fitness
data_fitness=filter_fitness(fitness_all_pd,backgrounds=keys,goi=["BEM1","BEM3","NRP1"],discard=["Not enough flanking regions"],set2zero=["Not enough reads",
    "Not enough insertions"],cols=["fitness_gene","fitness_domains_corrected"],essentiality=False) # drop genes that do not have enough reads or insertions

# +
## Plot the distribution of fitness values for the wild type for the corrected and non corrected fitness for all genes and essential genes


data_wt_non_corrected=data_fitness.loc["wt_merged","fitness_gene"]
data_wt_corrected=data_fitness.loc["wt_merged","fitness_domains_corrected"]

data_wt_non_corrected_essentials=data_wt_non_corrected[data_wt_non_corrected.index.isin(standard_essentials)]
data_wt_corrected_essentials=data_wt_corrected[data_wt_corrected.index.isin(standard_essentials)]

# +
# select the genes that are not essential

data_wt_non_corrected_non_essentials=data_wt_non_corrected[~data_wt_non_corrected.index.isin(standard_essentials)]
data_wt_corrected_non_essentials=data_wt_corrected[~data_wt_corrected.index.isin(standard_essentials)]

# +
plt.subplots(1,1,figsize=(8,5))


sns.kdeplot(data_wt_non_corrected_non_essentials,alpha=0.4,label="All genes",color="black",shade=True);
sns.kdeplot(data_wt_non_corrected_essentials,alpha=0.4,label="Essential genes",color="pink",shade=True);

plt.xlabel("Fitness non corrected by domains ")

plt.vlines(data_wt_non_corrected_non_essentials.mean(),0,2.5,linestyle="--",color="black",label="Mean fitness non corrected")
plt.vlines(data_wt_non_corrected_essentials.mean(),0,2.5,linestyle="--",color="pink",label="Mean fitness non corrected essentials")

#plt.savefig("../figures/fitness_distribution_non_corrected_wt_essentials.png",dpi=300)

# -

data_wt_non_corrected_essentials.mean(),data_wt_non_corrected_non_essentials.mean()

# +
plt.subplots(1,1,figsize=(8,5))


sns.kdeplot(data_wt_corrected_non_essentials,alpha=0.4,label="All genes",color="black",shade=True);
sns.kdeplot(data_wt_corrected_essentials,alpha=0.4,label="Essential genes",color="pink",shade=True);

plt.xlabel("Fitness corrected by domains ")

plt.vlines(data_wt_corrected_non_essentials.mean(),0,2.5,linestyle="--",color="black",label="Mean fitness corrected")
plt.vlines(data_wt_corrected_essentials.mean(),0,2.5,linestyle="--",color="pink",label="Mean fitness corrected essentials")

#plt.savefig("../figures/fitness_distribution_corrected_wt_essentials.png",dpi=300)
# -

data_wt_corrected_non_essentials.mean(),data_wt_corrected_essentials.mean()

essentials_pred=data_wt_non_corrected[data_wt_non_corrected<data_wt_non_corrected.mean()-data_wt_non_corrected.std()].index.tolist()

data_wt_corrected.max()/5

# +
## Construct a ROC curve for predicting essential genes based on the fitness .values()


from sklearn.metrics import roc_curve, auc

y=np.zeros(len(data_wt_non_corrected)) ## labels for ROC curve
y_pred=np.zeros(len(data_wt_non_corrected)) ## predicted values for ROC curve
#essentials_pred=data_wt_non_corrected[data_wt_non_corrected<data_wt_non_corrected.mean()-0.5*2*data_wt_non_corrected.std()].index.tolist()
essentials_pred=data_wt_non_corrected[data_wt_non_corrected<0.42].index.tolist()
y_corrected=np.zeros(len(data_wt_corrected)) ## labels for ROC curve
y_pred_corrected=np.zeros(len(data_wt_corrected)) ## predicted values for ROC curve
#essentials_pred_corrected=data_wt_corrected[data_wt_corrected<data_wt_corrected.mean()-0.5*2*data_wt_corrected.std()].index.tolist()
essentials_pred_corrected=data_wt_corrected[data_wt_corrected<0.42].index.tolist()

j=0
for i in data_wt_non_corrected.index:
    if i in standard_essentials.tolist():
        y[j]=1
    j+=1
#
j=0
for i in data_wt_non_corrected.index:
    if i in essentials_pred:
        y_pred[j]=1
    j+=1


j=0
for i in data_wt_corrected.index:
    if i in standard_essentials.tolist():
        y_corrected[j]=1
    j+=1
#
j=0
for i in data_wt_corrected.index:
    if i in essentials_pred_corrected:
        y_pred_corrected[j]=1
    j+=1
# -

data_wt_non_corrected.mean()-1*data_wt_non_corrected.std(),data_wt_corrected.mean()-1*data_wt_corrected.std()

len(y),len(y_corrected),len(standard_essentials)

# +
## to construc a confusiion matrix 
from sklearn.metrics import confusion_matrix

# confusion_matrix(y, y_pred)/len(y)*100


cm = confusion_matrix(y, y_pred)

ax= plt.subplot()
sns.heatmap(cm, annot=True, ax = ax,fmt="d",cmap="Blues"); #annot=True to annotate cells

ax.set_xlabel('Predicted labels');ax.set_ylabel('True labels');
ax.set_title('Confusion Matrix');
ax.xaxis.set_ticklabels(['Non-essential', 'Essential']); ax.yaxis.set_ticklabels(['Non-essential', 'Essential']);

#plt.savefig("../figures/confusion_matrix_non_corrected_fitness.png",dpi=300)


# +
## to construc a confusiion matrix 
from sklearn.metrics import confusion_matrix

# confusion_matrix(y, y_pred)/len(y)*100


cm = confusion_matrix(y_corrected, y_pred_corrected)

ax= plt.subplot()
sns.heatmap(cm, annot=True, ax = ax,fmt="d",cmap="Blues"); #annot=True to annotate cells

ax.set_xlabel('Predicted labels');ax.set_ylabel('True labels');
ax.set_title('Confusion Matrix');
ax.xaxis.set_ticklabels(['Non-essential', 'Essential']); ax.yaxis.set_ticklabels(['Non-essential', 'Essential']);

#plt.savefig("../figures/confusion_matrix_corrected_fitness.png",dpi=300)

# +
from sklearn import metrics

#fpr, tpr, thresholds = metrics.roc_curve(y, fitness2rocprob)
fpr, tpr, thresholds = metrics.roc_curve(y, y_pred)
area=metrics.auc(fpr,tpr)

figure,ax=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate",fontsize=16)
plt.ylabel("True Positive Rate",fontsize=16)
#plt.title("ROC of the  corrected fitness to predict essentiality",fontsize=16)

plt.plot(fpr, tpr ,label=f"AUC={area:.2f}",color="darkorange",lw=2)
ax.tick_params(axis="both",labelsize=16)
ax.legend(loc="lower right",fontsize=16)

#figure.savefig("../figures/fig_non_corrected_fitness_ROC_curve.png",dpi=300)

# +
from sklearn import metrics

#fpr, tpr, thresholds = metrics.roc_curve(y, fitness2rocprob)
fpr, tpr, thresholds = metrics.roc_curve(y_corrected, y_pred_corrected)
area=metrics.auc(fpr,tpr)

figure,ax=plt.subplots(nrows=1,ncols=1,figsize=(8,5))

plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate",fontsize=16)
plt.ylabel("True Positive Rate",fontsize=16)
#plt.title("ROC of the  corrected fitness to predict essentiality",fontsize=16)

plt.plot(fpr, tpr ,label=f"AUC={area:.2f}",color="darkorange",lw=2)
ax.tick_params(axis="both",labelsize=16)
ax.legend(loc="lower right",fontsize=16)

#figure.savefig("../figures/fig_corrected_fitness_ROC_curve.png",dpi=300)

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

keys= ['wt_merged','wt_a','wt_b','dnrp1_1','dnrp1_2']
data=list_data_pd.loc[keys] # Take only data from targeted genotypes


# +
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
figure.savefig("../figures/fig_differences_biological_replicates.png",dpi=300)



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
"dnrp1_2": "../data/dnrp1_2/dnrp1-2_merged-DpnII-NlaIII-a_trimmed.sorted.bam_clean.wig"}


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

plt.legend(fontsize=16)
plt.grid(linewidth=0.5)
plt.tight_layout()

plt.savefig("../figures/fig_reads_per_insertion_all_genes_along_gene_length.png",dpi=300,format="png")



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
plt.legend(fontsize=16)
plt.grid(linewidth=0.5)
plt.tight_layout()

plt.savefig("../figures/fig_reads_per_insertion_essential_genes_along_gene_length.png",dpi=300,format="png")
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

flanking_regions_data_pd

len(discarded_genes2fitness),discarded_genes2fitness

# +

np.where(discarded_genes2fitness=="OSW1")
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
  

# +
## import protein domains

## data generated in the notebook "analysis_proteins_domains.ipynb"
data_domains=pd.read_excel("../postprocessed-data/genomic-domains-wt.xlsx",index_col="Unnamed: 0")
#data_reads_domains=pd.read_excel("../postprocessed-data/reads-per-domain-all-backgrounds-new.xlsx",index_col="Unnamed: 0")
#data_insertions_domains=pd.read_excel("../postprocessed-data/insertions-per-domain-all-backgrounds-new.xlsx",index_col="Unnamed: 0")


## data from yeastmine
domains_names=pd.read_csv('../data/Domains_all_genes_protein_coordinates_yeastmine.tsv',sep="\t")
domains_names.index=domains_names["Gene Name"]


# +
protein_domains_data=defaultdict(dict)

for i in data.index:

    gene=data.loc[i,"Gene name"]
    
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

    if gene in discarded_genes2fitness:
        protein_domains_data[gene]["classification"]="Not enough flanking regions"
    else:
        protein_domains_data[gene]["classification"]="OK"
        


# -

protein_domains_data_pd=pd.DataFrame.from_dict(protein_domains_data,orient="index")
protein_domains_data_pd.head(3)

# +
# how many genes have domains?

protein_domains_data_pd[protein_domains_data_pd["domains coordinates"].isna()==False].shape[0]

# +
# compute the reads in each protein domain

domain_coordinates=protein_domains_data_pd["domains coordinates"]
domain_coordinates=domain_coordinates[domain_coordinates.isna()==False]

for i in domain_coordinates.index: # over all genes 
    i_index=protein_domains_data_pd.index.get_loc(i)
    b=data[data["Gene name"]==i]["Insertion locations"]
    b=np.array(from_excel_to_list(np.array(b)[0]))
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
                
    protein_domains_data_pd.loc[i,"reads_domain"]=totalreads
    protein_domains_data_pd.loc[i,"insertions_domain"]=totalinsertions

    


# -

protein_domains_data_pd

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

reads_domains=[]
for i in protein_domains_data_pd.index:
    tmp=protein_domains_data_pd.loc[i,"reads_domain"] # total reads of all domains in the gene 
    if type(tmp)==list: # the list indicates that that protein has annotated domains
        reads_domains.append(tmp[0])
    

# +
# r=np.delete(r,5805,axis=0) # deleting the ade2 gene
# r=np.delete(r,1927,axis=0) # deleting the ura3 gene
# r[5805,:]=np.zeros_like(r[5805,:])
# r[1927,:]=np.zeros_like(r[1927,:])
ref=np.log2(np.median(np.sum(r,axis=1))) # reference fitness, assumption: most genes are neutral in the wild type
# r are the reads per insertions of every gene divided in chunks . Here I am summing all the reads per gene and taking its median 

ref_domains=np.log2(np.median(reads_domains)) # reference fitness, assumption: most genes are neutral in the wild type
# -

ref,ref_domains,len(reads_domains),len(protein_domains_data_pd)

# +
fitness_models=defaultdict(dict)

data=list_data_pd.loc["wt_merged"]

for i in np.arange(0,len(data)):
    gene=data.loc[i,"Gene name"]
    if gene in discarded_genes2fitness:
        fitness_models[gene]["fitness_gene"]="Not enough flanking regions"
        fitness_models[gene]["fitness_gene_std"]="Not enough flanking regions"
        fitness_models[gene]["fitness_domains_vector"]="Not enough flanking regions"
        fitness_models[gene]["fitness_domains_average"]="Not enough flanking regions"
        fitness_models[gene]["fitness_domains_corrected"]="Not enough flanking regions"
        fitness_models[gene]["domains"]="Not enough flanking regions"
    
    else:
        if np.sum(r[i,1:9])!=0:
            fitness_models[gene]["fitness_gene"]=np.log2(np.sum(r[i,1:9]))/ref # getting the 80% central part of the reads per insertions
            fitness_models[gene]["fitness_gene_std"]=np.log2(np.std(r[i,1:9]))/ref
            if np.array(protein_domains_data_pd.loc[gene,"domains coordinates"]).size>1:
                nume=np.array(protein_domains_data_pd.loc[gene,"reads_domain"])
                deno=np.array(protein_domains_data_pd.loc[gene,"insertions_domain"])
                if np.sum(deno)==0 or np.sum(nume)==0:
                    deno=1
                    nume=1
                elif deno=='':
                    deno=1
                else:
                    tmp=np.log2(nume/deno)
                # replace nan values with zeros
                if type(tmp)!=np.float64:
                    tmp[np.isnan(tmp)] = 0
                else:
                    tmp=0
                fitness_models[gene]["fitness_domains_vector"]=tmp/ref_domains
            
                fitness_models[gene]["fitness_domains_average"]=np.mean(fitness_models[gene]["fitness_domains_vector"])
                fitness_models[gene]["fitness_domains_std"]=np.std(fitness_models[gene]["fitness_domains_vector"])
                fitness_models[gene]["domains"]="annotated"
                ## computing the corrected fitness as the fitness domain that has the maximum difference with the average fitnes
                tmp=np.absolute(fitness_models[gene]["fitness_domains_vector"]-fitness_models[gene]["fitness_gene"])
                index_max=np.where(tmp==np.max(tmp))[0]
                if index_max.size>1:
                    index_max=index_max[0]
                    fitness_models[gene]["fitness_domains_corrected"]=fitness_models[gene]["fitness_domains_vector"][index_max]
                else:
                    fitness_models[gene]["fitness_domains_corrected"]=fitness_models[gene]["fitness_domains_vector"][index_max]
            else:
                fitness_models[gene]["fitness_domains_vector"]=fitness_models[gene]["fitness_gene"]
                fitness_models[gene]["fitness_domains_std"]=fitness_models[gene]["fitness_gene_std"]
                fitness_models[gene]["fitness_domains_average"]=fitness_models[gene]["fitness_gene"]
                fitness_models[gene]["fitness_domains_corrected"]=fitness_models[gene]["fitness_gene"]
                fitness_models[gene]["domains"]="non annotated domains"
    

# +
fitness_models_pd=pd.DataFrame.from_dict(fitness_models,orient="index")
fitness_models_pd["fitness_domains_average"].replace(-np.inf,0,inplace=True)
fitness_models_pd["fitness_domains_corrected"].replace(-np.inf,0,inplace=True)
fitness_models_pd["fitness_domains_average"].replace(np.inf,1.5,inplace=True)
fitness_models_pd["fitness_domains_corrected"].replace(np.inf,1.5,inplace=True)
fitness_models_pd["fitness_gene"].replace(-np.inf,0,inplace=True)
fitness_models_pd["fitness_gene"].replace(np.inf,1.5,inplace=True)


# -

fitness_models_pd.loc["NRP1"],protein_domains_data_pd.loc["NRP1"]

fitness_annotated_domains=fitness_models_pd[fitness_models_pd["domains"]=="annotated"]
fitness_annotated_domains

fitness_annotated_domains.loc[:,"fitness_domains_corrected"][0][0]
f_corrected=[]
for i in fitness_annotated_domains.index:
    if type(fitness_annotated_domains.loc[i,"fitness_domains_corrected"])!=np.float64:
        f_corrected.append(fitness_annotated_domains.loc[i,"fitness_domains_corrected"][0])
    else:
        f_corrected.append(0)


# +
# plt.figure(figsize=(5,5))
plt.hist(fitness_annotated_domains.loc[:,"fitness_gene"],bins=100,label="whole gene");
plt.hist(fitness_annotated_domains.loc[:,"fitness_domains_average"],bins=100,label="average domains",alpha=0.6);
plt.hist(f_corrected,bins=100,label="fitness from domain with strongest effect",alpha=0.6);

plt.legend()

plt.xlabel("Fitness")

plt.ylabel("Number of genes")

plt.title("Fitness of genes with annotated domains")

#plt.savefig("../figures/fitness_annotated_domains_hist.png")
# -

plt.scatter(fitness_annotated_domains.loc[:,"fitness_gene"],f_corrected,alpha=0.3)
plt.xlabel("fitness whole gene")
plt.ylabel("fitness domain with strongest effect")
plt.xlim(0,1.5)
plt.ylim(0,1.5)
plt.plot([0,1.5],[0,1.5],color="black",marker="None",linestyle="--")
plt.savefig("../figures/fitness_whole_gene_vs_fitness_domain_with_strongest_effect.png",dpi=300)

# +
## Fitness of essentials genes

fitness_annotated_domains["Essential"]=np.zeros_like(fitness_annotated_domains.index)

for i in standard_essentials:
    if i in fitness_annotated_domains.index:
        fitness_annotated_domains.loc[i,"Essential"]=1
    else:
        continue

fitness_essential=fitness_annotated_domains[fitness_annotated_domains["Essential"]==1]

# -

fitness_essential

# +


f_corrected_essential=[]
for i in fitness_essential.index:
    if type(fitness_annotated_domains.loc[i,"fitness_domains_corrected"])!=np.float64:
        f_corrected_essential.append(fitness_annotated_domains.loc[i,"fitness_domains_corrected"][0])
    else:
        f_corrected_essential.append(fitness_annotated_domains.loc[i,"fitness_domains_corrected"])

plt.hist(f_corrected,bins=100,label="fitness from domain with strongest effect",alpha=0.6);
#plt.hist(fitness_essential.loc[:,"fitness_gene"],bins=100,label="whole gene essentials");
#plt.hist(fitness_essential.loc[:,"fitness_domains_average"],bins=100,label="average domains essentials",alpha=0.6);
plt.hist(f_corrected_essential,bins=100,label="fitness (essentials)from domain with strongest effect",alpha=0.6);

plt.legend()
# -

plt.scatter(fitness_essential.loc[:,"fitness_gene"],f_corrected_essential,alpha=0.3)
plt.xlabel("fitness whole gene")
plt.ylabel("fitness domain with strongest effect")
plt.xlim(0,1.5)
plt.ylim(0,1.5)
plt.plot([0,1.5],[0,1.5],color="black",marker="None",linestyle="--")


# +

plt.hist(fitness_annotated_domains.loc[:,"fitness_gene"],bins=100,label="whole gene",alpha=0.6);
plt.hist(fitness_essential.loc[:,"fitness_gene"],bins=100,label="whole gene essentials",alpha=0.6);

plt.legend()

# +
index_positive_domains=np.where(f_corrected_essential>fitness_essential.loc[:,"fitness_gene"])[0]

fitness_essential.iloc[index_positive_domains,:]


# +
plt.scatter(fitness_annotated_domains.loc[:,"fitness_gene"],fitness_annotated_domains.loc[:,"fitness_gene_std"],alpha=0.3)
plt.scatter(fitness_essential.loc[:,"fitness_domains_average"],fitness_essential.loc[:,"fitness_domains_std"],alpha=0.3)
plt.xlabel("fitness whole gene")
plt.ylabel("fitness standard deviation")

plt.plot([0,1.5],[0,1.5],color="black",marker="None",linestyle="--")
plt.xlim(0,1.5)
plt.ylim(0,1.5)

# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3.8.10 ('satay-dev')
#     language: python
#     name: python3
# ---

# +
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys
from collections import defaultdict
from ast import literal_eval

from from_excel_to_list import from_excel_to_list
from transposonmapper.statistics import volcano

# +


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


# +
## import essential genes used in transposonmapper

essentials_satay=pd.read_csv("../postprocessed-data/Cerevisiae_AllEssentialGenes_List.txt",header=0,sep="\t")

essentials_satay.columns=["gene name"]

# import conversion file from systematic names to standard names 
conversion=pd.read_csv("../postprocessed-data/from_systematic_genenames2standard_genenames.csv",
header=0,sep=",")

conversion.columns=["systematic name","standard name"]

# save the standard names of the essential genes in a systematic format
standard_essentials=[]
for names in essentials_satay.loc[:,"gene name"]:
    
    if names in conversion["systematic name"].values:
        standard_essentials.append(conversion.loc[conversion["systematic name"]==names]["standard name"].values[0])


# +
# import excel file with the normalized data

data_norm_pd=pd.read_excel("../postprocessed-data/data_norm_linear_transformation_per_background.xlsx",
engine='openpyxl',index_col="background")
data_norm_pd.drop(columns=["Unnamed: 0","Unnamed: 1"],inplace=True)

polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")
polarity_genes.fillna(0,inplace=True)

# +
# Importing fitness data from the intergenic model 

fitness_data=pd.read_excel("../postprocessed-data/fitness_coarse_grained_all_pd.xlsx",index_col="Unnamed: 0")

# +
satay_wt=fitness_data[fitness_data.loc[:,"background"]=="wt_merged"]
satay_wt.index=satay_wt.loc[:,"Gene name"]

#satay_wt_ho=satay_wt[satay_wt.loc[:,"Gene name"]=="HO"]
satay_wt_ho=satay_wt.loc["HO"]
satay_wt["fitness2HO"]=satay_wt["fitness"]/satay_wt_ho["fitness"]

# -

satay_nrp1=fitness_data[fitness_data.loc[:,"background"]=="dnrp1_merged"]
satay_nrp1.index=satay_nrp1.loc[:,"Gene name"]
satay_nrp1_ho=satay_nrp1.loc["HO"]#HO locus seems to not interact with Nrp1 according SGD so I could normalize the fitness values of the double mutants to this value as well 
## the normalization of fitnesses value sin the dnrp1 background should be to what is considered the "WT" which is the HO locus in the wt background
satay_nrp1["fitness2HO"]=satay_nrp1["fitness"]/satay_wt_ho["fitness"]
#satay_nrp1["fitness2wt"]=satay_nrp1["fitness"]/satay_wt["fitness"]# fitness of the double mutants in nrp1 background

fitness_wt_nrp1=satay_wt.loc["NRP1","fitness2HO"]

satay_nrp1.loc["POL31"],satay_wt.loc["POL31"]

# +
## Genetic interactions on nrp1 
## mutants such that the fitness of nrp1 double mutants-f(nrp1)f(mutant) is different than zero 
gi_score=defaultdict(dict)

for i in satay_wt.index:
    if satay_nrp1.loc[i,"fitness2HO"]!=np.inf or satay_nrp1.loc[i,"fitness2wt"]!=-np.inf:
        
        tmp=satay_nrp1.loc[i,"fitness2HO"]-fitness_wt_nrp1*satay_wt.loc[i,"fitness2HO"]
        gi_score[i]["score"]=tmp
    else:
        tmp=0 # No interaction by default because they cant be measured 
    
gi_score_pd=pd.DataFrame.from_dict(gi_score,orient="index")      

# +
## Removing inf from the score 

gi_score_pd[gi_score_pd.loc[:,"score"]==np.inf]=0
gi_score_pd[gi_score_pd.loc[:,"score"]==-np.inf]=0
gi_score_pd[gi_score_pd.loc[:,"score"]==np.nan]=0
#gi_score_pd.loc[gi_score_pd.loc[:,"score"]==np.inf,"score"]=0
# -

# Probable interval for GI
min_bound=float(np.min(gi_score_pd)/4)
max_bound=float(np.max(gi_score_pd)/2)
print("The minimum value of the interaction score is",float(np.min(gi_score_pd)),
"and the maximum is",float(np.max(gi_score_pd)))
print("A nice interval to search for nrp1 interactors are genes whose score is below",min_bound,
"and above",max_bound)


# +
## selecting same number of interactors 

N_gi=100
gi_score_pd_sorted_pos=gi_score_pd.sort_values(by="score",ascending=False)
gi_score_pd_sorted_pos.fillna(0,inplace=True)



gi_score_pd_sorted_neg=gi_score_pd.sort_values(by="score",ascending=True)
gi_score_pd_sorted_neg.fillna(0,inplace=True)

pos_interactors=gi_score_pd_sorted_pos[0:N_gi]
neg_interactors=gi_score_pd_sorted_neg[0:N_gi]



# -

gi_score_pd.loc["MEC1","score"]

# +
## selecting interactors based on their values with respect the interaction scores 

neg_interactors=gi_score_pd[gi_score_pd.loc[:,"score"]<min_bound] # presumably negative interactions
pos_interactors=gi_score_pd[gi_score_pd.loc[:,"score"]>max_bound] # presumably positive interactions

# -

gi_score_pd_sorted_pos[0:50]["score"].iloc[-1]

# +
## exploring the scores and interactors 

fig,ax=plt.subplots(1,1,figsize=(5,5))

N,bins,patches=ax.hist(gi_score_pd.loc[:,"score"],bins=50,color="red",alpha=0.4);

ax.set_xlabel("GI score")
ax.set_ylabel("Frequency")
ax.set_xlim(-1,1)
ax.set_title("Distribution of GI scores for Nrp1")

# ax.vlines(x=min_bound,ymin=0,ymax=1000,color="black",linestyle="--")
# ax.vlines(x=max_bound,ymin=0,ymax=1000,color="black",linestyle="--")
ax.text(x=-0.9,y=np.max(N)-30,s="Negative\ninteractors",rotation=0,fontsize=12)
ax.text(x=0.5,y=np.max(N)-30,s="Positive\ninteractors",rotation=0,fontsize=12)

ax.text(x=-0.9,y=np.max(N)-100,s=str(len(neg_interactors))+" genes",rotation=0,fontsize=12)
ax.text(x=0.5,y=np.max(N)-100,s=str(len(pos_interactors))+" genes",rotation=0,fontsize=12)

## if the interactors are selected based on their scores
# ni_patches=len(bins)-np.where(bins[bins>min_bound])[0]-1 # patches in the neutral part
# pi_patches=len(bins)-np.where(bins[bins>max_bound])[0]-1 # patches in the positive part

## if the interactors are sleected based on their number 
ni_patches=len(bins)-np.where(bins[bins>gi_score_pd_sorted_neg[0:N_gi]["score"].iloc[-1]])[0]-1 # patches in the neutral part
pi_patches=len(bins)-np.where(bins[bins>gi_score_pd_sorted_pos[0:N_gi]["score"].iloc[-1]])[0]-1 # patches in the positive part

for i in ni_patches:
    patches[i-1].set_facecolor('gray')
for i in pi_patches:
    patches[i-1].set_facecolor('green')

genes=["BEM1","BEM3","MEC1"]
i=0
j=0
colors=["black","black","black","purple","black"]
for gene in genes:
    ax.vlines(x=gi_score_pd.loc[gene,"score"],ymin=0,ymax=np.max(N)-200,color=colors[j],linestyle="--",alpha=0.5)
    ax.text(x=gi_score_pd.loc[gene,"score"],y=100+i,s=gene,rotation=0,fontsize=12,color=colors[j])
    i=i+50
    j=j+1


plt.tight_layout()

fig.savefig("../figures/fig_distribution_of_GI_scores_for_Nrp1_genes_annotated_same_number_interactors.png",dpi=300)

# +
## Plot the gene enrichment pie plot of the negative and positive interactors 

neg_interactors_name=neg_interactors.index
pos_interactors_name=pos_interactors.index

## Investigate how them correlate with our current knowledge about it 



# +
import gseapy as gp
from gseapy import barplot, dotplot
type_gi="Neg_GI_from_fitness_same_number_interactors"
goi=neg_interactors_name.tolist()
yeast = gp.get_library_name(organism='Yeast')

sets=[yeast[2],yeast[5],yeast[8],yeast[11],yeast[16] ] 
#['GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
#'Gene_Interaction_Hubs_BioGRID_2018','Pfam_Domains_2019']
# #%% enrichment 

for i in np.arange(0,len(sets)): 

  enr=gp.enrichr(gene_list=goi,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                
                  outdir='../postprocessed-data/enrich-analysis_dnrp1/'+type_gi+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
      
# to save your figure, make sure that ``ofname`` is not None
  ax = dotplot(enr.res2d, title=sets[i],cmap='viridis_r', size=20, figsize=(3,5),ofname=sets[i]+type_gi)
  
    
# -

sets

# +
# simple plotting function
from gseapy import barplot, dotplot
for i in np.arange(0,len(sets)):
    
# to save your figure, make sure that ``ofname`` is not None
    ax = dotplot(enr[i].res2d, title=sets[i],cmap='viridis_r', size=20, figsize=(3,5),ofname=sets[i])
# -

enr[3].res2d.head(2)

keys= ['wt_merged','dnrp1_merged','dbem3_merged','dbem1dbem3_a']

# ## Volcano plots

# +
path_a = r"../data/"
filelist_a = ["wt_a/WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_pergene.txt",
"wt_b/WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_pergene.txt"]
path_b = r"../data/"
filelist_b = ["dnrp1_a/dnrp1-1_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt",
"dnrp1_b/dnrp1-2_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt"]


variable = 'read_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.001 #set threshold above which p-values are regarded significant
normalize=True

trackgene_list = ['nrp1','bem3','bem1','bem2',"mec1"] # ["cdc42"]


figure_title = "$\Delta$nrp1 vs WT "

fc_interval=[1.5,-1.5]
pv_values=[3,3]

volcano_df_nrp1_wt = volcano(path_a=path_a, filelist_a=filelist_a,
            path_b=path_b, filelist_b=filelist_b, 
            fold_change_interval=fc_interval,p_value_interval=pv_values,
            variable=variable,
            significance_threshold=significance_threshold,
            normalize=normalize,
            trackgene_list=trackgene_list,
            figure_title=figure_title,savefigure=True)

# +
volcano_df_nrp1_wt.sort_values(by=['fold_change'], inplace=True)

dnrp1_genes_positive_enriched=volcano_df_nrp1_wt[volcano_df_nrp1_wt["significance"]==True][0:50]

dnrp1_genes_negative_enriched=volcano_df_nrp1_wt[volcano_df_nrp1_wt["significance"]==True][-50:]

volcano_df_nrp1_wt[volcano_df_nrp1_wt.loc[:,"gene_names"]=="CLA4"]
#dnrp1_genes_positive_enriched[dnrp1_genes_positive_enriched.loc[:,"gene_names"]=="whi3"]

# +
## Plot number of normalized insertions from WT/ total number of insertions of the library 
# vs the same for dnrp1

keys= ['wt_merged','dnrp1_merged']
data=data_norm_pd.loc[keys]

variable=["tr_normalized_windows" ,"reads_normalized_windows"]
var_norm=["Insertions" ,"Reads"]
x_1=data.loc["wt_merged"][variable[0]]/data.loc["wt_merged"][var_norm[0]].sum()
y_1=data.loc["dnrp1_merged"][variable[0]]/data.loc["dnrp1_merged"][var_norm[0]].sum()

x_0=data.loc["wt_merged"][variable[1]]/data.loc["wt_merged"][var_norm[1]].sum()
y_0=data.loc["dnrp1_merged"][variable[1]]/data.loc["dnrp1_merged"][var_norm[1]].sum()


fig , ax = plt.subplots(ncols=2,figsize=(8,3))

plt.subplots_adjust(wspace=0.4)

# sns.regplot(x_0,y_0,fit_reg=True,color="black",marker="o",ax=ax[0],label="Normalized \n Insertions",
#             scatter_kws={"s":30,"alpha":0.2})
# sns.regplot(x_1,y_1,fit_reg=True,color="black",marker="o",ax=ax[1],label="Normalized \n Reads",
#             scatter_kws={"s":30,"alpha":0.2})

ax[0].scatter(x_0,y_0,color="black",marker="o",s=30,alpha=0.2,label="Normalized \n Insertions")

ax[1].scatter(x_1,y_1,color="black",marker="o",s=30,alpha=0.2,label="Normalized \n Reads")

for axes in ax:
    axes.set_xscale("log")
    axes.set_yscale("log")    
    xmin,xmax = axes.get_xlim()
    axes.set_ylim(xmin,xmax)
    axes.set_xlim(xmin,xmax)
    axes.set_xlabel("WT-log",fontsize=14)

    axes.set_ylabel("$\Delta$nrp1-log",fontsize=14)
    axes.tick_params(axis='both', which='major', labelsize=16)
    axes.legend()

# -

test=["wt_merged","bem1-aid_b"]
data_test=data_norm_pd.loc[test]

# +
variable=["tr_normalized_windows" ,"reads_normalized_windows"]
var_norm=["Insertions" ,"Reads"]
data=data_test

x_1=data.loc["wt_merged"][variable[0]]/data.loc["wt_merged"][var_norm[0]].sum()
y_1=data.loc["bem1-aid_b"][variable[0]]/data.loc["bem1-aid_b"][var_norm[0]].sum()

x_0=data.loc["wt_merged"][variable[1]]/data.loc["wt_merged"][var_norm[1]].sum()
y_0=data.loc["bem1-aid_b"][variable[1]]/data.loc["bem1-aid_b"][var_norm[1]].sum()
# x_1=data_test.loc["wt_merged"]["reads_over_windows"]
# y_1=data_test.loc[test[1]]["reads_over_windows"]

# x_0=data_test.loc["wt_merged"]["insertions_over_windows"]
# y_0=data_test.loc[[test[1]]]["insertions_over_windows"]
# x=data.loc["wt_merged"]["Reads"]
# y=data.loc["dnrp1_merged"]["Reads"]

fig , ax = plt.subplots(ncols=2,figsize=(6,2))

plt.subplots_adjust(wspace=0.5)

# sns.regplot(x_0,y_0,fit_reg=True,color="black",marker="o",ax=ax[0])
# sns.regplot(x_1,y_1,fit_reg=True,color="black",marker="o",ax=ax[1])

ax[0].scatter(x_0,y_0,color="black",marker="o",alpha=0.5,label="Insertions")

ax[1].scatter(x_1,y_1,color="black",marker="o",alpha=0.5,label="Reads")


for axes in ax:
    # axes.set_xscale("log")
    # axes.set_yscale("log")    
    # axes.set_xlim(0,0.01)
    # axes.set_ylim(0,0.01)
    axes.set_xlabel("WT")
    axes.set_ylabel("$\Delta$bem1")
    axes.legend()


# +

def pie_chart_enrichment(data,process,type,savefig=False):

    
    

    terms=[]
    for i in data.loc[:,"Term"].tolist():
        terms.append(i.split(" (")[0])

    ## Data to plot
    
    data2plot = data.loc[:,"Combined Score"][0:10]
    labels = terms[0:10]

    #define Seaborn color palette to use
    colors = sns.color_palette('pastel')[0:10]


    #create pie chart
    fig,axes=plt.subplots(1,1,figsize=(10,10))
    #patches,texts=plt.pie(data2plot, labels = labels, colors = colors, autopct='%.0f%%',textprops={'fontsize': 25});
    
    explode = list()
    for k in labels:
        explode.append(0.1)
        
    pie = plt.pie(data2plot, explode=explode, shadow=True, autopct='%1.1f%%', colors=colors)
    plt.legend(pie[0], labels, loc="best",fontsize=12)
    
    plt.tight_layout()
    if savefig==True:
        plt.savefig("../postprocessed-data/enrich-analysis_dnrp1/pie_chart_enrichment_"+process+"_"+type+".png")

    return 

# +
## Enrichment analysis with gseapy library of genes that have different enrichment score in the volcano plot:
import gseapy as gp 
type_gi="Negative"
goi=dnrp1_genes_negative_enriched.loc[:,"gene_names"].tolist()
yeast = gp.get_library_name(organism='Yeast')
yeast
sets=[yeast[2],yeast[5],yeast[8],yeast[11],yeast[16] ] #['GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
#'Gene_Interaction_Hubs_BioGRID_2018','Pfam_Domains_2019']
# #%% enrichment 
for i in np.arange(0,len(sets)): 

  enr = gp.enrichr(gene_list=goi,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                  description=type_gi +'_enrichment_dnrp1',
                  outdir='../postprocessed-data/enrich-analysis_dnrp1/'+type_gi+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                )



# +
## Import analyzed data and plot it :
type_gi="Negative"
out_dir='../postprocessed-data/enrich-analysis_dnrp1/'+type_gi+'/'
filename_bio="GO_Biological_Process_2018.Yeast.enrichr.reports.txt"
filename_cell="GO_Cellular_Component_2018.Yeast.enrichr.reports.txt"
filename_mol="GO_Molecular_Function_2018.Yeast.enrichr.reports.txt"

path=out_dir+filename_bio
biological_process=pd.read_csv(path,sep="\t")

path=out_dir+filename_cell
cellular_component=pd.read_csv(path,sep="\t")

path=out_dir+filename_mol
molecular_function=pd.read_csv(path,sep="\t")

pie_chart_enrichment(biological_process,"biological_process",type=type_gi,savefig=True)

# +
goi_systematic=[]

goi=dnrp1_genes_negative_enriched.loc[:,"gene_names"].tolist()

for i in goi:
    tmp=conversion[conversion.loc[:,"standard name"]==i]["systematic name"].values
    if len(tmp)>0:
        goi_systematic.append(tmp[0])
    else:
        goi_systematic.append(i)

# +
type_int="NRP1 GI NG" # NRP1 GI NG
cell_map_goi=pd.read_excel("../postprocessed-data/cell-map-data-on-NRP1.xlsx",header=0,
sheet_name=type_int)

cell_map_goi=cell_map_goi.loc[:,"ORF"].tolist()

goi_standard=[]

for i in cell_map_goi:
    tmp=conversion[conversion.loc[:,"systematic name"]==i]["standard name"].values
    if len(tmp)>0:

        if type(tmp[0])!=float:
        
            goi_standard.append(tmp[0])
            
        else:
            goi_standard.append(i)
       
    
    else:
        goi_standard.append(i)
# -

yeast = gp.get_library_name(organism='Yeast')
yeast
sets=[yeast[2],yeast[5],yeast[8],yeast[11],yeast[16] ] #['GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
#'Gene_Interaction_Hubs_BioGRID_2018','Pfam_Domains_2019']
# #%% enrichment
enr=[] 
for i in np.arange(0,len(sets)): 

  enr.append( gp.enrichr(gene_list=goi_standard,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                  outdir='../postprocessed-data/enrich-analysis_dnrp1/cellmap/'+type_int+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                ))

# +
# simple plotting function
from gseapy import barplot, dotplot
for i in np.arange(0,len(sets)):
    
# to save your figure, make sure that ``ofname`` is not None
    ax = dotplot(enr[i].res2d, title=sets[i],cmap='viridis_r', size=20, figsize=(3,5),ofname=sets[i])
# -

biological_process=pd.read_csv("../postprocessed-data/enrich-analysis_dnrp1/cellmap/"+type_int+"/"+ "GO_Biological_Process_2018.Yeast.enrichr.reports.txt",sep="\t")
pie_chart_enrichment(biological_process,"biological_process_cell_map",type=type_int,savefig=True)

# +
## Analysis of essential genes in WT and dnrp1
# - import the genes that are duplicated 
# - import the genes that have low insertions in the flanking regions

# To extract potential essential genes take :

# - genes with low normalized transposons counts (given by a 10kb window) NI
# - genes with the longest region void  of transposons (FI) 
# - Highest domain likelihood scores (Benoit score) 

# Compare low NI with known essential genes.
# -

data_wt=data_norm_pd.loc["wt_merged"]
data_wt.reset_index(inplace=True)
data_wt[data_wt.loc[:,"Gene name"]==standard_essentials[0]].index[0]

# +
data_wt=data_norm_pd.loc["wt_merged"].copy()
data_wt.reset_index(inplace=True)
data_wt["true essential"]=np.zeros(len(data_wt))

for true_genes in standard_essentials:
    
    index_essential=data_wt[data_wt.loc[:,"Gene name"]==true_genes].index
    if len(index_essential)>0:
        if type(index_essential[0])==np.int64:
            #scores_wt.loc[index_essential[0],"true essential"]=1
            tmp=data_wt.index[index_essential]

            data_wt.loc[tmp,"true essential"]=1
    
    
data_wt.fillna(0,inplace=True)
# -

data_wt.groupby(["true essential"]).count()

# +


g=sns.pairplot(data_wt,hue="true essential",vars=["Reads","tr_normalized_windows",
"reads_normalized_windows","tr-density","tr_normalized_windows"],
diag_kind="kde",diag_kws={"shade":True},plot_kws = {'alpha': 0.6, 's': 80, 'edgecolor': 'k'},
size = 4)

g.axes[0,0].set_xlim(0,10000)
# g.axes[0,1].set_xlim(0,100)
g.axes[0,1].set_xlim((0,0.05))
g.axes[0,2].set_xlim((0,0.05))
g.axes[0,3].set_xlim((0,0.2))

plt.tight_layout()
# -



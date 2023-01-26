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
#     display_name: Python 3.9.7 ('transposonmapper')
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

from scipy import stats

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

# -

keys,len(keys)

# +
import pickle
with open("../postprocessed-data/fitness_models_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

fitness_all_pd=pd.concat(b,axis=0,keys=keys)

# -

standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str)

# +

polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")
polarity_genes.fillna(0,inplace=True)

# +
data=[]
for backg in keys:
    f=fitness_all_pd.loc[backg]
    f=f[f.loc[:,"fitness_gene"]!="Not enough flanking regions"]
    data.append(f)

data_fitness=pd.concat(data,axis=0,keys=keys)

# -

data_fitness.loc["wt_merged","NRP1"]

# +
## Computing the interactors of NRP1 using the whole fitness of the gene 
gene="NRP1"
nrp1_f_a=data_fitness.loc["wt_merged",gene]["fitness_gene"] # in wt_A , nrp1 is not found
nrp1_f_b=data_fitness.loc["wt_b",gene]["fitness_gene"]
data_b=data_fitness.loc["wt_b"]
data_a=data_fitness.loc["wt_a"]

intersection_genes=list((set(data_b.index)&set(data_a.index)))

significance_threshold = 0.05 #set significance threshold
gi=defaultdict(dict)
ttest_tval_list = [np.nan]*2 #initialize list for storing t statistics
ttest_pval_list = [np.nan]*2 #initialize list for storing p-values
signif_thres_list = False #initialize boolean list for indicating datapoints with p-value above threshold
fc_list = [np.nan]*2
for gene in intersection_genes :
    geneX=gene
    
    geneX_f_a=data_fitness.loc["wt_a",geneX]["fitness_gene"]
    if geneX in data_fitness.loc["wt_b"].index:
        geneX_f_b=data_fitness.loc["wt_b",geneX]["fitness_gene"]
        if geneX in data_fitness.loc["dnrp1_1"].index and geneX in data_fitness.loc["dnrp1_2"].index:
            geneXnrp1_f_a=data_fitness.loc["dnrp1_1",geneX]["fitness_gene"]
            geneXnrp1_f_b=data_fitness.loc["dnrp1_2",geneX]["fitness_gene"]
            
            variable_a_array=[geneXnrp1_f_a,geneXnrp1_f_b]
            variable_b_array=[geneX_f_a*nrp1_f_a,geneX_f_b*nrp1_f_b]
            
            ttest_val = stats.ttest_ind(variable_a_array, variable_b_array) #T-test
            gi[gene]["p_statistic"]=ttest_val[1]
            ttest_tval_list = ttest_val[0]
            gi[gene]["t_statistic"]=ttest_tval_list
            if not ttest_val[1] == 0: #prevent p=0 to be inputted in log
                ttest_pval_list = -1*np.log10(ttest_val[1])
                
            else:
                ttest_pval_list = 0
            gi[gene]["p_value"]=ttest_pval_list
            if ttest_pval_list >= -1*np.log10(significance_threshold):
                gi[gene]["significance"]=True
            else:
                gi[gene]["significance"]=False
            #DETERMINE FOLD CHANGE PER GENE
            if np.mean(variable_b_array) == 0 or np.mean(variable_a_array) == 0:
                fc_list = 0
            else:
                #fc_list = np.log2(np.mean(variable_a_array) / np.mean(variable_b_array))
                fc_list=np.mean(variable_a_array)-np.mean(variable_b_array)
            gi[gene]["fold_change"]=fc_list

        


# +
gi_pd_fitness_gene=pd.DataFrame.from_dict(gi,orient="index")

gi_pd_fitness_gene["gene_names"]=gi_pd_fitness_gene.index

gi_pd_fitness_gene
# -

gi_pd_fitness_gene[gi_pd_fitness_gene.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)

# +
## Now use the fitness correction by domains to find GI 
gene="NRP1"
nrp1_f_a=data_fitness.loc["wt_merged",gene]["fitness_domains_average"] # in wt_A , nrp1 is not found
nrp1_f_b=data_fitness.loc["wt_b",gene]["fitness_domains_average"]
data_b=data_fitness.loc["wt_b"]
data_a=data_fitness.loc["wt_a"]

intersection_genes=list((set(data_b.index)&set(data_a.index)))



significance_threshold = 0.05 #set significance threshold
gi=defaultdict(dict)
ttest_tval_list = [np.nan]*2 #initialize list for storing t statistics
ttest_pval_list = [np.nan]*2 #initialize list for storing p-values
signif_thres_list = False #initialize boolean list for indicating datapoints with p-value above threshold
fc_list = [np.nan]*2
for gene in intersection_genes :
    geneX=gene
    
    geneX_f_a=data_fitness.loc["wt_a",geneX]["fitness_domains_corrected"]
    if type(geneX_f_a)!=np.float64:
        geneX_f_a=geneX_f_a[0]
        if geneX_f_a==0:
            geneX_f_a=data_fitness.loc["wt_a",geneX]["fitness_domains_average"]
            if geneX_f_a==0:
                geneX_f_a=data_fitness.loc["wt_a",geneX]["fitness_gene"]
    else:
        geneX_f_a=geneX_f_a
        if geneX_f_a==0:
            geneX_f_a=data_fitness.loc["wt_a",geneX]["fitness_domains_average"]
            if geneX_f_a==0:
                geneX_f_a=data_fitness.loc["wt_a",geneX]["fitness_gene"]
    
    if geneX in data_fitness.loc["wt_b"].index:
        geneX_f_b=data_fitness.loc["wt_b",geneX]["fitness_domains_corrected"]
        if type(geneX_f_b)!=np.float64:
            geneX_f_b=geneX_f_b[0]
            if geneX_f_b==0:
                geneX_f_b=data_fitness.loc["wt_b",geneX]["fitness_domains_average"]
                if geneX_f_b==0:
                    geneX_f_b=data_fitness.loc["wt_b",geneX]["fitness_gene"]
        else:
            geneX_f_b=geneX_f_b
            if geneX_f_b==0:
                geneX_f_b=data_fitness.loc["wt_b",geneX]["fitness_domains_average"]
                if geneX_f_b==0:
                    geneX_f_b=data_fitness.loc["wt_b",geneX]["fitness_gene"]

        if geneX in data_fitness.loc["dnrp1_1"].index and geneX in data_fitness.loc["dnrp1_2"].index:
            geneXnrp1_f_a=data_fitness.loc["dnrp1_1",geneX]["fitness_domains_corrected"]
            geneXnrp1_f_b=data_fitness.loc["dnrp1_2",geneX]["fitness_domains_corrected"]
            if type(geneXnrp1_f_a)!=np.float64:
                geneXnrp1_f_a=geneXnrp1_f_a[0]
            if type(geneXnrp1_f_b)!=np.float64:
                geneXnrp1_f_b=geneXnrp1_f_b[0]
            
            
            variable_a_array=[geneXnrp1_f_a,geneXnrp1_f_b]
            variable_b_array=[geneX_f_a*nrp1_f_a,geneX_f_b*nrp1_f_b]
            
            ttest_val = stats.ttest_ind(variable_a_array, variable_b_array) #T-test
            gi[gene]["p_statistic"]=ttest_val[1]
            ttest_tval_list = ttest_val[0]
            gi[gene]["t_statistic"]=ttest_tval_list
            if not ttest_val[1] == 0: #prevent p=0 to be inputted in log
                ttest_pval_list = -1*np.log10(ttest_val[1])
                
            else:
                ttest_pval_list = 0
            gi[gene]["p_value"]=ttest_pval_list
            if ttest_pval_list >= -1*np.log10(significance_threshold):
                gi[gene]["significance"]=True
            else:
                gi[gene]["significance"]=False
            #DETERMINE FOLD CHANGE PER GENE
            if np.mean(variable_b_array) == 0 and np.mean(variable_a_array) == 0:
                fc_list = 0
            else:
                # fc_list = np.log2(np.mean(variable_a_array) / np.mean(variable_b_array))
                fc_list=np.mean(variable_a_array)-np.mean(variable_b_array)
            gi[gene]["fold_change"]=fc_list

        



# +
gi_pd=pd.DataFrame.from_dict(gi,orient="index")

gi_pd["gene_names"]=gi_pd.index

gi_pd[gi_pd.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)
# -

gi_pd.loc["BEM1"]

# +
## Distribution of fitness effects

data_fitness_wt=data_fitness.loc[["dnrp1_1","dnrp1_2"]]
data_fitness_wt.loc[:,"fitness_domains_corrected"].astype(float).hist(bins=100)
plt.title("Distribution of fitness effects in $\Delta$nrp1")
plt.xlabel("Fitness effect")
plt.ylabel("Number of genes")
plt.tight_layout()
#plt.savefig("../figures/fig_fitness_effects_dnrp1.png",dpi=300,transparent=True)

# +
from annotate_volcano import annotate_volcano   #import annotate_volcano function
volcano_df=gi_pd
fig=annotate_volcano(volcano_df,[0.7,-0.25],[1.5,0.5],figure_title="Interactors of nrp1 in WT",variable="fitness corrected")

#plt.savefig("../figures/fig_volcano_interactors_nrp1.png",dpi=300,transparent=True)


# +
gi_pd_significant=gi_pd[gi_pd.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)

neg_satay_signif=gi_pd_significant[gi_pd_significant.loc[:,"fold_change"]<0].loc[:,"gene_names"].index
neg_satay=gi_pd[gi_pd.loc[:,"fold_change"]<0].loc[:,"gene_names"].index

pos_satay_signif=gi_pd_significant[gi_pd_significant.loc[:,"fold_change"]>0].loc[:,"gene_names"].index
pos_satay=gi_pd[gi_pd.loc[:,"fold_change"]>0].loc[:,"gene_names"].index

# +
import gseapy as gp
from gseapy import barplot, dotplot
type_gi="Pos_GI_satay"
goi=pos_satay_signif.tolist()
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

neg_nrp1=pd.read_excel("../postprocessed-data/NRP1_genetic_interactions_filtered_by_negat.xlsx",skiprows=8)
pos_nrp1=pd.read_excel("../postprocessed-data/NRP1_genetic_interactions_filtered_by_pos.xlsx",skiprows=8)

neg_nrp1_genes=neg_nrp1.loc[:,"Interactor.1"]
pos_nrp1_genes=pos_nrp1.loc[:,"Interactor.1"]

# +
common_neg_signif=[]
for gene in neg_satay_signif:
    if gene in neg_nrp1_genes.unique():
        common_neg_signif.append(gene)

common_neg=[]
for gene in neg_satay:
    if gene in neg_nrp1_genes.unique():
        common_neg.append(gene)

common_pos_signif=[]
for gene in pos_satay_signif:
    if gene in pos_nrp1_genes.unique() :
        common_pos_signif.append(gene)

common_pos=[]
for gene in pos_satay:
    if gene in pos_nrp1_genes.unique() :
        common_pos.append(gene)

# +
consistent_neg_interactors_nrp1=list(set(neg_satay_signif) & set(neg_nrp1_genes) )

consistent_pos_interactors_nrp1=list(set(pos_satay_signif) & set(pos_nrp1_genes))

consistent_pos_interactors_nrp1,consistent_neg_interactors_nrp1

# +
from matplotlib_venn import venn2, venn2_circles

plt.figure(figsize=(10,5))
venn2(subsets = (len(pos_satay_signif), len(pos_nrp1_genes), len(consistent_pos_interactors_nrp1)), 
set_labels = ('SATAy', 'existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(pos_satay_signif), len(pos_nrp1_genes), len(consistent_pos_interactors_nrp1)));
plt.title("Consistency in positive interactions")
#
plt.savefig("../figures/venn_pos_nrp1.png",dpi=300)
# -

plt.figure(figsize=(10,5))
venn2(subsets = (len(neg_satay_signif), len(neg_nrp1_genes), len(consistent_neg_interactors_nrp1)), 
set_labels = ('SATAy', 'existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(neg_satay_signif), len(neg_nrp1_genes), len(consistent_neg_interactors_nrp1)));
plt.title("Consistency in negative interactions")
plt.savefig("../figures/venn_neg_nrp1.png",dpi=300)

# +
pos_satay_neg_costanzo=list(set(pos_satay_signif) & set(neg_nrp1_genes) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(pos_satay_signif), len(neg_nrp1_genes), len(pos_satay_neg_costanzo)), 
set_labels = ('PI SATAy', 'NI existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(pos_satay_signif), len(neg_nrp1_genes), len(pos_satay_neg_costanzo)));
plt.title("Misclasification I")

#plt.savefig("../figures/venn_misclass1.png",dpi=300)

# +
neg_satay_pos_costanzo=list(set(neg_satay_signif) & set(pos_nrp1_genes) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(neg_satay_signif), len(pos_nrp1_genes), len(neg_satay_pos_costanzo)), 
set_labels = ('NI SATAy', 'PI existing'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(neg_satay_signif), len(pos_nrp1_genes), len(neg_satay_pos_costanzo)));
plt.title("Misclasification II")

#plt.savefig("../figures/venn_misclass2.png",dpi=300)

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



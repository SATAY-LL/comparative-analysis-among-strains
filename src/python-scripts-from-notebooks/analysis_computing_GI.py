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

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys
from collections import defaultdict
from ast import literal_eval
from scipy.stats import norm
from scipy import stats

keys=['dbem3_b',
 'WT_1-Benoit',
 'dnrp1_2',
 'bem1-aid_a',
 'dbem1dbem3_b',
 'wt_merged',
 'dbem3_a-trimmed',
 'dbem1dbem3_a',
 'bem1-aid-dbem3_a',
 'bem1-aid_merged',
 'bem1-aid-dbem3_b',
 'dnrp1_1',
 'wt_b',
 'wt_a',
 'dnrp1_merged',
 'WT_2-Benoit',
 'bem1-aid_b',
 'dbem3_merged',
 'dbem3_a']

# +
import pickle
with open("../postprocessed-data/fitness_models_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

fitness_all_pd=pd.concat(b,axis=0,keys=keys)


with open("../postprocessed-data/discarded_genes_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

discarded_all=b

# +
discarded_genes_bem1=np.loadtxt("../postprocessed-data/discarded_genes_dbem1.txt",dtype=str)
discarded_genes_wt=np.loadtxt("../postprocessed-data/discarded_genes_wt.txt",dtype=str)
discarded_genes_bem13=np.loadtxt("../postprocessed-data/discarded_genes_dbem13.txt",dtype=str)


standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str)

polarity_genes=pd.read_csv("../postprocessed-data/polarity_genes_venn_Werner.txt",index_col="Gene")
polarity_genes.fillna(0,inplace=True)
# -

backgrounds=["wt_a","wt_b","bem1-aid_a","bem1-aid_b","dbem1dbem3_a","dbem1dbem3_b","dbem3_a","dbem3_b","dnrp1_1","dnrp1_2"]

fitness_backg=fitness_all_pd.loc[backgrounds]

# +
data=[]
for keys in backgrounds:
    f=fitness_backg.loc[keys]
    f=f[f.loc[:,"fitness_gene"]!="Not enough flanking regions"]
    data.append(f)

data_fitness=pd.concat(data,axis=0,keys=backgrounds)


# +
## Distribution of fitness effects
plt.figure(figsize=(8,5))
data_fitness_wt=data_fitness.loc[["wt_a","wt_b"]]
data_fitness_wt.loc[:,"fitness_gene"].hist(bins=100,color="gray",alpha=0.7)
#data_fitness_wt.loc[:,"fitness_domains_corrected"].astype(float).hist(bins=100)
#plt.title("Distribution of fitness effects in WT",fontsize=20)


plt.xlabel("Fitness effect",fontsize=14)
plt.ylabel("Number of genes",fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14);
plt.title("Distribution of fitness effects in wild type",fontsize=14)
plt.grid(linewidth=0.3)
plt.tight_layout()
plt.savefig("../figures/fig_fitness_effects_wt.png",dpi=300,transparent=True)

# +
## Distribution of fitness effects

data_fitness_wt=data_fitness.loc[["bem1-aid_a","bem1-aid_b"]]
data_fitness_wt.loc[:,"fitness_domains_corrected"].astype(float).hist(bins=100)
plt.title("Distribution of fitness effects in $\Delta$bem1")
plt.xlabel("Fitness effect")
plt.ylabel("Number of genes")
plt.tight_layout()
#plt.savefig("../figures/fig_fitness_effects_dbem1.png",dpi=300,transparent=True)

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
## Distribution of fitness effects

data_fitness_wt=data_fitness.loc[["dbem1dbem3_a","dbem1dbem3_b"]]
data_fitness_wt.loc[:,"fitness_domains_corrected"].astype(float).hist(bins=100)
plt.title("Distribution of fitness effects in $\Delta$bem1$\Delta$bem3")
plt.xlabel("Fitness effect")
plt.ylabel("Number of genes")
plt.tight_layout()
# -

## some fitness values 
data_fitness.loc["wt_a","CDC24"],data_fitness.loc["wt_a","NRP1"]
#data_fitness_wt.loc["BEM1","fitness_gene"]

standard_essentials


# +
datawtessntials_a=data_fitness.loc["wt_a"].copy()

datawtessntials_b=data_fitness.loc["wt_b"].copy()

datawtessntials_a.loc[:,"Essential"]=False
datawtessntials_b.loc[:,"Essential"]=False


for i in standard_essentials:
    if i in datawtessntials_a.index:
        datawtessntials_a.loc[i,"Essential"]=True
    if i in datawtessntials_b.index:
        datawtessntials_b.loc[i,"Essential"]=True



# -

x=datawtessntials_a.loc[datawtessntials_a.loc[:,"Essential"]==True,"fitness_domains_corrected"]
x_1=datawtessntials_a.loc[datawtessntials_a.loc[:,"Essential"]==True,"fitness_gene"]
x_2=datawtessntials_a.loc[datawtessntials_a.loc[:,"Essential"]==True,"fitness_domains_average"]
y=datawtessntials_b.loc[datawtessntials_b.loc[:,"Essential"]==True,"fitness_domains_corrected"]
y_1=datawtessntials_b.loc[datawtessntials_b.loc[:,"Essential"]==True,"fitness_gene"]
y_2=datawtessntials_b.loc[datawtessntials_b.loc[:,"Essential"]==True,"fitness_domains_average"]

# +
x_float=[]
for i in x:
    if type(i)!=np.float64:
        x_float.append(i[0])
    else:
        x_float.append((i))

y_float=[]
for i in y:
    if type(i)!=np.float64:
        y_float.append(i[0])
    else:
        y_float.append((i))
# -

plt.figure(figsize=(8,5))
#plt.hist(x_1,bins=100,alpha=0.5,label="fitness of the whole gene",color="pink");
plt.hist(x_float,bins=100,alpha=0.8,label="fitness corrected by domains",color="pink");
#plt.hist(y_float,bins=100,alpha=0.5,label="fitness domains correction b",color="red");
#plt.hist(y_1,bins=100,alpha=0.2,label="fitness_gene b",color="blue");
#plt.vlines(np.mean(x_1),0,40,linestyles="dashed",color="black",label="mean")
plt.vlines(np.mean(x_float),0,200,linestyles="dashed",color="black",label="mean")
# plt.vlines(np.mean(x_float)-np.std(x_float),0,200,linestyles="dashed",color="black",
# label="std",linewidth=0.5)
# plt.vlines(np.mean(x_float)+np.std(x_float),0,200,linestyles="dashed",color="black",
# linewidth=0.5)
plt.legend()
plt.title("Distribution of fitness effects in essential genes using the domain correction",fontsize=14)
plt.xlabel("Fitness effect",fontsize=14)
plt.ylabel("Number of genes",fontsize=14)
plt.grid(linewidth=0.3)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14);
plt.tight_layout()
#plt.savefig("../figures/fig_fitness_effects_wt_essential_gene.png",dpi=300,transparent=True)
plt.savefig("../figures/fig_fitness_effects_wt_essential_domains_correction.png",dpi=300,transparent=True)
#plt.savefig("../figures/fig_fitness_effects_wt_essential_gene_vs_domains.png",dpi=300,transparent=True)

np.std(x_float)

# +
x_3=datawtessntials_a.loc[datawtessntials_a.loc[:,"Essential"]==False,"fitness_domains_corrected"]
x_4=datawtessntials_a.loc[datawtessntials_a.loc[:,"Essential"]==False,"fitness_gene"]
x_5=datawtessntials_a.loc[datawtessntials_a.loc[:,"Essential"]==False,"fitness_domains_average"]

y_3=datawtessntials_b.loc[datawtessntials_b.loc[:,"Essential"]==False,"fitness_domains_corrected"]
y_4=datawtessntials_b.loc[datawtessntials_b.loc[:,"Essential"]==False,"fitness_gene"]
y_5=datawtessntials_b.loc[datawtessntials_b.loc[:,"Essential"]==False,"fitness_domains_average"]

# +
fig,axes=plt.subplots(1,2,figsize=(15,5))
plt.subplot(1,2,1)
plt.title("essentials")
plt.errorbar(["a-gene","b-gene"],[x_1.mean(),y_1.mean()],yerr=[x_1.std(),y_1.std()],
fmt="o",color="red",capsize=10)
plt.errorbar(["a-domains","b-domains"],[x_2.mean(),y_2.mean()],yerr=[x_2.std(),y_2.std()],
fmt="o",color="blue",capsize=10)
plt.errorbar(["a-corrected","b-corrected"],[x.mean(),y.mean()],yerr=[x.std(),y.std()],
fmt="o",color="black",capsize=10)

plt.ylim(0,1.2)
plt.grid("on")
plt.subplot(1,2,2)

plt.title("non essentials")

plt.errorbar(["a-gene","b-gene"],[x_4.mean(),y_4.mean()],yerr=[x_4.std(),y_4.std()],
fmt="o",color="red",capsize=10,label="Non essentials")
plt.errorbar(["a-domains","b-domains"],[x_5.mean(),y_5.mean()],yerr=[x_5.std(),y_5.std()],
fmt="o",color="blue",capsize=10,label="Non essentials")
plt.errorbar(["a-corrected","b-corrected"],[x_3.mean(),y_3.mean()],yerr=[x_3.std(),y_3.std()],
fmt="o",color="black",capsize=10,label="Non essentials")

plt.ylim(0,1.2)
plt.grid("on")
# -

a=data_fitness.loc["wt_a","fitness_gene"]
b=data_fitness.loc["wt_b","fitness_gene"]
c=np.minimum(len(a),len(b))


# +
from sklearn.linear_model import LinearRegression

#initiate linear regression model
model = LinearRegression()

tmp_0=data_fitness.loc["wt_a","fitness_gene"][0:c]
tmp_1=data_fitness.loc["wt_b","fitness_gene"][0:c]



X, y = tmp_0.values.reshape(-1,1), tmp_1.values.reshape(-1,1)

#fit regression model
model.fit(X, y)

#calculate R-squared of regression model
r_squared = model.score(X, y)

figure,ax=plt.subplots(1,1,figsize=(8,5))
plt.scatter(tmp_0,tmp_1,s=1)
ax.plot([0,1.7],[0,1.7],color="black",linestyle="--")
ax.text(1.3, 1.55, '$R^2=%.3f$' % (r_squared),fontsize=12)

# +
bem1_f_a=data_fitness.loc["wt_a","BEM1"]["fitness_gene"]
bem1_f_b=data_fitness.loc["wt_b","BEM1"]["fitness_gene"]
data=data_fitness.loc["wt_a"]

significance_threshold = 0.05 #set significance threshold
gi=defaultdict(dict)
ttest_tval_list = [np.nan]*2 #initialize list for storing t statistics
ttest_pval_list = [np.nan]*2 #initialize list for storing p-values
signif_thres_list = False #initialize boolean list for indicating datapoints with p-value above threshold
fc_list = [np.nan]*2
for gene in data.index :
    geneX=gene
    
    geneX_f_a=data_fitness.loc["wt_a",geneX]["fitness_gene"]
    if geneX in data_fitness.loc["wt_b"].index:
        geneX_f_b=data_fitness.loc["wt_b",geneX]["fitness_gene"]
        if geneX in data_fitness.loc["bem1-aid_a"].index and geneX in data_fitness.loc["bem1-aid_b"].index:
            geneXbem1_f_a=data_fitness.loc["bem1-aid_a",geneX]["fitness_gene"]
            geneXbem1_f_b=data_fitness.loc["bem1-aid_b",geneX]["fitness_gene"]
            
            variable_a_array=[geneXbem1_f_a,geneXbem1_f_b]
            variable_b_array=[geneX_f_a*bem1_f_a,geneX_f_b*bem1_f_b]
            
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



#gi_pd[gi_pd.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)
# -

gi_pd_fitness_gene.loc["BEM3"]

# +
from annotate_volcano import annotate_volcano   #import annotate_volcano function
volcano_df=gi_pd_fitness_gene
fig=annotate_volcano(volcano_df,[2.2,-0.5],[1.5,1.3],figure_title="test",variable="fitness_per_gene")


# -

gi_pd_significant=gi_pd_fitness_gene[gi_pd_fitness_gene.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)

gi_pd_significant

# +
## comparison with existing data

sl_bem1=pd.read_excel("../postprocessed-data/BEM1_genetic_interactions_filtered_by_synt.xlsx",skiprows=8)
neg_bem1=pd.read_excel("../postprocessed-data/BEM1_genetic_interactions_filtered_by_negat.xlsx",skiprows=8)
pos_bem1=pd.read_excel("../postprocessed-data/BEM1_genetic_interactions_filtered_by_pos.xlsx",skiprows=8)
# -

sl_bem1_genes=sl_bem1.loc[:,"Interactor.1"]
neg_bem1_genes=neg_bem1.loc[:,"Interactor.1"]
pos_bem1_genes=pos_bem1.loc[:,"Interactor.1"]


len(sl_bem1_genes.unique())+len(neg_bem1_genes.unique()),len(pos_bem1_genes.unique())

# +
neg_satay_signif=gi_pd_significant[gi_pd_significant.loc[:,"fold_change"]<0].loc[:,"gene_names"].index
neg_satay=gi_pd_fitness_gene[gi_pd_fitness_gene.loc[:,"fold_change"]<0].loc[:,"gene_names"].index

pos_satay_signif=gi_pd_significant[gi_pd_significant.loc[:,"fold_change"]>0].loc[:,"gene_names"].index
pos_satay=gi_pd_fitness_gene[gi_pd_fitness_gene.loc[:,"fold_change"]>0].loc[:,"gene_names"].index
# -

len(neg_satay)+len(pos_satay),len(neg_satay_signif)+len(pos_satay_signif)

# +
common_neg_signif=[]
for gene in neg_satay_signif:
    if gene in neg_bem1_genes.unique() or gene in sl_bem1_genes.unique():
        common_neg_signif.append(gene)

common_neg=[]
for gene in neg_satay:
    if gene in neg_bem1_genes.unique() or gene in sl_bem1_genes.unique():
        common_neg.append(gene)

common_pos_signif=[]
for gene in pos_satay_signif:
    if gene in pos_bem1_genes.unique() :
        common_pos_signif.append(gene)

common_pos=[]
for gene in pos_satay:
    if gene in pos_bem1_genes.unique() :
        common_pos.append(gene)
# -

len(common_neg)/len(neg_satay),len(common_neg_signif)/len(neg_satay_signif),len(common_pos)/len(pos_satay),len(common_pos_signif)/len(pos_satay_signif)

len(common_neg)/(len(sl_bem1_genes.unique())+len(neg_bem1_genes.unique())),len(common_neg_signif)/(len(sl_bem1_genes.unique())+len(neg_bem1_genes.unique())),len(common_pos)/len(pos_bem1_genes.unique()),len(common_pos_signif)/len(pos_bem1_genes.unique())

# +
## Compare only with constanzo data
neg_bem1_costanzo_genes=neg_bem1[neg_bem1.loc[:,"Note"]=="Costanzo M"]["Interactor.1"].reset_index(drop=True)
neg_bem1_costanzo_pvalue=neg_bem1[neg_bem1.loc[:,"Note"]=="Costanzo M"]["Source"].reset_index(drop=True)
neg_bem1_costanzo_score=neg_bem1[neg_bem1.loc[:,"Note"]=="Costanzo M"]["P-value"].reset_index(drop=True)

pos_bem1_costanzo_genes=pos_bem1[pos_bem1.loc[:,"Note"]=="Costanzo M"]["Interactor.1"].reset_index(drop=True)
pos_bem1_costanzo_pvalue=pos_bem1[pos_bem1.loc[:,"Note"]=="Costanzo M"]["Source"].reset_index(drop=True)
pos_bem1_costanzo_score=pos_bem1[pos_bem1.loc[:,"Note"]=="Costanzo M"]["P-value"].reset_index(drop=True)

# +
## Now use the fitness correction by domains to find GI 
bem1_f_a=data_fitness.loc["wt_a","BEM1"]["fitness_domains_average"] # to accomodate for zero value of the corrected one
bem1_f_b=data_fitness.loc["wt_b","BEM1"]["fitness_domains_average"]# to accomodate for zero value
data=data_fitness.loc["wt_a"]

significance_threshold = 0.05 #set significance threshold
gi=defaultdict(dict)
ttest_tval_list = [np.nan]*2 #initialize list for storing t statistics
ttest_pval_list = [np.nan]*2 #initialize list for storing p-values
signif_thres_list = False #initialize boolean list for indicating datapoints with p-value above threshold
fc_list = [np.nan]*2
for gene in data.index :
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

        if geneX in data_fitness.loc["bem1-aid_a"].index and geneX in data_fitness.loc["bem1-aid_b"].index:
            geneXbem1_f_a=data_fitness.loc["bem1-aid_a",geneX]["fitness_domains_corrected"]
            geneXbem1_f_b=data_fitness.loc["bem1-aid_b",geneX]["fitness_domains_corrected"]
            if type(geneXbem1_f_a)!=np.float64:
                geneXbem1_f_a=geneXbem1_f_a[0]
            if type(geneXbem1_f_b)!=np.float64:
                geneXbem1_f_b=geneXbem1_f_b[0]
            
            
            variable_a_array=[geneXbem1_f_a,geneXbem1_f_b]
            variable_b_array=[geneX_f_a*bem1_f_a,geneX_f_b*bem1_f_b]
            
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

## replace inf values with zeros

# gi_pd.loc[gi_pd.loc[:,"fold_change"]==-np.inf,"fold_change"]=np.min(gi_pd.loc[gi_pd.loc[:,"fold_change"]!=-np.inf,"fold_change"])
# gi_pd.loc[gi_pd.loc[:,"fold_change"]==np.inf,"fold_change"]=np.max(gi_pd.loc[gi_pd.loc[:,"fold_change"]!=np.inf,"fold_change"])
#gi_pd[gi_pd.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)
# -

gi_pd[gi_pd.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)

# +
gi_pd_significant=gi_pd[gi_pd.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)

neg_satay_signif=gi_pd_significant[gi_pd_significant.loc[:,"fold_change"]<0].loc[:,"gene_names"].index
neg_satay=gi_pd[gi_pd.loc[:,"fold_change"]<0].loc[:,"gene_names"].index

pos_satay_signif=gi_pd_significant[gi_pd_significant.loc[:,"fold_change"]>0].loc[:,"gene_names"].index
pos_satay=gi_pd[gi_pd.loc[:,"fold_change"]>0].loc[:,"gene_names"].index

# +
common_neg_signif=[]
for gene in neg_satay_signif:
    if gene in neg_bem1_costanzo_genes.unique() or gene in sl_bem1_genes.unique():
        common_neg_signif.append(gene)

common_neg=[]
for gene in neg_satay:
    if gene in neg_bem1_costanzo_genes.unique() or gene in sl_bem1_genes.unique():
        common_neg.append(gene)

common_pos_signif=[]
for gene in pos_satay_signif:
    if gene in pos_bem1_costanzo_genes.unique() :
        common_pos_signif.append(gene)

common_pos=[]
for gene in pos_satay:
    if gene in pos_bem1_costanzo_genes.unique() :
        common_pos.append(gene)

# +
consistent_neg_interactors_bem1=list(set(neg_satay) & set(neg_bem1_costanzo_genes) )

consistent_pos_interactors_bem1=list(set(pos_satay) & set(pos_bem1_costanzo_genes))

consistent_pos_interactors_bem1,consistent_neg_interactors_bem1

# +
from matplotlib_venn import venn2, venn2_circles

plt.figure(figsize=(10,5))
venn2(subsets = (len(pos_satay), len(pos_bem1_costanzo_genes), len(consistent_pos_interactors_bem1)), 
set_labels = ('SATAy', 'SGA'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(pos_satay), len(pos_bem1_costanzo_genes), len(consistent_pos_interactors_bem1)));
plt.title("Consistency in positive interactions")
#
#plt.savefig("../figures/venn_pos.png",dpi=300)

# +

plt.figure(figsize=(10,5))
venn2(subsets = (len(neg_satay), len(neg_bem1_costanzo_genes), len(consistent_neg_interactors_bem1)), 
set_labels = ('SATAy', 'SGA'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(neg_satay), len(neg_bem1_costanzo_genes), len(consistent_neg_interactors_bem1)));
plt.title("Consistency in negative interactions")
#plt.savefig("../figures/venn_neg.png",dpi=300)

# +
pos_satay_neg_costanzo=list(set(pos_satay) & set(neg_bem1_costanzo_genes) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(pos_satay), len(neg_bem1_costanzo_genes), len(pos_satay_neg_costanzo)), 
set_labels = ('PI SATAy', 'NI SGA'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(pos_satay), len(neg_bem1_costanzo_genes), len(pos_satay_neg_costanzo)));
plt.title("Misclasification I")

#plt.savefig("../figures/venn_misclass1.png",dpi=300)

# +
neg_satay_pos_costanzo=list(set(neg_satay) & set(pos_bem1_costanzo_genes) )
plt.figure(figsize=(10,5))
venn2(subsets = (len(neg_satay), len(pos_bem1_costanzo_genes), len(neg_satay_pos_costanzo)), 
set_labels = ('NI SATAy', 'PI SGA'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(neg_satay), len(pos_bem1_costanzo_genes), len(neg_satay_pos_costanzo)));
plt.title("Misclasification II")

#plt.savefig("../figures/venn_misclass2.png",dpi=300)
# -

## Plot the volcano plot
from annotate_volcano import annotate_volcano
volcano_df=gi_pd
fig=annotate_volcano(volcano_df,[4,-2.3],[2,1.5],figure_title="Interactors of BEM1 in WT",variable="fitness_corrected")
plt.savefig("../figures/volcano_bem1_wt.png",dpi=300,transparent=True)


# +
# import gseapy as gp
# from gseapy import barplot, dotplot
# type_gi="Neg_GI_satay"
# goi=neg_satay_signif.tolist()
# yeast = gp.get_library_name(organism='Yeast')

# sets=[yeast[2],yeast[5],yeast[8],yeast[11],yeast[16] ] 
# #['GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
# #'Gene_Interaction_Hubs_BioGRID_2018','Pfam_Domains_2019']
# # #%% enrichment 

# for i in np.arange(0,len(sets)): 

#   enr=gp.enrichr(gene_list=goi,
#                   gene_sets=sets[i],
#                   organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                
#                   outdir='../postprocessed-data/enrich-analysis_dbem1/'+type_gi+'/',
#                   # no_plot=True,
#                   cutoff=0.5 # test dataset, use lower value from range(0,1)
#                 )
      
# # to save your figure, make sure that ``ofname`` is not None
#   ax = dotplot(enr.res2d, title=sets[i],cmap='viridis_r', size=20, figsize=(3,5),ofname=sets[i]+type_gi)
  
    
# -

gi_pd[gi_pd.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=True)

geneX="BEM3"
fbem1bem3=data_fitness.loc["bem1-aid_a",geneX]["fitness_domains_corrected"]
fbem1bem3b=data_fitness.loc["bem1-aid_b",geneX]["fitness_domains_corrected"]

# +
bem1_f_a=data_fitness.loc["wt_a","BEM1"]["fitness_domains_average"] # to accomodate for zero value of the corrected one
bem1_f_b=data_fitness.loc["wt_b","BEM1"]["fitness_domains_average"]

bem1_f_a,bem1_f_b

# +
fgenex=data_fitness.loc["wt_a",geneX]["fitness_domains_corrected"] 
fgenexb=data_fitness.loc["wt_b",geneX]["fitness_domains_corrected"]


np.mean([fbem1bem3-fgenex*bem1_f_a,fbem1bem3b-fgenexb*bem1_f_b])
# -

backgrounds

# +
## to know if nrp1 deletion is beneficial in the dbem1dbem3 background


gene="NRP1"

f_bem13a=data_fitness.loc["dbem1dbem3_a"]
f_bem13b=data_fitness.loc["dbem1dbem3_b"]

f_bem1a=data_fitness.loc["bem1-aid_a"]
f_bem1b=data_fitness.loc["bem1-aid_b"]


f_bem3a=data_fitness.loc["dbem3_a"]
f_bem3b=data_fitness.loc["dbem3_b"]

f_wta=data_fitness.loc["wt_a"]
f_wtb=data_fitness.loc["wt_b"]


f_nrp1_b13a=f_bem13a.loc[gene,"fitness_domains_corrected"]

f_nrp1_b13b=f_bem13b.loc[gene,"fitness_domains_corrected"]

f_nrp1_b1a=f_bem1a.loc[gene,"fitness_domains_corrected"]
f_nrp1_b1b=f_bem1b.loc[gene,"fitness_domains_corrected"]

f_nrp1_b3a=f_bem3a.loc[gene,"fitness_domains_corrected"]
f_nrp1_b3b=f_bem3b.loc[gene,"fitness_domains_corrected"]

f_nrp1_wta=f_wta.loc[gene,"fitness_domains_corrected"]
f_nrp1_wtb=f_wtb.loc[gene,"fitness_domains_corrected"]


score_a=f_nrp1_b13a-f_nrp1_b1a*f_nrp1_b3a-f_nrp1_wta
score_b=f_nrp1_b13b-f_nrp1_b1b*f_nrp1_b3b-f_nrp1_wtb

score_a,score_b,np.mean([score_a,score_b])
# -

f_dbem1a=data_fitness.loc["bem1-aid_a",:]["fitness_domains_corrected"].astype(float) 
f_dbem1a.loc["OPT2"]

# +
## interactors in dbem1dbem3

## Now use the fitness correction by domains to find GI 
f_dbem1a=data_fitness.loc["bem1-aid_a",:]["fitness_domains_corrected"].astype(float) # to accomodate for zero value of the corrected one
f_dbem1b=data_fitness.loc["bem1-aid_b",:]["fitness_domains_corrected"].astype(float)  # to accomodate for zero value

f_dbem3a=data_fitness.loc["dbem3_a",:]["fitness_domains_corrected"].astype(float) # to accomodate for zero value of the corrected one
f_dbem3b=data_fitness.loc["dbem3_b",:]["fitness_domains_corrected"].astype(float) # to accomodate for zero value of the corrected one

f_dbem1dbem3a=data_fitness.loc["dbem1dbem3_a",:]["fitness_domains_corrected"].astype(float) # to accomodate for zero value of the corrected one
f_dbem1dbem3b=data_fitness.loc["dbem1dbem3_b",:]["fitness_domains_corrected"].astype(float) # to accomodate for zero value of the corrected one

gene_intersection=list(set(f_dbem1a.index)&set(f_dbem3a.index)&set(f_dbem1dbem3a.index)
                       &set(f_dbem1b.index)&set(f_dbem3b.index)&set(f_dbem1dbem3b.index))

significance_threshold = 0.05 #set significance threshold
gi=defaultdict(dict)
ttest_tval_list = [np.nan]*2 #initialize list for storing t statistics
ttest_pval_list = [np.nan]*2 #initialize list for storing p-values
signif_thres_list = False #initialize boolean list for indicating datapoints with p-value above threshold
fc_list = [np.nan]*2
for gene in gene_intersection :

            
    variable_a_array=[f_dbem1dbem3a.loc[gene],f_dbem1dbem3b.loc[gene]]
    variable_b_array=[f_dbem1a.loc[gene]*f_dbem3a.loc[gene],f_dbem1b.loc[gene]*f_dbem3b.loc[gene]]
    
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
gibem1bem3=pd.DataFrame.from_dict(gi,orient="index")

gibem1bem3["gene_names"]=gibem1bem3.index

data=gibem1bem3.dropna(how='all')
data=data[~data.isin([np.nan, np.inf, -np.inf]).any(1)]

data_signif=data[data["significance"]==True]
# -

data_signif.sort_values(by="fold_change",ascending=False)

## Plot the volcano plot
from annotate_volcano import annotate_volcano
volcano_df=data
fig=annotate_volcano(volcano_df,[4,-2.3],[2,1.5],figure_title="Interactors of BEM1 in bem3$\Delta$",variable="fitness_corrected")
#plt.savefig("../figures/volcano_bem1_dbem3.png",dpi=300,transparent=True)


# +
gi_pd_significant_bem13=data[data.loc[:,"significance"]==True].sort_values(by="fold_change",ascending=False)

neg_satay_signif_bem13=gi_pd_significant_bem13[gi_pd_significant_bem13.loc[:,"fold_change"]<0].loc[:,"gene_names"].index
neg_satay_bem13=data[data.loc[:,"fold_change"]<0].loc[:,"gene_names"].index

pos_satay_signif_bem13=gi_pd_significant_bem13[gi_pd_significant_bem13.loc[:,"fold_change"]>0].loc[:,"gene_names"].index
pos_satay_bem13=data[data.loc[:,"fold_change"]>0].loc[:,"gene_names"].index

# +
neg_consistent_interactors_bem1_andbem1bem3=list(set(neg_satay_signif)&set(neg_satay_signif_bem13))

pos_consistent_interactors_bem1_andbem1bem3=list(set(pos_satay_signif)&set(pos_satay_signif_bem13))

# +

plt.figure(figsize=(10,5))
venn2(subsets = (len(pos_satay_signif), len(pos_satay_signif_bem13), len(pos_consistent_interactors_bem1_andbem1bem3)), 
set_labels = ('PI BEM1 WT', 'PI BEM1 $\Delta$bem3'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = (len(pos_satay_signif), len(pos_satay_signif_bem13), len(pos_consistent_interactors_bem1_andbem1bem3)));
plt.title("How many PI are the same in BEM1 and BEM1$\Delta$bem3")

#plt.savefig("../figures/pos_venn_bem1_bem1bem3.png",dpi=300,transparent=True)

# +
plt.figure(figsize=(10,5))
a=(len(neg_satay_signif), len(neg_satay_signif_bem13), len(neg_consistent_interactors_bem1_andbem1bem3))
venn2(subsets = a, set_labels = ('NI BEM1 WT', 'NI BEM1 $\Delta$bem3'), set_colors=('purple', 'g'), alpha = 0.5);
venn2_circles(subsets = a);
plt.title("How many NI are the same in BEM1 and BEM1$\Delta$bem3")

#plt.savefig("../figures/venn_ni_bem1_bem1bem3.png",dpi=300,transparent=True)

# +
## Essential genes

# in WT threshold f<0.4 (from domains corrected)

data_fitness_wta=data_fitness.loc["wt_a","fitness_domains_corrected"]
data_fitness_wtb=data_fitness.loc["wt_b","fitness_domains_corrected"]

data_fitness_dbem1a=data_fitness.loc["bem1-aid_a","fitness_domains_corrected"]
data_fitness_dbem1b=data_fitness.loc["bem1-aid_b","fitness_domains_corrected"]

data_fitness_dbem1dbem3a=data_fitness.loc["dbem1dbem3_a","fitness_domains_corrected"]
data_fitness_dbem1dbem3b=data_fitness.loc["dbem1dbem3_b","fitness_domains_corrected"]

# in WT threshold f<0.4(mean)-0.28(std) (from domains corrected)
essential_genes_wta=data_fitness_wta[data_fitness_wta<0.12].index
essential_genes_wtb=data_fitness_wtb[data_fitness_wtb<0.12].index

essential_genes_dbem1a=data_fitness_dbem1a[data_fitness_dbem1a<0.12].index
essential_genes_dbem1b=data_fitness_dbem1b[data_fitness_dbem1b<0.12].index

essential_genes_dbem1dbem3a=data_fitness_dbem1dbem3a[data_fitness_dbem1dbem3a<0.12].index
essential_genes_dbem1dbem3b=data_fitness_dbem1dbem3b[data_fitness_dbem1dbem3b<0.12].index

# +
number_e_genes=[0.5*(len(essential_genes_wta)+len(essential_genes_wtb)),
0.5*(len(essential_genes_dbem1a)+len(essential_genes_dbem1b)),
0.5*(len(essential_genes_dbem1dbem3a)+len(essential_genes_dbem1dbem3b))]

std_number_e_genes=[np.std([len(essential_genes_wta),len(essential_genes_wtb)]),
np.std([len(essential_genes_dbem1a),len(essential_genes_dbem1b)]),
np.std([len(essential_genes_dbem1dbem3a),len(essential_genes_dbem1dbem3b)])]


# +
plt.figure(figsize=(5,5))

number_e_genes_refactor=[number_e_genes[1],number_e_genes[2],number_e_genes[0]]
std_number_e_genes_refactor=[std_number_e_genes[1],std_number_e_genes[2],std_number_e_genes[0]]
plt.errorbar(["BEM1$\Delta$","BEM1$\Delta$bem3","WT"],number_e_genes_refactor,
yerr=std_number_e_genes_refactor,fmt="o",color="black",capsize=10)
plt.ylabel("Predicted number of essential genes")

plt.plot(number_e_genes_refactor,"--",color="black",linewidth=0.5)
plt.tight_layout()
#plt.savefig("../figures/predicted_number_essential_genes_trajectory.png",dpi=300,transparent=True)

# +
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.errorbar(["WT","$\Delta$bem1","$\Delta$bem1$\Delta$bem3"],number_e_genes,
yerr=std_number_e_genes,fmt="o",color="black",capsize=10,markersize=15)
ax1.set_ylabel("Predicted number of essential genes")

ax1.plot(number_e_genes,"--",color="black",linewidth=0.5)

ax2.errorbar(["WT","$\Delta$bem1","$\Delta$bem1$\Delta$bem3"],[0.0082,0.001,0.006,],fmt="o",
color="orange",capsize=10,markersize=15)
ax2.plot([0.0082,0.001,0.006],"--",color="black",linewidth=0.5)

ax2.set_ylabel("Growth rate (min$^{-1}$)")
plt.tight_layout()

plt.savefig("../figures/predicted_number_essential_genes_trajectory_growth_rate.png",dpi=300,transparent=True) 
#plt.savefig("../figures/predicted_number_essential_genes_trajectory.png",dpi=300,transparent=True)

# +
## venn diagram common essential genes : are the esssential genes from bem1 also in WT? 

common_wt_bem1_e_genes=list(set(essential_genes_wta)&set(essential_genes_wtb)&set(essential_genes_dbem1a)&set(essential_genes_dbem1b))
diff_wt_bem1_e_genesa=np.setdiff1d(essential_genes_dbem1a,essential_genes_wta)
diff_wt_bem1_e_genesb=np.setdiff1d(essential_genes_dbem1b,essential_genes_wtb)
diff_wt_bem1_e_genes=np.unique(diff_wt_bem1_e_genesa.tolist()+diff_wt_bem1_e_genesb.tolist())

e_wt_and_no_dbem1=np.unique(np.setdiff1d(essential_genes_wta,essential_genes_dbem1a).tolist()+np.setdiff1d(essential_genes_wtb,essential_genes_dbem1b).tolist())
# diff_wt_bem1_e_genes=np.unique(diff_wt_bem1_e_genesa+diff_wt_bem1_e_genesb)


common_bem1_bem1bem3_e_genes=list(set(essential_genes_dbem1a)&set(essential_genes_dbem1b)&set(essential_genes_dbem1dbem3a)&set(essential_genes_dbem1dbem3b))
diff_bem1_bem1bem3_e_genes=np.unique(np.setdiff1d(essential_genes_dbem1dbem3a,essential_genes_dbem1a).tolist()+np.setdiff1d(essential_genes_dbem1dbem3b,essential_genes_dbem1b).tolist())


common_wt_bem1bem3_e_genes=list(set(essential_genes_wta)&set(essential_genes_wtb)&set(essential_genes_dbem1dbem3a)&set(essential_genes_dbem1dbem3b))
common_wt_bem1_bem1bem3_e_genes=list(set(essential_genes_wta)&set(essential_genes_wtb)&set(essential_genes_dbem1a)&set(essential_genes_dbem1b)&set(essential_genes_dbem1dbem3a)&set(essential_genes_dbem1dbem3b))
# -

len(e_wt_and_no_dbem1)

# +
plt.figure(figsize=(10,5))


a=(number_e_genes[0], number_e_genes[1], len(common_wt_bem1_e_genes))
venn2(subsets = a, set_labels = ('Essential WT', 'Essentials $\Delta$bem1'),
 set_colors=('gray', 'blue'), alpha = 0.5);
venn2_circles(subsets = a);
plt.title("How many essential genes are the same in WT and $\Delta$bem1")

plt.savefig("../figures/venn_essential_genes_wt_bem1.png",dpi=300,transparent=True)

# +
plt.figure(figsize=(10,5))


a=(number_e_genes[1], number_e_genes[2], len(common_bem1_bem1bem3_e_genes))
venn2(subsets = a, set_labels = ('Essentials $\Delta$bem1', 'Essentials $\Delta$bem1$\Delta$bem3'),
 set_colors=('blue', 'yellow'), alpha = 0.5);
venn2_circles(subsets = a);
plt.title("How many essential genes are the same in $\Delta$bem1 and $\Delta$bem1$\Delta$bem3")

plt.savefig("../figures/venn_essential_genes_bem1_bem1bem3.png",dpi=300,transparent=True)

# +
plt.figure(figsize=(10,5))


a=(number_e_genes[0], number_e_genes[2], len(common_wt_bem1bem3_e_genes))
venn2(subsets = a, set_labels = ('Essentials WT', 'Essentials $\Delta$bem1$\Delta$bem3'), 
set_colors=('gray', 'yellow'), alpha = 0.5);
venn2_circles(subsets = a);
plt.title("How many essential genes are the same in WT and $\Delta$bem1$\Delta$bem3")

plt.savefig("../figures/venn_essential_genes_wt_bem1bem3.png",dpi=300,transparent=True)

# +
from matplotlib_venn import venn3, venn3_circles

a=(np.round(number_e_genes[0],0).astype(int), np.round(number_e_genes[1],0).astype(int), len(common_wt_bem1_e_genes), number_e_genes[2], 
len(common_wt_bem1bem3_e_genes), len(common_bem1_bem1bem3_e_genes), len(common_wt_bem1_bem1bem3_e_genes))
venn3(subsets = a, set_labels = ('Essentials WT', 'Essentials $\Delta$bem1','Essentials $\Delta$bem1$\Delta$bem3'), alpha = 0.5,
set_colors=('gray', 'blue', 'yellow'));
#venn3_circles(subsets = a);

#plt.savefig("../figures/venn_essential_genes_wt_bem1_bem1bem3.png",dpi=300,transparent=True)

# +
# save to txt file

np.savetxt("../data/essential_genes_wt_dbem1_dbem13.txt",common_wt_bem1bem3_e_genes,fmt="%s")

# +
import gseapy as gp
from gseapy import barplot, dotplot
type_gi="essentials_wt_and no_in_dbem1"
goi=e_wt_and_no_dbem1.tolist()
yeast = gp.get_library_name(organism='Yeast')

sets=[yeast[2],yeast[5],yeast[8],yeast[11],yeast[16] ] 
#['GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
#'Gene_Interaction_Hubs_BioGRID_2018','Pfam_Domains_2019']
# #%% enrichment 

for i in np.arange(0,len(sets)): 

  enr=gp.enrichr(gene_list=goi,
                  gene_sets=sets[i],
                  organism='Yeast', # don't forget to set organism to the one you desired! e.g. Yeast
                
                  outdir='../postprocessed-data/enrich-analysis_essential_genes/'+type_gi+'/',
                  # no_plot=True,
                  cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
      
# to save your figure, make sure that ``ofname`` is not None
  ax = dotplot(enr.res2d, title=sets[i],cmap='viridis_r', size=20, figsize=(3,5),ofname=sets[i]+type_gi)
  
    

# +
# how many essential genes are also NI genes?

e_genes_neg_bem1=list(set(essential_genes_dbem1a)&set(essential_genes_dbem1b)&set(neg_satay_signif))

e_genes_neg_bem1bem3=list(set(essential_genes_dbem1dbem3a)&set(essential_genes_dbem1dbem3b)&set(neg_satay_signif_bem13))

# +
plt.figure(figsize=(10,5))


a=(number_e_genes[1], len(neg_satay_signif), len(e_genes_neg_bem1))
venn2(subsets = a, set_labels = ('Essentials $\Delta$bem1', 'NI $\Delta$bem1'),
 set_colors=('blue', 'purple'), alpha = 0.5);
venn2_circles(subsets = a);
plt.title("How many essential genes are in $\Delta$bem1 are also NI genes")

#plt.savefig("../figures/venn_essential_genes_bem1_bem1bem3.png",dpi=300,transparent=True)

# +
plt.figure(figsize=(10,5))


a=(number_e_genes[2], len(neg_satay_signif_bem13), len(e_genes_neg_bem1bem3))
venn2(subsets = a, set_labels = ('Essentials $\Delta$bem1$\Delta$bem3', 'NI $\Delta$bem1$\Delta$bem3'),
 set_colors=('yellow', 'purple'), alpha = 0.5);
venn2_circles(subsets = a);
plt.title("How many essential genes are in $\Delta$bem1$\Delta$bem3 are also NI genes")

#plt.savefig("../figures/venn_essential_genes_bem1_bem1bem3.png",dpi=300,transparent=True)
# -

neg_satay_signif_bem13

# +
plt.plot(f_dbem1dbem3a.loc[neg_satay_signif_bem13],label="dbem1dbem3a")
plt.plot(f_dbem1dbem3b.loc[neg_satay_signif_bem13],label="dbem1dbem3b")

plt.plot(f_dbem1a.loc[neg_satay_signif_bem13],label="dbem1a")
plt.plot(f_dbem1b.loc[neg_satay_signif_bem13],label="dbem1b")

plt.plot(f_dbem3a.loc[neg_satay_signif_bem13],label="dbem3a")
plt.plot(f_dbem3b.loc[neg_satay_signif_bem13],label="dbem3b")

plt.xticks("")
plt.legend()

# +


plt.hist(f_dbem1dbem3a.loc[pos_satay_signif_bem13],label="dbem1dbem3a")
plt.hist(f_dbem1dbem3b.loc[pos_satay_signif_bem13],label="dbem1dbem3b")

plt.hist(f_dbem1a.loc[pos_satay_signif_bem13],label="dbem1a")
plt.hist(f_dbem1b.loc[pos_satay_signif_bem13],label="dbem1b")

plt.hist(f_dbem3a.loc[pos_satay_signif_bem13],label="dbem3a")
plt.hist(f_dbem3b.loc[pos_satay_signif_bem13],label="dbem3b")

#plt.xticks("")
plt.legend()
# -

pos_bem1_e_bem1bem3=list(set(pos_satay_signif)&set(essential_genes_dbem1dbem3a)&set(essential_genes_dbem1dbem3b))

f_dbem1dbem3b.loc[pos_bem1_e_bem1bem3]

f_dbem1b.loc[pos_bem1_e_bem1bem3]

polarity_genes.index

backgrounds

tmp=data_fitness.loc["wt_b"]
tmp.loc["BEM1","fitness_domains_corrected"].astype(float)

# +
## make an array with the fitness of every polarity gene in WT, dbem1 and dbem1dbem3

polarity_genes_fitness=np.zeros((len(polarity_genes.index),6))

for back in backgrounds[0:6]:
    tmp=data_fitness.loc[back]
    for i in range(len(polarity_genes.index)):
        if polarity_genes.index[i] in tmp.index:
        
            polarity_genes_fitness[i,backgrounds.index(back)]=tmp.loc[polarity_genes.index[i],"fitness_domains_corrected"].astype(float)
        else:
            polarity_genes_fitness[i,backgrounds.index(back)]=np.nan



# +
array2heatmap=polarity_genes_fitness
chosen_polarity_genes=polarity_genes.index

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 15))
array2heatmap[~np.isfinite(array2heatmap)] = 0
sns.heatmap(array2heatmap,cmap="seismic",vmin=0,vmax=1,xticklabels=backgrounds[0:6],
yticklabels=chosen_polarity_genes ,cbar=True,annot=True,cbar_kws={'label': 'Fitness corrected'})

labels=[]
for i in np.arange(0,len(chosen_polarity_genes)):
    labels.append(chosen_polarity_genes[i]+"-"+str(polarity_genes.loc[chosen_polarity_genes[i],:].unique()))
ax.set_yticklabels(labels);

plt.savefig("../figures/heatmap_polarity_genes.png",dpi=300,transparent=True)

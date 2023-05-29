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

# +
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys
from collections import defaultdict
from ast import literal_eval
from scipy.stats import norm
from cmath import isfinite
from scipy.stats import pearsonr
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


from from_excel_to_list import from_excel_to_list

## standard for plots
plt.rc('font', family='serif',size=14)
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)



# +
#### compute the dataframe of differences between SATAY fitness and each of the rest####

def compute_fitness_diff(satay_dataset,other_techique_dataset,name_other_technique):

    fitness_difference=defaultdict(dict)
    for i in other_techique_dataset.index:
        if name_other_technique=="sga":
            if i.upper() in satay_dataset.index and satay_dataset.loc[i.upper(),"fitness"]!="Not enough flanking regions":
                if type(other_techique_dataset.loc[i,"fitness"])!=np.float64:
                    fitness_difference[i]["difference"]=(other_techique_dataset.loc[i,"fitness"].unique()[0]-satay_dataset.loc[i.upper(),"fitness"])
                    fitness_difference[i][name_other_technique]=other_techique_dataset.loc[i,"fitness"].unique()[0]
                    fitness_difference[i]["SATAY"]=satay_dataset.loc[i.upper(),"fitness"]
                else:
                    fitness_difference[i]["difference"]=(other_techique_dataset.loc[i,"fitness"]-satay_dataset.loc[i.upper(),"fitness"])
                    fitness_difference[i][name_other_technique]=other_techique_dataset.loc[i,"fitness"]
                    fitness_difference[i]["SATAY"]=satay_dataset.loc[i.upper(),"fitness"]
            else:
                fitness_difference[i]["difference"]=np.nan
                fitness_difference[i][name_other_technique]=np.nan
                fitness_difference[i]["SATAY"]=np.nan
        else:
            if i in satay_dataset.index and satay_dataset.loc[i,"fitness"]!="Not enough flanking regions":
                if type(other_techique_dataset.loc[i,"fitness"])!=np.float64:
                    fitness_difference[i]["difference"]=(other_techique_dataset.loc[i,"fitness"].unique()[0]-satay_dataset.loc[i.upper(),"fitness"])
                    fitness_difference[i][name_other_technique]=other_techique_dataset.loc[i,"fitness"].unique()[0]
                    fitness_difference[i]["SATAY"]=satay_dataset.loc[i,"fitness"]
                else:
                    fitness_difference[i]["difference"]=(other_techique_dataset.loc[i,"fitness"]-satay_dataset.loc[i.upper(),"fitness"])
                    fitness_difference[i][name_other_technique]=other_techique_dataset.loc[i,"fitness"]
                    fitness_difference[i]["SATAY"]=satay_dataset.loc[i,"fitness"]
            else:
                fitness_difference[i]["difference"]=np.nan
                fitness_difference[i][name_other_technique]=np.nan
                fitness_difference[i]["SATAY"]=np.nan


    fitness_difference=pd.DataFrame(fitness_difference).T
    fitness_difference=fitness_difference.dropna(how='all')
    fitness_difference=fitness_difference[~fitness_difference.isin([np.nan, np.inf, -np.inf]).any(1)]

    return fitness_difference

### Normalization function####

def normalize_data(data):
    data_norm=(data-data.min())/(data.max()-data.min())
    return data_norm

### Plotting function####

def plot_comparison_techniques(data_difference,data_satay,data_technique,name_technique,savefig=True,domain_correction=False):

    corr, p = pearsonr(data_technique, data_satay)
    data2fit=data_difference
    
    mu, std = norm.fit(data2fit)

    fig,ax=plt.subplots(1,2,figsize=(10,5))
    plt.subplots_adjust(wspace=0.3)
    if name_technique=="SGA":
        ax[0].hist(data_technique,bins=100,density=True,label=name_technique,
    color="orange",histtype="step");
    elif name_technique=="QIAN":
        ax[0].hist(data_technique,bins=100,density=True,label=name_technique,
    color="green",histtype="step");
    elif name_technique=="BRESLOW":
        ax[0].hist(data_technique,bins=100,density=True,label=name_technique,
    color="purple",histtype="step");
    ax[0].hist(data_satay,bins=100,density=True,label="SATAY",
    color="grey",histtype="step");

    # axins0 = ax[0].inset_axes([0.25, 0.26, 0.3, 0.3])
    
    # axins0.errorbar(["SGA"],data_technique.mean(),yerr=data_technique.std(),fmt='o',
    #             capsize=5,capthick=2,elinewidth=2,color="black");
    # axins0.errorbar(["SATAY"],data_satay.mean(),yerr=data_satay.std(),fmt='o',
    #                 capsize=5,capthick=2,elinewidth=2,color="gray");
    
    # axins0.set_ylim(0,1.2)
    # axins0.set_xticks([0,1])
    # axins0.set_xticklabels(["",""])
    # axins0.set_ylabel("$\mu \pm \sigma$")
    ax[0].set_xlabel("Normalized fitness")
    ax[0].set_ylabel("Density")
    ax[0].set_title("Fitness distributions")
    ax[0].legend(loc="upper left")



    ax[1].hist2d(data_technique,data_satay,bins=(100,100),vmin=0,
    vmax=1,density=True,cmap=plt.cm.Greys);

    fig.colorbar(ax[1].collections[0],ax=ax[1])

    ax[1].plot([0,1],[0,1],color="black",alpha=0.8,linestyle="--",linewidth=1)

    #ax[1].scatter(data.loc[:,"sga"],data.loc[:,"satay"],s=10,alpha=0.2,color="gray");
    #sns.jointplot(data.loc[:,"sga"],data.loc[:,"satay"],kind="scatter",s=10,alpha=0.2,ax=ax[1]);
    ax[1].set_xlabel("Normalized " + name_technique + " fitness")
    ax[1].set_ylabel("Normalized SATAY fitness")

    ax[1].set_xlim(0,1)
    ax[1].set_ylim(0,1)
    ax[1].text(0,0.8,"Pearson corr={:.2f}".format(corr),color="black")



    # ax[2].hist(data_difference,bins=100,alpha=0.8,density=True,color="gray");
    # xmin, xmax = ax[2].get_xlim()
    # x = np.linspace(xmin, xmax, 100)
    # p = norm.pdf(x, mu, std)
    # ax[2].plot(x, p, 'k', linewidth=2)
    # ax[2].text(-0.8, 1, r'$\mu=%.2f,\ \sigma=%.2f$' % (mu, std))

    # ax[2].set_xlabel("Norm Fitness "+ name_technique+" - Norm Fitness SATAY")
    # ax[2].set_ylabel("Density")
    # ax[2].set_xlim(-1,1)
        
    for axes in ax:
        
        axes.tick_params(bottom=True, top=True, left=True, right=True)
        axes.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
        

        
    plt.tight_layout()

    if savefig and domain_correction==True:
        print("Domain correction applied on fitness and saving figure...")
        fig.savefig("../figures/fig_corrected_Fitness_"+name_technique+"_vs_SATAY.png",dpi=300,transparent=True)
    if savefig and domain_correction==False:
        fig.savefig("../figures/fig_Fitness_"+name_technique+"_vs_SATAY.png",dpi=300,transparent=True)
        print("Saving figure...")
    
# -

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
## Fitness datasets to compare with SATAY

#### SGA #####

sga_fitness=pd.read_excel('../data/fitness-SGD.xlsx',index_col=0)

sga_fitness=sga_fitness.loc[:,"fitness"]
sga_fitness=pd.DataFrame(sga_fitness,columns=["fitness"])

#### Qian #####
conversion=pd.read_csv('../postprocessed-data/from_systematic_genenames2standard_genenames.csv',delimiter=',')
conversion.columns=['systematic_name','standard_name']

qian_dataset=pd.read_excel('../postprocessed-data/fitness_Qian.xlsx',
header=1, engine='openpyxl')



#%% Replacing the systematic names from QIan to standard names
qian_dataset_SC=qian_dataset.loc[:,['ORF','SC fitness']]
j=0
for i in qian_dataset_SC['ORF']:
    value=conversion[conversion.loc[:,'systematic_name']==i]['standard_name']
    if len(value)!=0 :
        qian_dataset_SC.loc[j,'ORF_standard']=value.tolist()[0]
    j=j+1
#%%
qian_dataset_SC.columns=["ORF","fitness","Gene"]

qian_dataset_SC.index=qian_dataset_SC.loc[:,'Gene']
qian_dataset_SC.drop(columns=["ORF","Gene"],inplace=True)
qian_dataset_SC.dropna(inplace=True)


##### Breslow #####
breslow_dataset=pd.read_excel('../postprocessed-data/fitness_breslow.xlsx',
header=0, engine='openpyxl',sheet_name="Deletion growth rates")
breslow_dataset=breslow_dataset.loc[:,['Gene','Median Growth Rate']]
breslow_dataset.index=breslow_dataset['Gene']
breslow_dataset.drop(columns=["Gene"],inplace=True)
breslow_dataset.columns=['fitness']

# +
#######  Import fitness from SATAY #########
import pickle
with open("../postprocessed-data/fitness_models_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

satay_fitness=pd.concat(b,axis=0,keys=keys)

standard_essentials=np.loadtxt("../postprocessed-data/standard_essentials.txt",dtype=str)


# -

def get_the_right_dataframe(satay_fitness_dataset,background,type="average",essential_genes=standard_essentials):
    
    x=satay_fitness_dataset.loc[background]

    satay_wt2compare=x.loc[~x.index.isin(standard_essentials)] ## only use non essential genes 
    #satay_wt2compare=satay_wt2compare.loc[:,"fitness_domains_corrected"]
    if type == "average":
        satay_wt2compare=satay_wt2compare.loc[:,"fitness_gene"]
    if type == "domains":
        satay_wt2compare=satay_wt2compare.loc[:,"fitness_domains_corrected"]


    satay_wt2compare=pd.DataFrame(satay_wt2compare)
    satay_wt2compare.columns=["fitness"]

    satay_wt2compare.index.name="Gene"
    satay_wt2compare.dropna(inplace=True)

    satay_wt2compare=satay_wt2compare[satay_wt2compare.loc[:,"fitness"]!= "Not enough flanking regions"]

    satay_wt2compare=satay_wt2compare[satay_wt2compare.loc[:,"fitness"]!= "Not enough reads"]
    satay_wt2compare=satay_wt2compare[satay_wt2compare.loc[:,"fitness"]!= "Not enough insertions"]

    return satay_wt2compare

satay_wt2compare_average=get_the_right_dataframe(satay_fitness,"wt_merged",type="average",essential_genes=standard_essentials)
satay_wt2compare_domains=get_the_right_dataframe(satay_fitness,"wt_merged",type="domains",essential_genes=standard_essentials)

# +
# satay_wt2compare_average_norm=satay_wt2compare_average.apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))

# qian_dataset_SC_norm=qian_dataset_SC.apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
# breslow_dataset_norm=breslow_dataset.apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
# sga_fitness_norm=sga_fitness.apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))

# satay_wt2compare_domains_norm=satay_wt2compare_domains.apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
# -

# ### Scaling the data using the robust scaler to avoid outliers 

# +
# from sklearn.preprocessing import RobustScaler

# scaler = RobustScaler().fit(satay_wt2compare)
# satay_wt2compare_scaled = scaler.transform(satay_wt2compare)

# scaler = RobustScaler().fit(breslow_dataset)
# breslow_dataset_scaled = scaler.transform(breslow_dataset)

# scaler = RobustScaler().fit(qian_dataset_SC)
# qian_dataset_scaled = scaler.transform(qian_dataset_SC)


# scaler = RobustScaler().fit(sga_fitness)
# sga_dataset_scaled = scaler.transform(sga_fitness)

# #%% convert to  dataframe the scalers

# satay_wt2compare_scaled=pd.DataFrame(satay_wt2compare_scaled,index=satay_wt2compare.index,columns=["fitness"])

# breslow_dataset_scaled=pd.DataFrame(breslow_dataset_scaled,index=breslow_dataset.index,columns=["fitness"])

# qian_dataset_scaled=pd.DataFrame(qian_dataset_scaled,index=qian_dataset_SC.index,columns=["fitness"])

# sga_dataset_scaled=pd.DataFrame(sga_dataset_scaled,index=sga_fitness.index,columns=["fitness"])




# +
array_sga_index=sga_fitness.index

array_sga_index_upper=[i.upper() for i in array_sga_index]

sga_fitness.index=array_sga_index_upper



# +
## take the intercep of all index 
all_index_arrays=np.array([satay_wt2compare_average.index,breslow_dataset.index,
qian_dataset_SC.index,sga_fitness.index,satay_wt2compare_domains.index])

d=set.intersection(*map(set,all_index_arrays))

# satay_wt2compare_average_normalize=satay_wt2compare_average.loc[d].apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
# breslow_dataset_normalize=breslow_dataset.loc[d].apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
# qian_dataset_normalize=qian_dataset_SC.loc[d].apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
# sga_dataset_normalize=sga_fitness.loc[d].apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
# satay_wt2compare_domains_normalize=satay_wt2compare_domains.loc[d].apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))

satay_wt2compare_average_normalize=satay_wt2compare_average.loc[d].apply(lambda x: (x - np.mean(x)) / (np.std(x)))
breslow_dataset_normalize=breslow_dataset.loc[d].apply(lambda x: (x - np.mean(x)) / (np.std(x)))
qian_dataset_normalize=qian_dataset_SC.loc[d].apply(lambda x: (x - np.mean(x)) / (np.std(x)))
sga_dataset_normalize=sga_fitness.loc[d].apply(lambda x: (x - np.mean(x)) / (np.std(x)))
satay_wt2compare_domains_normalize=satay_wt2compare_domains.loc[d].apply(lambda x: (x - np.mean(x)) / (np.std(x)))

## remove index that are duplicates 

sga_dataset_normalize=sga_dataset_normalize.loc[~sga_dataset_normalize.index.duplicated(keep='first')]


# +
## Remove nan and Inf values

satay_average=satay_wt2compare_average_normalize.fitness/satay_wt2compare_average_normalize.fitness.median()
sga=sga_dataset_normalize.fitness/sga_dataset_normalize.fitness.median()
qian=qian_dataset_normalize.fitness/qian_dataset_normalize.fitness.median()
breslow=breslow_dataset_normalize.fitness/breslow_dataset_normalize.fitness.median()
satay_domains=satay_wt2compare_domains_normalize.fitness/satay_wt2compare_domains_normalize.fitness.median()


satay_average=satay_average[satay_average.apply(lambda x: isfinite(x))]
sga=sga[sga.apply(lambda x: isfinite(x))]
qian=qian[qian.apply(lambda x: isfinite(x))]
breslow=breslow[breslow.apply(lambda x: isfinite(x))]
satay_domains=satay_domains[satay_domains.apply(lambda x: isfinite(x))]

## select the common index_intersection

index_arrays=np.array([satay_average.index,sga.index,qian.index,breslow.index,satay_domains.index])

d1=set.intersection(*map(set,index_arrays))

satay_average=satay_average.loc[d1]
sga=sga.loc[d1]
qian=qian.loc[d1]
breslow=breslow.loc[d1]
satay_domains=satay_domains.loc[d1]


# +
corr_sga, p = pearsonr(sga,satay_average)
corr_qian,p=pearsonr(qian,satay_average)  
corr_breslow,p=pearsonr(breslow,satay_average)

corr_sga_2,p=pearsonr(sga,satay_domains)
corr_qian_2,p=pearsonr(qian,satay_domains)
corr_breslow_2,p=pearsonr(breslow,satay_domains)

corr_satay,p=pearsonr(satay_average,satay_domains)

print("Pearson correlation coefficient SGA vs SATAY_average: ",corr_sga)
print("Pearson correlation coefficient Qian vs SATAY_average: ",corr_qian)
print("Pearson correlation coefficient Breslow vs SATAY_average: ",corr_breslow)


print("Pearson correlation coefficient SGA vs Qian: ",pearsonr(sga,qian)[0])
print("Pearson correlation coefficient SGA vs Breslow: ",pearsonr(sga,breslow)[0])
print("Pearson correlation coefficient Qian vs Breslow: ",pearsonr(qian,breslow)[0])

print("Pearson correlation coefficient SGA vs SATAY_domains: ",corr_sga_2)
print("Pearson correlation coefficient Qian vs SATAY_domains: ",corr_qian_2)
print("Pearson correlation coefficient Breslow vs SATAY_domains: ",corr_breslow_2)

print("Pearson correlation coefficient SATAY_average vs SATAY_domains: ",corr_satay)

# +
plt.figure(figsize=(5,5))

plt.scatter(satay_average,sga,marker='o',color='orange',label="SGA",s=2)
plt.scatter(satay_average,qian,marker='o',color='g',label="Qian",s=2)
plt.scatter(satay_average,breslow,marker='o',color='purple',label="Breslow",s=2)

#
plt.legend()


plt.xlabel("SATAY fitness scaled")
plt.ylabel("Scaled fitness from other studies")

plt.tight_layout()

#plt.savefig("../figures/fig_fitness_comparison_in_one_plot.png",dpi=300)

# +
plt.subplots(figsize=(5,5))
satay_fitness_wt=satay_fitness.loc["wt_merged"]
plt.hist(satay_fitness_wt.loc[d1,"fitness_gene"],bins=100,color="k",histtype="step",linewidth=2);
#plt.hist(satay_wt.loc[d1,"fitness_domains_corrected"].astype(float),bins=100,histtype="step",linewidth=2,color="k");

sga_fitness.index=sga_fitness.index.str.upper()
plt.hist(sga_fitness.loc[d1,"fitness"],bins=100,color="orange",alpha=0.5,histtype="step",linewidth=2);
plt.hist(breslow_dataset.loc[d1,"fitness"],bins=100,color="purple",alpha=0.5,histtype="step",linewidth=2);
plt.hist(qian_dataset_SC.loc[d1,"fitness"],bins=100,color="green",alpha=0.5,    histtype="step",linewidth=2);



plt.ylabel("Number of genes")
plt.xlabel("Fitness")



plt.legend(["SATAY","SGA","Breslow","Qian"],loc="upper left")

#plt.title("SATAY fitness correction")

plt.tight_layout()

plt.savefig("../figures/fig_fitness_comparison_in_one_plot_hist.png",dpi=300)


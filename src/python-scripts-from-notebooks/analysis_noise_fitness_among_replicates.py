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
#     display_name: 'Python 3.9.7 64-bit (''transposonmapper'': conda)'
#     language: python
#     name: python3
# ---

# +
## Importing the required python libraries 
import os, sys
import warnings
import timeit
import numpy as np
import pandas as pd 
import pkg_resources
import matplotlib.pyplot as plt
import re
import seaborn as sns
from collections import defaultdict

from from_excel_to_list import from_excel_to_list
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression

plt.rc('font', family='serif',size=14)
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

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
## Importing fitness values 

import pickle
with open("../postprocessed-data/fitness_models_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

fitness_all_pd=pd.concat(b,axis=0,keys=keys)


# +
backgrounds=["wt_a","wt_b","bem1-aid_a","bem1-aid_b","dbem1dbem3_a","dbem1dbem3_b","dbem3_a","dbem3_b","dnrp1_1","dnrp1_2"]
fitness_backg=fitness_all_pd.loc[backgrounds]
data=[]
for keys in backgrounds:
    f=fitness_backg.loc[keys]
    f=f[f.loc[:,"fitness_gene"]!="Not enough flanking regions"]
    f=f[f.loc[:,"fitness_gene"]!="Not enough reads"]
    f=f[f.loc[:,"fitness_gene"]!="Not enough insertions"]
    f=f[f.loc[:,"fitness_domains_corrected"]!="Not enough insertions"]
    data.append(f)

data_fitness=pd.concat(data,axis=0,keys=backgrounds)


# -

data_fitness.columns

# +
replicate_1=data_fitness.loc["wt_a","fitness_gene"]
replicate_2=data_fitness.loc["wt_b","fitness_gene"]

## take the same index from both replicates
replicate_1=replicate_1[replicate_1.index.isin(replicate_2.index)]
replicate_2=replicate_2[replicate_2.index.isin(replicate_1.index)]


# +


#initiate linear regression model
model = LinearRegression()
## Fit a regression model

X, y = replicate_1.values.reshape(-1,1), replicate_2.values.reshape(-1,1)

#fit regression model
model.fit(X, y)

#calculate R-squared of regression model
r_squared = model.score(X, y)

# +
plt.scatter(replicate_1,replicate_2,c="black",s=1)

plt.xlabel("tech replicate 1")
plt.ylabel("tech replicate 2")

plt.plot([0,1.5],[0,1.5],color="black",linestyle="--")
plt.text(0, 1, '$R^2=%.3f$' % (r_squared),fontsize=12)
plt.grid(True,linewidth=0.2)
plt.title("Fitness whole gene")



# +
replicate_1_gene=data_fitness.loc["wt_a","fitness_domains_corrected"]
replicate_2_gene=data_fitness.loc["wt_b","fitness_domains_corrected"]

## take the same index from both replicates
replicate_1_gene=replicate_1_gene[replicate_1_gene.index.isin(replicate_2_gene.index)]
replicate_2_gene=replicate_2_gene[replicate_2_gene.index.isin(replicate_1_gene.index)]



#initiate linear regression model
model = LinearRegression()
## Fit a regression model

X, y = replicate_1_gene.values.reshape(-1,1), replicate_2_gene.values.reshape(-1,1)

#fit regression model
model.fit(X, y)

#calculate R-squared of regression model
r_squared_gene = model.score(X, y)

# +

figure,ax=plt.subplots(1,2,figsize=(10,5))

plt.subplots_adjust(wspace=0.4)
#plt.suptitle("Fitness values among technical replicates",fontsize=18)

ax[0].scatter(replicate_1,replicate_2,c="black",s=1)

ax[0].set_xlabel("tech replicate 1")
ax[0].set_ylabel("tech replicate 2")

ax[0].plot([0,2],[0,2],color="black",linestyle="--")
ax[0].text(0.1, 1.7, '$R^2=%.3f$' % (r_squared),fontsize=12)
ax[0].grid(True,linewidth=0.2)
ax[0].set_title("Non corrected fitness by domains")




plt.scatter(replicate_1_gene,replicate_2_gene,c="black",s=1)

ax[1].set_xlabel("tech replicate 1")
ax[1].set_ylabel("tech replicate 2")

ax[1].plot([0,2],[0,2],color="black",linestyle="--")
ax[1].text(0.1, 1.7, '$R^2=%.3f$' % (r_squared_gene),fontsize=12)
ax[1].grid(True,linewidth=0.2)
ax[1].set_title("Domain corrected fitness")

for ax in ax:
    ax.set_xlim(0,2)
    ax.set_ylim(0,2)
    ax.set_aspect('equal', 'box')

plt.tight_layout()
#plt.savefig("../figures/fig_fitness_replicates.png",dpi=300)
# -

from functions_analysis_frompergene2fitness import reads_per_insertion_along_gene_length
r_replicates=[]
replicates_backg=['wt_b',  'wt_a',]
for background in replicates_backg:
    
    r,gene_coordinates,reads_location,insertion_locations=reads_per_insertion_along_gene_length(list_data_pd,background,number_of_parts=10)
    
    r_replicates.append(r)

# +
## analyzing the reads per insertion per gene across replicates 

from sklearn.linear_model import LinearRegression
gene_parts=np.linspace(0,1,11)

r_sum_0=np.sum(r_replicates[0],axis=1)
r_sum_1=np.sum(r_replicates[1],axis=1)


#initiate linear regression model
model = LinearRegression()
X, y = r_sum_0.reshape(-1,1), r_sum_1.reshape(-1,1)

#fit regression model
model.fit(X, y)

#calculate R-squared of regression model
r_squared_sum = model.score(X, y)

ax=plt.figure(figsize=(5,5))
plt.scatter(r_sum_0,r_sum_1,s=1,color="k")
plt.yscale("log")
plt.xscale("log")
plt.plot([0,100000],[0,100000],color="k",linestyle="--")
plt.text(10, 1000, '$R^2=%.3f$' % (r_squared_sum),fontsize=12)

plt.xlabel("technical replicate 1",fontsize=12)
plt.ylabel("technical replicate 2",fontsize=12)

plt.title("Reads per insertions along gene length",fontsize=12)
plt.tight_layout()
plt.savefig("../figures/technical_replicates_reads_per_insertion.png",dpi=300)

# +
## select 100 numbers from 0 to 6599

random_numbers=np.random.randint(0,6599,10)

# +
ax=plt.figure(figsize=(5,5))

colors=plt.inferno()


for i in random_numbers:
    plt.scatter(r_replicates[0][i],r_replicates[1][i],s=50,color=colors)
    

#initiate linear regression model
model = LinearRegression()
r_0=r_replicates[0][random_numbers]
r_1=r_replicates[1][random_numbers]
X, y = r_0.reshape(-1,1), r_1.reshape(-1,1)

#fit regression model
model.fit(X, y)

#calculate R-squared of regression model
r_squared = model.score(X, y)

plt.yscale("log")
plt.xscale("log")
plt.ylabel("technical replicate 2",fontsize=12)
plt.xlabel("technical replicate 1",fontsize=12)
plt.title("Reads per insertions along gene length",fontsize=12)

plt.plot([0,100],[0,100],color="k",linestyle="--")
plt.text(0.1, 50, '$R^2=%.3f$' % (r_squared),fontsize=12)
plt.tight_layout()
#plt.savefig("../figures/technical_replicates_reads_per_insertion_10genes.png",dpi=300)

# +
## Make a subplot of the previous plots 

fig,ax=plt.subplots(1,2,figsize=(10,5))
plt.subplots_adjust(wspace=0.4)

plt.subplot(1,2,1)
plt.scatter(r_sum_0,r_sum_1,s=1,color="k")
plt.yscale("log")
plt.xscale("log")
plt.plot([0,100000],[0,100000],color="k",linestyle="--")
plt.text(10, 1000, '$R^2=%.3f$' % (r_squared_sum),fontsize=12)

plt.xlabel("technical replicate 1",fontsize=12)
plt.ylabel("technical replicate 2",fontsize=12)

#plt.title("Reads per insertions along gene length",fontsize=12)


plt.subplot(1,2,2)

colors=plt.inferno()


for i in random_numbers:
    plt.scatter(r_replicates[0][i],r_replicates[1][i],s=50,color=colors)

plt.yscale("log")
plt.xscale("log")
plt.ylabel("technical replicate 2",fontsize=12)
plt.xlabel("technical replicate 1",fontsize=12)
#plt.title("Reads per insertions along gene length",fontsize=12)

plt.plot([0,100],[0,100],color="k",linestyle="--")
plt.text(0.1, 50, '$R^2=%.3f$' % (r_squared),fontsize=12)

plt.savefig("../figures/technical_replicates_reads_per_insertion_subplots.png",dpi=300)

# +
## Exclude outliers of fitness .values()


IQR_replicate_1=np.percentile(replicate_1,75)-np.percentile(replicate_1,25)
IQR_replicate_2=np.percentile(replicate_2,75)-np.percentile(replicate_2,25)


replicate_1=replicate_1[replicate_1.values>1.5*IQR_replicate_1]
replicate_2=replicate_2[replicate_2.values>1.5*IQR_replicate_2]




## take the same index from both replicates
replicate_1=replicate_1[replicate_1.index.isin(replicate_2.index)]
replicate_2=replicate_2[replicate_2.index.isin(replicate_1.index)]

## initiate regression model
model = LinearRegression()
## Fit a regression model

X, y = replicate_1.values.reshape(-1,1), replicate_2.values.reshape(-1,1)

#fit regression model
model.fit(X, y)

#calculate R-squared of regression model
r_squared = model.score(X, y)

plt.scatter(replicate_1,replicate_2,c="black",s=1)
plt.text(0.1, 1.3, '$R^2=%.3f$' % (r_squared),fontsize=12)


# -

## Plot the relation between the mean and the std of the reads per insertion location
reads_mean=[]
reads_std=[]
a=list_data_pd.loc["wt_merged"]
reads=[]
insertions=[]
for i in np.arange(0,len(a)):
    
    tmp=from_excel_to_list(a.loc[i,"Reads per insertion location"])
    reads_mean.append(np.mean(tmp))
    reads_std.append(np.std(tmp))
    reads.append(tmp)
    if type(tmp)!=int:
        insertions.append(len(tmp)+1)
    else:
        insertions.append( 1)

plt.hist(insertions,bins=100);
1/np.mean(insertions)**2

plt.hist(reads_mean,bins=100);

# +
reads_mean=np.array(reads_mean)
reads_std=np.array(reads_std)
reads=np.array(reads)


Y=np.zeros_like(reads_mean)
X=np.zeros_like(reads_mean)

for i in np.arange(0,len(reads_mean)):
    Y[i]= (np.linalg.norm(reads[i]-reads_mean[i]))**2/np.array(reads_mean[i])-1
    X[i]=  reads_mean[i]


# +
## Use OLS to fit a linear regression model
import statsmodels.formula.api as smf

data = pd.DataFrame({'independent_var': X, 'dependient_var': Y})

model = smf.ols(formula='dependient_var ~ independent_var-1', data=data).fit()

## Print the model
print(model.summary())


# +
# Plotting
alpha = model.params

c = np.linspace(min(reads_mean),max(reads_mean),200)
d = c+alpha[0]*c**2 
# d=c+c-40

fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(111)
ax1.scatter(reads_mean,reads_std,alpha=0.4,s=1,c="k",label="Data");
ax1.plot(c,d,'--r',label='Neg. Binomial');
ax1.plot(c,c,color='limegreen',linestyle='--',label='Poisson')

# ax1.set_ylim(0,1000)
# ax1.set_xlim(0,1000)

ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.legend()
plt.xlabel('Mean')
plt.ylabel('Variance')
plt.tight_layout()
plt.savefig("../figures/fig_variance_mean_poisson_distribution_reads.png",dpi=300)

# +
plt.figure(figsize=(10,5))

plt.scatter(reads_mean,reads_std,c="black",s=1)


x_cords = range(0,10000)
y_cords= [(1.5*x**2) for x in x_cords]
p = np.poly1d( np.polyfit(reads_mean, reads_std, 1) )
plt.plot(x_cords, p(x_cords), color="green",linestyle="--")

plt.plot([0,10000],[0,10000],color="red",linestyle="--",label="$\sigma=\mu$")
#plt.plot(x_cords, y_cords, color="green",linestyle="--",label="$\sigma=\mu+0.5\mu$")
#plt.plot(x_cords, p(x_cords), color="green",linestyle="--")
plt.grid()
plt.xlabel("mean of reads per gene")
plt.ylabel("std of reads per gene")



plt.legend()
plt.yscale("log")
plt.xscale("log")
plt.xlim(0.1,10000)
plt.ylim(0.1,10000)

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

from from_excel_to_list import from_excel_to_list
from transposonmapper.statistics import volcano

from scipy import stats

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


# +
import pickle
with open("../postprocessed-data/fitness_models_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

fitness_all_pd=pd.concat(b,axis=0,keys=keys)

# +
data=[]
for backg in keys:
    f=fitness_all_pd.loc[backg]
    f=f[f.loc[:,"fitness_gene"]!="Not enough flanking regions"]
    for i in f.index:
        if f.loc[i,"fitness_gene"]=="Not enough reads":
            f.loc[i,"fitness_gene"]=0
        elif f.loc[i,"fitness_gene"]=="Not enough insertions":
            f.loc[i,"fitness_gene"]=0
        elif f.loc[i,"fitness_domains_corrected"]=="Not enough insertions":
            f.loc[i,"fitness_domains_corrected"]=0

    # f[f.loc[:,"fitness_gene"]=="Not enough reads"]=0
    # f[f.loc[:,"fitness_gene"]=="Not enough insertions"]=0
    # f[f.loc[:,"fitness_domains_corrected"]=="Not enough insertions"]=0
    
    # f=f[f.loc[:,"fitness_gene"]!="Not enough reads"]
    # f=f[f.loc[:,"fitness_gene"]!="Not enough insertions"]
    # f=f[f.loc[:,"fitness_domains_corrected"]!="Not enough insertions"]
    data.append(f)

data_fitness=pd.concat(data,axis=0,keys=keys)

# -

keys

# +
fitness_dbem1=data_fitness.loc["bem1-aid_merged"]
fitness_wt=data_fitness.loc["wt_merged"]
fitness_dbem1dbem3=data_fitness.loc["bem1-aid-dbem3_a"]

genes2test=["BEM3","CLA4","RDI1","STE18","SEC8","NRP1","BEM2"]
column="fitness_gene"


fitness_dbem12test=[]
fitness_dbem1dbem32test=[]
for gene in genes2test:
    fitness_dbem12test.append(fitness_dbem1.loc[gene,column])
    fitness_dbem1dbem32test.append(fitness_dbem1dbem3.loc[gene,column])

fitness_bem1=fitness_wt.loc["BEM1",column]/np.max(fitness_wt.loc[:,column])
fitness_wt2test=fitness_wt.loc[genes2test,column]/np.max(fitness_wt.loc[:,column])

fitness_dbem12test=fitness_dbem12test/np.max(fitness_dbem1.loc[:,column])
fitness_dbem1dbem32test=fitness_dbem1dbem32test/np.max(fitness_dbem1dbem3.loc[:,column])

# -

# ## The Wright-Fisher model to simulate evolution 
#
# The Wright-Fisher model is a mathematical model in population genetics that describes the stochastic process of genetic drift in a finite population. It was proposed independently by Sewall Wright and Ronald Fisher in the early 20th century.
#
# The model assumes a constant population size and a non-overlapping generation structure, meaning that each generation is formed by a random sample of individuals from the previous generation. In each generation, the model simulates the process of genetic drift, which is the random fluctuation of allele frequencies due to sampling error. Specifically, the model assumes that the allele frequency in the next generation is determined by a binomial sampling process from the allele frequency in the current generation.
#
# The Wright-Fisher model is widely used in population genetics to study the dynamics of allele frequency change in a population, under various conditions such as neutral evolution, selection, migration, and mutation. It is also used to study the effects of demographic events such as population bottlenecks, founder effects, and population expansions. The model provides a baseline for understanding the expected patterns of genetic variation in natural populations, and it serves as a useful null model for testing hypotheses about the causes of genetic variation.
#
# Wright-Fisher model makes several simplifying assumptions in order to predict evolution. These assumptions are:
#
# - Non-overlapping generations: The model assumes that generations do not overlap, meaning that individuals in a given generation do not reproduce until after the previous generation has died off completely.
#
# - Random mating: The model assumes that individuals mate randomly with respect to their genotype, meaning that there is no selective preference for particular genotypes.
#
# - Infinite population size: The model assumes that the population size is infinitely large, so that random sampling errors do not significantly affect the genetic composition of the population.
#
# - No mutation: The model assumes that there is no new genetic variation introduced into the population through mutation.
#
# - No selection: The model assumes that there is no selection acting on the population, meaning that all genotypes have equal fitness.
#
# - Diploid inheritance: The model assumes that each individual has two copies of each gene (i.e., is diploid) and that each parent contributes one randomly chosen copy to their offspring.
#
# These assumptions allow the Wright-Fisher model to simplify the complex dynamics of genetic change in a finite population and make predictions about how allele frequencies will change over time. However, in reality, many of these assumptions are not fully met, and more complex models may be needed to accurately predict evolutionary outcomes.

# +
# One genotype evolution fitness

import numpy as np
import matplotlib.pyplot as plt

# Define fitness values
#fitness = np.array([0.35, 0.4, 0.4])
fitness=fitness_dbem12test

# Define population size
N = 100000

# Define number of generations to simulate
T = 12

# Initialize population with equal frequencies of all genotypes
p = np.ones(len(fitness)) / len(fitness)

# Initialize list to store trajectory data
trajectory = [p.copy()]

# Simulate Wright-Fisher model
for t in range(T):
    
    # Calculate selection coefficients
    s = fitness - np.mean(fitness)
    
    # Calculate expected frequencies of each genotype in next generation
    p_next = p * np.exp(s) / np.sum(p * np.exp(s))
    
    # Sample new generation from expected frequencies
    p = np.random.multinomial(N, p_next)
    p = p / np.sum(p)
    
    # Add current population frequencies to trajectory
    trajectory.append(p.copy())

# Plot fitness landscape as a line plot
fig = plt.figure(figsize=(15, 5))

# Plot fitness landscape as a line plot
ax1 = fig.add_subplot(121)
ax1.plot(fitness)
ax1.plot(fitness,"*", markersize=10, color="b")
ax1.set_xticks(np.arange(len(fitness)))
ax1.set_xticklabels(genes2test,rotation=45)
ax1.hlines(fitness_bem1, 0, len(fitness), linestyles='dashed', colors='k')
ax1.set_xlabel('Genotype')
ax1.set_ylabel('Fitness')
ax1.set_title('Fitness Landscape')

# Plot evolutionary trajectory
ax2 = fig.add_subplot(122)
trajectory = np.array(trajectory)
for i in range(trajectory.shape[1]):
    ax2.plot(trajectory[:, i], label='$\Delta$ {}'.format(genes2test[i]))
ax2.set_xlabel('Generations')
ax2.set_ylabel('Frequency')
ax2.set_title('Evolutionary Trajectory for dbem1')
ax2.legend()
ax2.set_xlim(0, 50)
plt.tight_layout()
plt.show()

# -

plt.imshow(trajectory, cmap='Greys', aspect='auto',vmax=0.8,vmin=0)
plt.xticks(np.arange(len(genes2test)),genes2test,rotation=45);

# +
# One genotype evolution fitness

import numpy as np
import matplotlib.pyplot as plt

# Define fitness values
#fitness = np.array([0.35, 0.4, 0.4])
fitness=fitness_dbem1dbem32test

# Define population size
N = 100000

# Define number of generations to simulate
T = 12

# Initialize population with equal frequencies of all genotypes
p = np.ones(len(fitness)) / len(fitness)

# Initialize list to store trajectory data
trajectory = [p.copy()]

# Simulate Wright-Fisher model
for t in range(T):
    # Calculate selection coefficients
    s = fitness - np.mean(fitness)
    
    # Calculate expected frequencies of each genotype in next generation
    p_next = p * np.exp(s) / np.sum(p * np.exp(s))
    
    # Sample new generation from expected frequencies
    p = np.random.multinomial(N, p_next)
    p = p / np.sum(p)
    
    # Add current population frequencies to trajectory
    trajectory.append(p.copy())

# Plot fitness landscape as a line plot
fig = plt.figure(figsize=(15, 5))

# Plot fitness landscape as a line plot
ax1 = fig.add_subplot(121)
ax1.plot(fitness)
ax1.plot(fitness,"*", markersize=10, color="b")
ax1.set_xticks(np.arange(len(fitness)))
ax1.set_xticklabels(genes2test,rotation=45)
#ax1.hlines(fitness_bem1, 0, len(fitness), linestyles='dashed', colors='k')
ax1.set_xlabel('Genotype')
ax1.set_ylabel('Fitness')
ax1.set_title('Fitness Landscape')

# Plot evolutionary trajectory
ax2 = fig.add_subplot(122)
trajectory = np.array(trajectory)
for i in range(trajectory.shape[1]):
    ax2.plot(trajectory[:, i], label='$\Delta$ {}'.format(genes2test[i]))
ax2.set_xlabel('Generations')
ax2.set_ylabel('Frequency')
ax2.set_title('Evolutionary Trajectory for dbem1-dbem3')
ax2.legend()
ax2.set_xlim(0, 50)
plt.tight_layout()
plt.show()

# -

plt.imshow(trajectory, cmap='Greys', aspect='auto',vmax=0.6,vmin=0)
plt.xticks(np.arange(len(genes2test)),genes2test,rotation=45);
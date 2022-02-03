## Read per gene insertions 
import re
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from scipy import stats

def read_pergene_insertions_file(filename): 
    """Read the pergene insertions file and converts to numbers the insertion
    and reads vector

    Parameters
    ----------
    filename : str 
        Path to the pergene insertion file 

    Returns
    -------
    DataFrame
        The data as a dataframe 
    """


    data=pd.read_csv(filename,sep="\t")
    lst=[data["Reads per insertion location"]]

    total=[]
    for j in lst:
        reads_per_location=[]
        for i in data.index:
            q1=[int(s) for s in re.findall(r'\d+', j[i])] # to extract numbers from a string 
            reads_per_location.append(q1)
        total.append(reads_per_location)
    data["Reads per insertion location"]=total[0]

    lst=[data["Insertion locations"]]    

    total=[]
    for j in lst:
        insertion_locations=[]
        for i in data.index:
            q1=[int(s) for s in re.findall(r'\d+', j[i])] # to extract numbers from a string 
            insertion_locations.append(q1)
        total.append(insertion_locations)

    data["Insertion locations"]=total[0]
        
    

    return data 

# Importing data 

def getting_pergene_data(datafile):
    """[summary]

    Parameters
    ----------
    datafile : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """    

    with open(datafile) as f:
        lines = f.readlines()[1:] #skip header

    genenames_list = [None]*len(lines)
    tnpergene_list = [None]*len(lines)
    readpergene_list = [None]*len(lines) 

    line_counter = 0
    for line in lines:
        l = re.split(','' |\t', line.strip('\n'))

        genenames_list[line_counter] = l[0]
        tnpergene_list[line_counter] = int(l[1])
        readpergene_list[line_counter] = int(l[2])

        line_counter += 1

    return genenames_list,tnpergene_list,readpergene_list





def scatter_replicates(data,names,tn_lim,reads_lim,save=False):
    """Scatter plots of replicates

    Parameters
    ----------
    data : numpy.array
        An array of 4 elements where the two firsts are the transposon insertion data from 
        replicate a and b and the last two the reads data from replicate a and b.
    names : numpy.array
        an array of two strings where the first one is the name of the 
        first replicate and the second one is the name of the second replicate.
    reads_lim : int
        Upper limit for the reads to show in the figure
    tn_lim : int
        Upper limit for the transposon insertion to show in the figure
    save : bool, optional
        Whether you wish to save the figure, by default True
    """

    fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(15,8))

    res_tn = stats.linregress(np.array(data[0]), np.array(data[1]))

    ax[0].scatter(data[0],data[1],s=10,c='b',alpha=0.5,label="original data")
    ax[0].plot(np.array(data[0]), res_tn.intercept + res_tn.slope*np.array(data[0]), 'r', label='fitted line')
    ax[0].set_xlabel(names[0])
    ax[0].set_ylabel(names[1])
    ax[0].set_xlim(0,tn_lim)
    ax[0].set_ylim(0,tn_lim)
    ax[0].set_title("Transposon insertions->" + f"std: {res_tn.stderr:.6f}")
    ax[0].legend()

    res_reads = stats.linregress(np.array(data[2]), np.array(data[3]))
    ax[1].scatter(data[2],data[3],s=10,c='b',alpha=0.5,label="original data")
    ax[1].plot(np.array(data[2]), res_reads.intercept + res_tn.slope*np.array(data[2]), 'r', label='fitted line')

    ax[1].set_xlim(0,reads_lim)
    ax[1].set_ylim(0,reads_lim)
    ax[1].set_xlabel(names[0])
    ax[1].set_ylabel(names[1])
    ax[1].set_title("Reads->"+ f"std: {res_reads.stderr:.6f}")
    ax[1].legend()

    
    if save==True:
        plt.savefig("../figures/scatter_replicates"+ names[0]+"_"+ names[1]+ ".png",dpi=300,transparent=False)


def array2frame(arrays,index):
    df=pd.DataFrame(arrays)
    df=df.transpose()
    df.index=index
    df.columns=["transposons", "reads"]
    return df


def create_pergene_insertions_excel_file(datafile,export=True):
    """[summary]

    Parameters
    ----------
    datafile : [type]
        [description]
    export : bool, optional
        [description], by default True

    Returns
    -------
    [type]
        [description]
    """

    data=[]
    for i in datafile:
        data.append(read_pergene_insertions_file(i))

    for i in data:
        for j in i.index:
            if not i.loc[j,"Reads per insertion location"]==[]: 
                i.loc[j,"Reads"]=np.sum(i.loc[j,"Reads per insertion location"])-np.max(i.loc[j,"Reads per insertion location"]) #REMOVE LARGEST VALUE TO REDUCE NOISE
                i.loc[j,"Insertions"]=len(i.loc[j,"Reads per insertion location"])
            else:
                i.loc[j,"Reads"]=0
                i.loc[j,"Insertions"]=0

    ## Export the data to excel
    if export==True:
        i=0
        for j in np.arange(0,len(data)):
            data[j].to_excel("../postprocessed-data/"+datafile[i].split("/")[2]+"_pergene_insertions.xlsx")
            i=i+1
    return data
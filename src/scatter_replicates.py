import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def scatter_replicates(data,names,tn_lim,reads_lim,save=True):
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

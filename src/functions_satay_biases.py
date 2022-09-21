import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict


def transposon_bias2centromeres(wigfile_path,centromeres_genomic_location_path,save=False):
    """This function plot a cumulative plot of the number of transposon vs the distance
    to centromeres to show that transposon cluster in the proximity of them for every chromosome.

    Parameters
    ----------
    wigfile_path : str
        Relative or absolute path of where is the wig file of the processed data
    centromeres_genomic_location_path : str
        Relative or absolute path of where is the file of the centromeres genomic location
    save : bool, optional
        _description_, by default False

    Returns
    -------
    fig, all distances to centromeres from every chromosome
        
    """

    import wiggelen as wg

    ## Import wig files 
    wig_reads=[]

    for x in wg.walk(open(wigfile_path)):
        wig_reads.append(x)

    wig_reads = np.array(wig_reads)

    ## reading wig file
    pos=[]
    chrom=[]
    for i in wig_reads:
        chrom.append(i[0]) ## chromosomes
        pos.append(float(i[1])) ## positions

    chrom=np.array(chrom)
    pos=np.array(pos)
    

    ## Importing centromeres genomic locations
    centromeres=pd.read_csv(centromeres_genomic_location_path,sep=",")

    chrom_names=["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI"]

    centromeres_names_translation={"chrI":"CEN1","chrII":"CEN2","chrIII":"CEN3",
    "chrIV":"CEN4","chrV":"CEN5","chrVI":"CEN6","chrVII":"CEN7","chrVIII":"CEN8",
    "chrIX":"CEN9","chrX":"CEN10","chrXI":"CEN11","chrXII":"CEN12","chrXIII":"CEN13",
    "chrXIV":"CEN14","chrXV":"CEN15","chrXVI":"CEN16"}

    ## Compute distance to centromeres for each chromosome
    distance2cent_all=defaultdict(dict)

    for i in chrom_names:
        pos_chrom=pos[np.where(chrom==i)]
        cen_data_chrom=centromeres[centromeres.loc[:,"centromere_name"]==centromeres_names_translation[i]]
        distance2cent_all[i]=(2*np.abs(pos_chrom-cen_data_chrom.loc[:,"start"].values[0]))

    ## Plotting

    fig,axes=plt.subplots(1,1,figsize=(5,5))

    x=np.arange(0,300000,10000)
    curves=[]
    for j in chrom_names:
        distance2cent=distance2cent_all[j]
        c=[]
        
    
        for i in x:
            y=np.sum(distance2cent<i)
            c.append(y)
            
        curves.append(c)  
        axes.plot(x,c,color="gray",linewidth=0.3,alpha=0.5)
    
        

    df = pd.DataFrame(curves)
    data2fit=df.mean()


    axes.plot(x,df.mean(),color="red",label="Mean")
    axes.set_xlabel("Distance to centromere[$\pm$kb]")
    axes.set_ylabel("Number of transposons")
    axes.set_xlim(0,300000)


    ## Linear fit
    x2fit=np.arange(200000,300000,10000)
    model = np.polyfit(x2fit, data2fit[20:30], 1)
    axes.plot(x,model[0]*x+model[1],color="black",label="Linear fit",alpha=0.8)
    #plt.plot(model[0]*x+model[1],color="blue",label="Linear fit")
    plt.text(10000, 10000, 'y=%.3fx+%.1f' % (model[0], model[1]), fontsize=10)
    axes.legend()
    if save==True:
        plt.tight_layout()
        fig.savefig("../figures/fig_centromere_bias.png",transparent=True,dpi=300)
    else:
        pass
    return fig,distance2cent_all
        
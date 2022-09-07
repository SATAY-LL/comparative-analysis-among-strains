import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def linear_transformation_per_background(pergene_insertions_all_data,background,chrom_length,
windows_size=10000):
    """Executes the linear transformation procedure for a given background.
    It will divide the insertions and reads per gene over the total amount in each chromosome.It 
    will also gives the insertions and reads per gene over the total amount in each 10kb window, as default.

    Parameters
    ----------
    pergene_insertions_all_data : pandas.dataframe
        The dataframe given info about the insertions and reads per gene per background. The background 
        is given as key of the dataframe.
    background : str
        the name of the background to be used as a string, for example "wt_merged"
    chrom_length : pandas.dataframe
        A dataframe where the length of each chromosome is given and the name of the chromosomes 
        is the index of the dataframe. 
    windows_size : int, optional

    Returns
    -------
    pandas dataframe 
        A copy of the input datframe with added columns representing 
        the linear transformation of the  reads and insertions. 
    """

    chrom=["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV",
    "XV","XVI"]

    data=pergene_insertions_all_data.copy()

    data_background=data.loc[background]
    

    reads_per_chrom=data_background.groupby(by=["Chromosome"]).sum()["Reads"]#Total number of reads per chromosome
    insertions_per_chrom=data_background.groupby(by=["Chromosome"]).sum()["Insertions"]#Total number of insertions per chromosome

    data_background.index=data_background["Chromosome"]#setting the index as the chromosome number 

    
    for i in chrom:

        data_background_chrom=data_background.loc[data_background["Chromosome"]==i] # filtering the dataframe by the chromosome number
        lengths_genes=data_background_chrom.loc[:,"End location"]-data_background_chrom.loc[:,"Start location"] # computing the length of every gene

        data_background.loc[data_background["Chromosome"]==i,"Linear-norm-reads"]=data_background_chrom.loc[:,"Reads"]/reads_per_chrom[i] # dividing the reads per gene by the total number of reads in the chromosome
        data_background.loc[data_background["Chromosome"]==i,"Linear-norm-insertions"]=data_background_chrom.loc[:,"Insertions"]/insertions_per_chrom[i] # dividing the insertions per gene by the total number of insertions in the chromosome
        
        data_background.loc[data_background["Chromosome"]==i,"Linear-norm-tr-density"]=data_background.loc[data_background["Chromosome"]==i,"Linear-norm-insertions"]*np.divide(np.array(chrom_length.loc[i]),np.array(lengths_genes)) # insertions/gene length/length of chrom=insertions*length of chrom/gene length

        data_background.loc[data_background["Chromosome"]==i,"tr-density"]=data_background_chrom.loc[:,"Insertions"]/lengths_genes
        data_background.loc[data_background["Chromosome"]==i,"length gene"]=lengths_genes
        
    data_without_index=data_background.copy()
    data_without_index.drop(columns=["Chromosome"],inplace=True)
    data_without_index.reset_index(inplace=True)

 #### Normalizing taking into account a windows size around the gene instead of the whole chromosome
    for gene in data_without_index["Gene name"].unique():
    
        location_gene=[data_without_index[data_without_index["Gene name"]==gene]["End location"].values[0],
                    data_without_index[data_without_index["Gene name"]==gene]["Start location"].values[0]]

       

        windows_location=[location_gene[0]-windows_size,location_gene[1]+windows_size]


        locations_ups=np.where((windows_location[0]<data_without_index["Start location"]) & (data_without_index["Start location"]<location_gene[0]))[0] #indexes of the genes that are upstream within the window of interest
        locations_down=np.where((location_gene[1]<data_without_index["End location"]) & (data_without_index["End location"] < windows_location[1]))[0] #indexes of the genes that are downstream the windows of interest

        locations_total=np.unique(np.concatenate((locations_ups,locations_down))) # total number of gene locations inside the window.setdefault()

        total_insertions=data_without_index.loc[locations_total,"Insertions"].sum() # total number of insertionsfor genes within the window
        total_reads=data_without_index.loc[locations_total,"Reads"].sum() # total number  of reads for genes within the window

        index_gene=np.where(data_without_index["Gene name"]==gene)[0][0]
         # filling the dataframe with the new values
        data_without_index.loc[index_gene,"tr_normalized_windows"]=data_without_index.loc[index_gene,"Insertions"]/total_insertions
        data_without_index.loc[index_gene,"reads_normalized_windows"]=data_without_index.loc[index_gene,"Reads"]/total_reads
        data_without_index.loc[index_gene,"insertions_over_windows"]=total_insertions
        data_without_index.loc[index_gene,"reads_over_windows"]=total_reads

    return data_without_index,reads_per_chrom,insertions_per_chrom



def linear_transformations_plots(normalized_data,type,background,saveFigure=False):
    """_summary_

    Parameters
    ----------
    normalized_data : _type_
        _description_
    type : _type_
        _description_
    background : _type_
        _description_
    saveFigure : bool, optional
        _description_, by default False
    """        


    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(8,2))
    plt.subplots_adjust(wspace=0.8)

    data=normalized_data.loc[background]

    data=data[~data.isin([np.nan, np.inf, -np.inf]).any(1)]      

    if type=="reads":
        
        
        cols=["Reads","Linear-norm-reads","reads_normalized_windows"]
        labels=["Raw reads","Norm. over chrom.","Norm. over 10kb windows"]

        fig.suptitle("Linear transformation of the reads",y=1.1,fontsize=18)
        
        for axes in ax:
            axes.set_ylim(0,3000)

        ax[0].set_xlim(0,np.max(data[cols[0]])/40)
        ax[1].set_xlim(0,np.max(data[cols[1]])/40)
        ax[2].set_xlim(0,np.max(data[cols[2]])/40)

        ax[0].set_ylabel("Frequency",fontsize=12)
        ax[1].set_ylabel(" ")   
        ax[2].set_ylabel(" ")

        sns.histplot(data[cols[0]],ax=ax[0],color="gray",binwidth=np.max(data[cols[0]])/1000)

        ax[0].set_xlabel(labels[0],fontsize=12)
        
        sns.histplot(data[cols[1]],ax=ax[1],color="blue",binwidth=np.max(data[cols[1]])/1000)
        
        ax[1].set_xlabel(labels[1],fontsize=12)

        sns.histplot(data[cols[2]],ax=ax[2],color="red", binwidth=np.max(data[cols[2]])/1000)
            
        ax[2].set_xlabel(labels[2],fontsize=12)
   

    elif type=="insertions":
        

        cols=["Insertions","Linear-norm-insertions","tr_normalized_windows"]
        labels=["Raw insertions","Norm. over chrom.","Norm. over 10kb windows"]
        
        fig.suptitle("Linear transformation of the insertions",y=1.1,fontsize=18)

        for axes in ax:
             axes.set_ylim(0,1000)

        ax[0].set_xlim(0,np.max(data[cols[0]])/15)
        ax[1].set_xlim(0,np.max(data[cols[1]])/15)
        ax[2].set_xlim(0,np.max(data[cols[2]])/15)

        ax[0].set_ylabel("Frequency",fontsize=12)
        ax[1].set_ylabel(" ")   
        ax[2].set_ylabel(" ")
        


        sns.histplot(data[cols[0]],ax=ax[0],color="gray",binwidth=np.max(data[cols[0]])/1000)

        ax[0].set_xlabel(labels[0],fontsize=12)
        
        sns.histplot(data[cols[1]],ax=ax[1],color="blue",binwidth=np.max(data[cols[1]])/1000)
        
        ax[1].set_xlabel(labels[1],fontsize=12)

        sns.histplot(data[cols[2]],ax=ax[2],color="red", binwidth=np.max(data[cols[2]])/1000)
            
        ax[2].set_xlabel(labels[2],fontsize=12)

    elif type=="transposon density":

        cols=["tr-density","Linear-norm-tr-density","length gene"]
        labels=["Insertion density (1/bp)","Norm. over chrom.(1/bp)","Gene length (bp)"]

        
        fig.suptitle("Transposon density normalization",y=1.1,fontsize=18)

        ax[0].set_xlim(0,np.max(data[cols[0]])/15)
        ax[1].set_xlim(0,np.max(data[cols[1]])/15)
        ax[2].set_xlim(0,np.max(data[cols[2]]/4))

        for axes in ax:
            axes.set_ylim(0,400)

        ax[0].set_ylabel("Frequency",fontsize=12)
        ax[1].set_ylabel(" ")   
        ax[2].set_ylabel(" ")
        


        sns.histplot(data[cols[0]],ax=ax[0],color="gray",binwidth=np.max(data[cols[0]])/1000)

        ax[0].set_xlabel(labels[0],fontsize=12)
        
        sns.histplot(data[cols[1]],ax=ax[1],color="blue",binwidth=np.max(data[cols[1]])/1000)
        
        ax[1].set_xlabel(labels[1],fontsize=12)

        sns.histplot(data[cols[2]],ax=ax[2],color="red", binwidth=np.max(data[cols[2]])/1000)
            
        ax[2].set_xlabel(labels[2],fontsize=12)

    elif type=="plot-genome-insertions": 
        
        coordinates_chrom=[]
        chrom=["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]

        y=[]
        y_linear=[]
        y_windows=[]

        cols=["Insertions","Linear-norm-insertions","tr_normalized_windows"]
        labels=["Raw insertions","Norm. over chrom.","Norm. over 10kb"]
        
        for i in chrom:
        
            coordinates_chrom.append((np.where(data["Chromosome"]==i)[0][0]))
            
            

            y.append(data[data["Chromosome"]==i][cols[0]])
            y_linear.append(data[data["Chromosome"]==i][cols[1]])
            y_windows.append(data[data["Chromosome"]==i][cols[2]])
        

        fig.suptitle("Genome wise data normalization",y=1.1,fontsize=18)

        sns.lineplot(data=np.concatenate(y),ax=ax[0],color="black")
        sns.lineplot(data=np.concatenate(y_linear),ax=ax[1],color="blue")
        sns.lineplot(data=np.concatenate(y_windows),ax=ax[2],color="red")

        ax[0].set_ylabel(labels[0],fontsize=12)
        ax[1].set_ylabel(labels[1],fontsize=12)
        ax[2].set_ylabel(labels[2],fontsize=12)

        for axes in ax:
            axes.set_xlabel("Gene positions",fontsize=12)
            

    elif type=="plot-genome-reads": 
    
        coordinates_chrom=[]
        chrom=["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]

        y=[]
        y_linear=[]
        y_windows=[]

        cols=["Reads","Linear-norm-reads","reads_normalized_windows"]
        labels=["Raw reads","Norm. over chrom.","Norm. over 10kb"]

        for i in chrom:
        
            coordinates_chrom.append((np.where(data["Chromosome"]==i)[0][0]))

            y.append(data[data["Chromosome"]==i][cols[0]])
            y_linear.append(data[data["Chromosome"]==i][cols[1]])
            y_windows.append(data[data["Chromosome"]==i][cols[2]])

        fig.suptitle("Genome wise data normalization",y=1.1,fontsize=18)

        sns.lineplot(data=np.concatenate(y),ax=ax[0],color="black")
        sns.lineplot(data=np.concatenate(y_linear),ax=ax[1],color="blue")
        sns.lineplot(data=np.concatenate(y_windows),ax=ax[2],color="red")

        ax[0].set_ylabel(labels[0],fontsize=12)
        ax[1].set_ylabel(labels[1],fontsize=12)
        ax[2].set_ylabel(labels[2],fontsize=12)

        for axes in ax:
            axes.set_xlabel("Gene positions",fontsize=12)


    elif type=="plot-genome-density": 
    
        coordinates_chrom=[]
        chrom=["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]

        y=[]
        y_linear=[]
        y_windows=[]

        cols=["tr-density","Linear-norm-tr-density","length gene"]
        labels=["Insertion density (1/bp)","Over chrom.(1/bp)","Gene length (bp)"]

        for i in chrom:
        
            coordinates_chrom.append((np.where(data["Chromosome"]==i)[0][0]))

            y.append(data[data["Chromosome"]==i][cols[0]])
            y_linear.append(data[data["Chromosome"]==i][cols[1]])
            y_windows.append(data[data["Chromosome"]==i][cols[2]])

        fig.suptitle("Genome wise data normalization",y=1.2,fontsize=18)

        sns.lineplot(data=np.concatenate(y),ax=ax[0],color="black")
        sns.lineplot(data=np.concatenate(y_linear),ax=ax[1],color="blue")
        sns.lineplot(data=np.concatenate(y_windows),ax=ax[2],color="red")

        ax[0].set_title(labels[0],fontsize=12)
        ax[1].set_title(labels[1],fontsize=12)
        ax[2].set_title(labels[2],fontsize=12)

        for axes in ax:
            axes.set_xlabel("Gene positions",fontsize=12)
            

            
    for axes in ax:
        axes.tick_params(axis='x', labelsize=16)
        axes.tick_params(axis='y', labelsize=16)

    

    

    if saveFigure:
        fig.savefig("../figures/linear_transformation_"+background+"_"+type+".png",
        bbox_inches="tight",dpi=400)


def linear_transformations_plots_per_chrom(normalized_data,background,type,chrom="I",
saveFigure=False):
    """_summary_

    Parameters
    ----------
    normalized_data : _type_
        _description_
    background : _type_
        _description_
    type : _type_
        _description_
    chrom : str, optional
        _description_, by default "I"
    saveFigure : bool, optional
        _description_, by default False
    """


    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,2))
    
    data=normalized_data.loc[background]
    data=data[data["Chromosome"]==chrom]
    data=data[~data.isin([np.nan, np.inf, -np.inf]).any(1)] 
    
    
    if type=="plot-genome-density": 
    
        cols=["tr-density","Linear-norm-tr-density"]
        labels=["Insertion density (1/bp)","Over chrom.(1/bp)"]

        

        fig.suptitle("Normalization over chrom:"+chrom,y=1.2,fontsize=18)

        sns.lineplot(data=data.loc[:,cols[0]],ax=ax[0],color="black")
        sns.lineplot(data=data.loc[:,cols[1]],ax=ax[1],color="red")

        ax[0].set_title(labels[0],fontsize=12)
        ax[1].set_title(labels[1],fontsize=12)

    elif type=="plot-genome-insertions":
            
            cols=["Insertions","tr_normalized_windows"]
            labels=["Raw insertions","Over 10kb on chrom."]
    
            fig.suptitle("Normalization over 10kb on chrom:"+chrom,y=1.2,fontsize=18)
    
            sns.lineplot(data=data.loc[:,cols[0]],ax=ax[0],color="black")
            sns.lineplot(data=data.loc[:,cols[1]],ax=ax[1],color="red")
    
            ax[0].set_title(labels[0],fontsize=12)
            ax[1].set_title(labels[1],fontsize=12)

    elif type=="plot-genome-reads":
        
        cols=["Reads","reads_normalized_windows"]
        labels=["Raw reads","Over 10kb on chrom."]

        fig.suptitle("Normalization over 10kb on chrom:"+chrom,y=1.2,fontsize=18)

        sns.lineplot(data=data.loc[:,cols[0]],ax=ax[0],color="black")
        sns.lineplot(data=data.loc[:,cols[1]],ax=ax[1],color="red")

        ax[0].set_title(labels[0],fontsize=12)
        ax[1].set_title(labels[1],fontsize=12)
        

    for axes in ax:
        axes.tick_params(axis='x', labelsize=16)
        axes.tick_params(axis='y', labelsize=16)
        axes.set_xlabel("Gene positions",fontsize=12)
        axes.set_ylabel(" ",fontsize=12)
    
    
            

    if saveFigure:
        fig.savefig("../figures/linear_transformation_"+background+"_"+type+"_"+chrom+".png",
        bbox_inches="tight",dpi=400)
            
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from from_excel_to_list import from_excel_to_list

## Plot the reads per gene with the domain locations underneath 

def plot_reads_per_gene_with_domains(data,gene,domain_genomics_coordinates,key): 
    """Make a figure to display the domains per gene underneath the reads per insertion location
    in a gene of interest. 

    Parameters
    ----------
    data : dataframe
        Data from one key that contains the information about domains and insertion locations
    gene : str
        gene of interest
    
    domain_genomics_coordinates: dataframe 
        Data from the genomic coordinates of the protein domains

    key: str
        Background to analyse , that is part of the keys of data 


    """



    reads_vector=from_excel_to_list(data.loc[gene,"Reads per insertion location"])

    insertions_vector=from_excel_to_list(data.loc[gene,"Insertion locations"])

    start_location=data.loc[gene,"Start location"]
    end_location=data.loc[gene,"End location"]

    

    length=end_location-start_location

    length_10_first=0.1*length
    length_10_last=0.9*length

    #figure=plt.figure(figsize=(9,3))
    grid = plt.GridSpec(15, 1, wspace=0.0, hspace=0.01)

    axc = plt.subplot(grid[14,0])

    ax = plt.subplot(grid[0:14,0])
    ax.set_xlabel("Insertion location")
    ax.set_ylabel("Reads")

    ax.set_title("Reads per insertion location for " + gene )
    # plt.xlim(start_location,end_location)

    ax.vlines(ymax=np.max(reads_vector)/2,ymin=0,x=start_location+length_10_first,color="red")
    ax.vlines(ymax=np.max(reads_vector)/2,ymin=0,x=start_location+length_10_last,color="red")

    ax.annotate("10%",xy=(start_location+length_10_first,np.max(reads_vector)/2),arrowprops=dict(facecolor='black', shrink=0.05))
    ax.annotate("90%",xy=(start_location+length_10_last,np.max(reads_vector)/2),arrowprops=dict(facecolor='black', shrink=0.05))

    noncoding_color = "#002538"
    essential_color = "#10e372"
    nonessential_color = "#d9252e"
    codingdna_color = '#29a7e6'
    textcolor = "#000000"

    textsize = 10


    axc.tick_params(labelsize=textsize)
    axc.set_yticklabels([])
    axc.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

    axc.tick_params(
        axis='y',          # changes apply to the y-axis
        which='both',      # both major and minor ticks are affected
        left=False,        # ticks along the bottom edge are off
        right=False,       # ticks along the top edge are off
        labelleft=False)   # labels along the bottom edge are off

    ax.grid(linestyle='-', alpha=1.0)
    ax.tick_params(labelsize=textsize)

    ax.tick_params(axis='x', which='major', pad=30)
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    ax.xaxis.get_offset_text().set_fontsize(textsize)
    

    ax.bar(insertions_vector,reads_vector)

    ax.set_xlim(start_location,end_location)

    facecolor=[essential_color,nonessential_color,noncoding_color,codingdna_color,
    "black","blue","red","purple"]


    i=0
    var=domain_genomics_coordinates.loc[[(key,gene)],:].dropna(axis=1)
    var=np.array(var).tolist()[0]
    axc.set_xlim(start_location,end_location)

    for width in var:
        axc.axvspan(width[0],width[1],facecolor=facecolor[i],alpha=0.3)
        if type(data.loc[gene,"protein domain"])==dict:
            data.loc[gene,"protein domain"]=np.array(data.loc[gene,"protein domain"])
            if type(data.loc[gene,"protein domain"].tolist())!=str:
                axc.annotate(str(data.loc[gene,"protein domain"].tolist()[i]),xy=(width[0],0),arrowprops=dict(facecolor='black', shrink=0.05))
            else:
                axc.annotate(str(data.loc[gene,"protein domain"]),xy=(width[0],0),arrowprops=dict(facecolor='black', shrink=0.05))

        else:
            if type(data.loc[gene,"protein domain"].tolist())!=str:
                axc.annotate(str(data.loc[gene,"protein domain"].tolist()[i]),xy=(width[0],0),arrowprops=dict(facecolor='black', shrink=0.05))
            else:
                axc.annotate(str(data.loc[gene,"protein domain"]),xy=(width[0],0),arrowprops=dict(facecolor='black', shrink=0.05))

        i+=1

    
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def annotate_volcano(volcano_df,fold_change_interval,p_value_interval):
    """Annotate genes according the fold change and p values of the volcano plot.

    Parameters
    ----------
    volcano_df : pandas.DataFrame
    data that contain the information from the volcano plot
        
    fold_change_interval : numpy.array

    a vector of two elements where the first element(positive value) set the upper threshold for the positive
    values and the second element(negative value) set the lower threshold for the negative values.


    p_value_interval : numpy.array
        a vector of two elements where the first element set the upper threshold for the positive
    values and the second element set the lower threshold for the negative values.
    """

    fig = plt.figure(figsize=(19.0,9.0))#(27.0,3))
    grid = plt.GridSpec(1, 1, wspace=0.0, hspace=0.0)
    ax = plt.subplot(grid[0,0])



    colors = {False:'black', True:'red'} # based on p-value significance 
    sc = ax.scatter(x=volcano_df['fold_change'], y=volcano_df['p_value'], alpha=0.4, marker='.', c=volcano_df['significance'].apply(lambda x:colors[x]))
    ax.grid(True, which='major', axis='both', alpha=0.4)
    ax.set_xlabel('Log2 FC')
    ax.set_ylabel('-1*Log10 p-value')



    target=volcano_df[(volcano_df["fold_change"]>fold_change_interval[0]) & (volcano_df["p_value"]>p_value_interval[0])]["gene_names"]
    target_left=volcano_df[(volcano_df["fold_change"]<fold_change_interval[1]) & (volcano_df["p_value"]>p_value_interval[1])]["gene_names"]

    if len(target)!=0:
        for i in np.arange(0,len(target)):

            index=np.where(volcano_df==target.iloc[i])[0][0]

        

            trackgene_annot = ax.annotate(volcano_df.iloc[index,:]['gene_names'], (volcano_df.iloc[index,:]['fold_change'],
                                        volcano_df.iloc[index,:]['p_value']),
                                        size=10, c='green', bbox=dict(boxstyle="round", fc="w"))
            trackgene_annot.get_bbox_patch().set_alpha(0.6)

    if len(target_left)!=0:
        for i in np.arange(0,len(target_left)):

            index=np.where(volcano_df==target_left.iloc[i])[0][0]

        

            trackgene_annot = ax.annotate(volcano_df.iloc[index,:]['gene_names'], (volcano_df.iloc[index,:]['fold_change'],
                                        volcano_df.iloc[index,:]['p_value']),
                                        size=10, c='green', bbox=dict(boxstyle="round", fc="w"))
            trackgene_annot.get_bbox_patch().set_alpha(0.6)
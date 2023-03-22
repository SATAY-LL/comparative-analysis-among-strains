import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def annotate_volcano(volcano_df,figure_title,variable=None,fold_change_interval=None,p_value_interval=None,
trackgene_list=None):
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

    fig = plt.figure(figsize=(5,5))#(27.0,3))
    grid = plt.GridSpec(1, 1, wspace=0.0, hspace=0.0)
    ax = plt.subplot(grid[0,0])
    colors = {False:'#D0D3D6', True:'pink'} # based on p-value significance 
    ax.grid(True, which='major', axis='both', alpha=0.4)
    ax.set_xlabel('Interaction score',fontsize=16)
    ax.set_ylabel('-Log$_{10}$p-value',fontsize=16)
    ax.set_title(figure_title,fontsize=20)
    ax.tick_params(labelsize=14)

    if variable!=None:

        
        sc = ax.scatter(x=volcano_df['fold_change_norm'], y=volcano_df['p_value'],
        s=100,marker='.', c=volcano_df['significance'].apply(lambda x:colors[x]),label= variable)
        
        ax.legend()

    else:
        sc = ax.scatter(x=volcano_df['fold_change_norm'], y=volcano_df['p_value'],
        s=100,marker='.', c=volcano_df['significance'].apply(lambda x:colors[x]))


    if fold_change_interval!=None and p_value_interval!=None:
        target=volcano_df[(volcano_df["fold_change_norm"]>fold_change_interval[0]) & (volcano_df["p_value"]>p_value_interval[0])]["gene_names"]
        target_left=volcano_df[(volcano_df["fold_change_norm"]<fold_change_interval[1]) & (volcano_df["p_value"]>p_value_interval[1])]["gene_names"]

        if len(target)!=0:
            for i in np.arange(0,len(target)):

                index=np.where(volcano_df==target.iloc[i])[0][0]

            

                trackgene_annot = ax.annotate(volcano_df.iloc[index,:]['gene_names'], (volcano_df.iloc[index,:]['fold_change_norm'],
                                            volcano_df.iloc[index,:]['p_value']),
                                            size=12, c='green', bbox=dict(boxstyle="round", fc="w"))
                # trackgene_annot.get_bbox_patch().set_alpha(0.6)

        if len(target_left)!=0:
            for i in np.arange(0,len(target_left)):

                index=np.where(volcano_df==target_left.iloc[i])[0][0]

            

                trackgene_annot = ax.annotate(volcano_df.iloc[index,:]['gene_names'], (volcano_df.iloc[index,:]['fold_change_norm'],
                                            volcano_df.iloc[index,:]['p_value']),
                                            size=12, c='green', bbox=dict(boxstyle="round", fc="w"))
                # trackgene_annot.get_bbox_patch().set_alpha(0.6)


    if  trackgene_list != None:
        # genenames_array = volcano_df['gene_names'].to_numpy()
        for trackgene in trackgene_list:
            trackgene = trackgene.upper()
            if trackgene in volcano_df.index:
                
                trackgene_annot = ax.annotate(volcano_df.loc[trackgene,'gene_names'], (volcano_df.loc[trackgene,'fold_change_norm'], 
                volcano_df.loc[trackgene,'p_value']),size=12, c='green', bbox=dict(boxstyle="round", fc="w"))
                # trackgene_annot.get_bbox_patch().set_alpha(0.6)
            else:
                print('WARNING: %s not found' % trackgene)
    
    return fig
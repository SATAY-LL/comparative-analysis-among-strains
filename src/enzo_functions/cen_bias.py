# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 11:11:46 2023

@author: ekingma1
"""
import pandas as pd
import numpy as np
from enzo_functions.data_import_tools import wiggle

# from insertion_bias.corr_func import power_law

def profile(tn_file,step_size=10000,
            chr_file = 'info_files/genetic/chromosome_sizes.txt',
            cen_file = 'info_files/genetic/centromere_locations.txt'):
    """
    

    Parameters
    ----------
    tn_file : TYPE
        DESCRIPTION.
    step_size : TYPE, optional
        DESCRIPTION. The default is 10000.
    chr_file : TYPE, optional
        DESCRIPTION. The default is 'info_files/genetic/chromosome_sizes.txt'.
    cen_file : TYPE, optional
        DESCRIPTION. The default is 'info_files/genetic/centromere_locations.txt'.

    Returns
    -------
    d_cen : TYPE
        DESCRIPTION.
    tn_csum : TYPE
        DESCRIPTION.

    """
    
    # Importing from a file if the file is in wiggle format
    if type(tn_file) == str and tn_file.split(".")[-1] == "wig":
        df_tn = wiggle(tn_file)
    
    # Use tn_file variable directly if data is presented as DataFrame    
    elif type(tn_file) == pd.DataFrame: 
        df_tn = tn_file
        
    # Exit if data is not presented in one of the correct formats
    else:
        exit(print("The variable tn_file is not properly defined"))
       
    tn_pos = np.array( df_tn["Position"] )
    tn_chrom = np.array( df_tn["Chrom"] )
     
    # Import information about the length of the chromosomes
    chr_sizes = pd.read_csv(chr_file,sep="\t")
    longest_chr = max(chr_sizes['Length (bp)'])
    
    # Import the location of the centromeres
    cen = pd.read_csv(cen_file,sep="\t")
    
    # Pre-initialize the array for the cumulatitive sum of insertions
    n_cols = int( np.ceil(longest_chr/step_size) )
    tn_csum = np.empty( ((len(chr_sizes)-1)*2,n_cols) )
    tn_csum[:] = np.nan
    
    # Determine the distances over which the cumulatative sum is calculated
    d_cen = np.arange( 0,n_cols*step_size,step_size )
    
    # Iterate over all chromosomes
    count = 0
    for i in cen.index:    
        # Get the size of the chromosome and the start and stop of the centromere
        chrom_size = chr_sizes.loc[i,'Length (bp)']
        cen_start = cen.loc[i,'start']
        cen_stop = cen.loc[i,'stop']
        
        # Create an array of the positions to bin over
        pos_upstr = np.arange( 0,cen_start,step_size )
        pos_down = np.arange( cen_stop,chrom_size,step_size ) - cen_stop
        
        # Get the insertions that map to the current chromosome
        ins_chrom = tn_pos[ tn_chrom == cen.loc[i,'chrom'] ]
        
        # Calculate the distance of each insertion site to the centromere
        d_upstr = cen_start - ins_chrom[ ins_chrom<cen_start ]
        d_downstr = ins_chrom[ ins_chrom>cen_stop ] - cen_stop
        
        # Sort the distance arrays in ascending order
        d_upstr.sort(); d_downstr.sort()
        
        # Calculate the cumulative number of insertions with the distance to the centromere
        csum_up = [ len( d_upstr[ 0:np.argmin(d_upstr<=n) ] ) for n in pos_upstr 
                       if len(d_upstr)!=0 and np.argmin(d_upstr<=n)!=0 ]
        
        csum_down = [ len( d_downstr[ 0:np.argmin(d_downstr<=n) ] ) for n in pos_down
                        if len(d_downstr)!=0 and np.argmin(d_downstr<=n)!=0 ] 
        
        # Store the sums of different chromosomes in a single array
        tn_csum[ count,0:len(csum_up)+1 ] = [0]+csum_up
        tn_csum[ count+1,0:len(csum_down) +1] = [0]+csum_down
        
        count+=2
        
    return d_cen, tn_csum, cen


def correct(d,tn_csum,maxdist=200000,method="poly",deg=4): 
    """
    

    Parameters
    ----------
    d : TYPE
        DESCRIPTION.
    tn_csum : TYPE
        DESCRIPTION.
    maxdist : TYPE, optional
        DESCRIPTION. The default is 300000.
    method : TYPE, optional
        DESCRIPTION. The default is "plaw".
    deg : TYPE, optional
        DESCRIPTION. The default is 4.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    # Remove columns with only nan values from the tn cumulative matrix
    last_nan = max( np.where( ~np.isnan(tn_csum) )[1] )+1
    tn_csum_nonNaN = tn_csum[:,0:last_nan]
    d = d[0:last_nan]
    
    # Initiallize the array for the average cumulative tn insertions    
    tn_csum_av = np.nanmean(tn_csum_nonNaN,axis=0)
    
    # Correct the centromere bias by fitting a power-law function
    if method == "plaw": 
        d_log = np.log( d[tn_csum_av>0] )
        tn_csum_log = np.log( tn_csum_av[tn_csum_av>0] ) 
        
        p,rss = np.polyfit(d_log, tn_csum_log, 1,full = True)[0:2] 
        tn_rate = power_law( p[0],np.exp(p[1]) )
    
        return tn_rate
    
    if method == "poly":
        p = np.polyfit(d[d<maxdist],tn_csum_av[d<maxdist],deg)
        s_fun = np.poly1d(p)
        
        return s_fun
    

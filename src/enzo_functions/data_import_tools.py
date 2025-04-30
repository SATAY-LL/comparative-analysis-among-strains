# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 15:16:28 2021

@author: ekingma1
"""
import pandas as pd
import numpy as np
import operator
from  from_excel_to_list import from_excel_to_list

#text2array = data_format_conversion.text2array

def per_gene(datafile,gene=[],remove_ade=True):
    
    data = pd.read_csv(datafile, sep='\t', header=0, index_col=0) 
    
    if remove_ade:
        
        data = data.drop("ADE2")
    
    if not gene:
        
        data['Reads per insertion location'] =\
            data['Reads per insertion location'].apply(from_excel_to_list)
            
        data['Insertion locations'] =\
            data['Insertion locations'].apply(from_excel_to_list)
    
    else:
        
        data = data.loc[gene]
        data['Reads per insertion location'] =\
            data['Reads per insertion location'].apply(text2array)
            
        data['Insertion locations'] =\
            data['Insertion locations'].apply(text2array)
            
    return data


def wiggle(wig_file,remove_ade=True,as_dict=False):
    
    # Predefine the lists
    N_reads=[]
    Chrom=[]
    pos=[]
    
    # Walk through the wig file to extract the data   
    w = open(wig_file)
    for line in w:
        
        if "chrom" in line:      
            chrm = line.split("=")[-1].strip()
            
        elif "".join( line.split() ).isdigit():
            pos.append( int(line.split()[0]) )
            N_reads.append( int(line.split()[1]) )
            Chrom.append( chrm )
    
    # Store the data in a dataframe
    out_data=pd.DataFrame({'Reads':N_reads,'Position':pos,'Chrom':Chrom})
    
    # (OPTIONAL) Remove all reads mapping to the ADE2 locus. Default is True
    if remove_ade: 
        ade2_reads= np.logical_and(out_data['Chrom']=='chrXV',
                    np.logical_and(out_data['Position']>564476,
                                   out_data['Position']<566191))
        out_data = out_data.loc[map(operator.not_, ade2_reads)]
    
    # (OPTIONAL) Convert the data format into a dictionary for each chromosome
    if as_dict:      
        out_data_dict = {}
        chrom_names = out_data['Chrom'].unique()
        
        for i in chrom_names:
            chrom_indx = out_data['Chrom'] == i
            chrom_insertions = out_data.loc[chrom_indx,['Reads','Position']]
            out_data_dict[i] = chrom_insertions
            
        out_data = out_data_dict
        
    return out_data
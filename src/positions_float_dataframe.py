## Convert to a function the assembly of the positions_float_pd 
import pandas as pd
from collections import defaultdict
import numpy as np

def positions_float_dataframe(data,keys):
    """[summary]

    Parameters
    ----------
    data : dataframe
        A dataframe whose index are the standard names of the postprocessed dataset, excluding the classified as 1.

    Returns
    -------
    dataframe
        A dataframe with the locations of each gene in the genome as a list of floats, for further processing. 
    """

    from from_excel_to_list import from_excel_to_list


    positions_float=defaultdict(dict)
    for k in keys:
        data_new=data[data.loc[:,"background"]==k]

        for i in data_new.index:

            if type(data_new.loc[i,"Position"])==pd.core.series.Series:
                tmp=[]
                tmp_insertions=[]
                for j in np.arange(0,len(data_new.loc[i,"Position"])):
                    tmp.append(from_excel_to_list(data_new.loc[i,"Position"][j]))
                    positions_float[i,k]["Positions_float"]=tmp
                    tmp_insertions.append(data_new.loc[i,"Ninsertions"][j])
                    positions_float[i,k]["Ninsertions"]=tmp_insertions
            else:
                positions_float[i,k]["Positions_float"]=(from_excel_to_list(data_new.loc[i,"Position"]))
                positions_float[i,k]["Ninsertions"]=data_new.loc[i,"Ninsertions"]
            
            
            positions_float[i,k]["Feature_type"]=np.unique(data_new.loc[i,"Feature_type"])[0]
    
    positions_float_pd=pd.DataFrame.from_dict(positions_float,orient='index')

    return positions_float_pd
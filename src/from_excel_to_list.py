import numpy as np
def from_excel_to_list(x):
    """Convert a string of numbers to a list of floats

    Parameters
    ----------
    x : str
        the values from the dataframe  on insertion locations and 
        reads per insertion location. It could be other values that are str but are 
        meant to be converted to floats.

    Returns
    -------
    list
        a list of floats
    """
   
    if x=='[]': # to avoid the error of the empty string
        x=0
    elif type(x)==np.int64:
        x=x
    else:
        x=x.replace('[', '')
        x=x.replace(']', '')
        x=x.split(',')
        x=[float(x) for x in x]
    return x
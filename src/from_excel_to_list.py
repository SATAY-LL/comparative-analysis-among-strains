def from_excel_to_list(x):
    """Convert a string of numbers to a list of floats

    Parameters
    ----------
    x : str
        the values from the dataframe  on insertion locations and reads per insertion location

    Returns
    -------
    list
        a list of floats
    """
    if len(x) < 3:# to avoid the error of the empty string
        x=0
    elif x=='[]':
        x=0
    else:
        x=x.replace('[', '')
        x=x.replace(']', '')
        x=x.split(',')
        x=[float(x) for x in x]
    return x
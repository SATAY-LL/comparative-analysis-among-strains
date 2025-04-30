import pandas as pd
import numpy as np
from src.module_intergenic_model import getting_r

data_1=pd.DataFrame({'reads-per-tr':[1,2,3,4,np.nan]})
data_2=pd.DataFrame({'reads-per-tr':[1,2,3,4,np.inf]})
data_3=pd.DataFrame({'reads-per-tr':[1,2,3,4,-np.inf]})

def test_getting_r():
    """Checking the formula for the fitness and that it replaces the inf for nans and nans 
    for the maximum of the dataset"""

    rates_1=getting_r(data_1)
    rates_2=getting_r(data_2)
    rates_3=getting_r(data_3)

    data_new=[1,2,3,4,4]
    T=90

    manual_rate=np.log(np.multiply(data_new,np.sum(data_new))/(np.sum(data_new)-data_new))/T
    manual_rate=manual_rate.tolist()

    
    for i in np.arange(0,len(rates_1[0])):
        assert rates_1[0][i]==manual_rate[i],"Check the formula for the fitness in the coarse grained model"
        assert rates_1[0][i]==rates_2[0][i]==rates_3[0][i],"The rates are not the same"
    
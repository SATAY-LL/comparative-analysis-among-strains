# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3.8.10 ('satay-dev')
#     language: python
#     name: python3
# ---

import pandas as pd
import numpy as np

cdc24=pd.read_csv('../postprocessed-data/CDC24_genetic_interactions.txt',sep='\t')
rga2=pd.read_csv('../postprocessed-data/RGA2_genetic_interactions.txt',sep='\t')

# +
cdc24=cdc24.Interactor
rga2=rga2.Interactor

common_interactors=set(cdc24).intersection(set(rga2))

# +
## save common interactors  to a txt f = open('myfile.txt', 'x')

f = open('../postprocessed-data/common_interactors_Rga2_Cdc24.txt', 'w')
for item in common_interactors:
    f.write("%s\n" % item)

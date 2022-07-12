import pandas as pd
import numpy as np
import plotly.express as px
from sqlalchemy import column
import streamlit as st

st.set_page_config(page_title="Example satay data",
page_icon=":partying_face:",layout="wide")

df=pd.read_excel(
    io="postprocessed-data/wt_a_pergene_insertions.xlsx",
    engine='openpyxl',
    index_col="Unnamed: 0")
 
df.rename(columns={"Gene name":"Gene_name"},inplace=True)

#-----sidebar-----

st.sidebar.header("Please Filter here:")

genes=st.sidebar.multiselect(
    "Select the gene of interest",
    options=df["Gene_name"].unique(),
    )

reads=st.sidebar.multiselect(
    "Select a value for the total reads per gene",
    options=np.arange(0,300,20)
)

insertions=st.sidebar.multiselect(
    "Select a value for the total insertions per gene",
    options=np.arange(0,300,20)
)

chrom=st.sidebar.multiselect(
    "Select a chromosome of interest",
    options=df["Chromosome"].unique()
)

df_selection=df.query(
    "Gene_name == @genes or Insertions == @insertions or Reads == @reads or Chromosome == @chrom" 
)

st.dataframe(df_selection)

# ----MainPage----

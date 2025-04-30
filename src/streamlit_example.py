from email.policy import default
from re import template
from matplotlib.pyplot import title
import pandas as pd
import numpy as np
import plotly.express as px
from sqlalchemy import column
from from_excel_to_list import from_excel_to_list
import streamlit as st

st.set_page_config(page_title="Example satay data",
page_icon=":partying_face:",layout="wide")

st.markdown("##")
st.title(":eight_spoked_asterisk: Satay example data ")

@st.cache
def get_data_from_excel():

    df=pd.read_excel(
        io="../postprocessed-data/wt_a_pergene_insertions.xlsx",
        engine='openpyxl',
        index_col="Unnamed: 0")
 
    df.rename(columns={"Gene name":"Gene_name"},inplace=True)
    return df

df=get_data_from_excel()
#-----sidebar-----

st.sidebar.header("Please Filter here:")

genes=st.sidebar.multiselect(
    "Select the gene of interest",
    options=df["Gene_name"].unique(),
    default=df["Gene_name"][0:10]
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

st.title(":game_die: Library numbers")

st.markdown("##")

## Some viz

total_insertions=df["Insertions"].sum()

total_reads=df["Reads"].sum()


left_column,middle_column=st.columns(2)

with left_column:
    st.subheader("Total insertions:")
    st.subheader(f"{total_insertions}")
with middle_column:
    st.subheader("Total reads:")
    st.subheader(f'{total_reads}')

st.markdown("---")

st.title(" :bar_chart: Insertions and reads distributions of the selection")
st.markdown("##")

fig_histogram_insertions=px.histogram(
    df_selection["Insertions"],title="Insertions histogram",
    template="plotly_white"
)
fig_histogram_reads=px.histogram(
    df_selection["Reads"],title="Reads histograms",
    template="plotly_white"
)

fig_histogram_insertions.update_layout(
    plot_bgcolor="rgba(0,0,0,0)",
    xaxis=(dict(showgrid=False))
)

fig_histogram_reads.update_layout(
    plot_bgcolor="rgba(0,0,0,0)",
    xaxis=(dict(showgrid=False))
)

left_column,right_column=st.columns(2)
left_column.plotly_chart(fig_histogram_insertions,use_container_width=True)
right_column.plotly_chart(fig_histogram_reads,use_container_width=True)


st.markdown("---")

st.title(" :bar_chart: Insertions and reads distributions of the library")
st.markdown("##")

fig_histogram_insertions=px.histogram(
    df["Insertions"],title="Insertions histogram",
    template="plotly_white"
)
fig_histogram_reads=px.histogram(
    df["Reads"],title="Reads histograms",
    template="plotly_white"
)

fig_histogram_insertions.update_layout(
    plot_bgcolor="rgba(0,0,0,0)",
    xaxis=(dict(showgrid=False))
)

fig_histogram_reads.update_layout(
    plot_bgcolor="rgba(0,0,0,0)",
    xaxis=(dict(showgrid=False))
)

left_column,right_column=st.columns(2)
left_column.plotly_chart(fig_histogram_insertions,use_container_width=True)
right_column.plotly_chart(fig_histogram_reads,use_container_width=True)





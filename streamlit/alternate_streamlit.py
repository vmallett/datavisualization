'''
Alternate streamlit test file 

Created by: Victoria Rachleff
Created on: 4/30/24 
'''

import streamlit as st
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import altair as alt
import sys

st.set_page_config(layout="wide")

st.title("INTERACTIVE Data Visualization: Quantitative Neuropathology")

# st.write('Testing plot')

# st.image('/Users/victoriarachleff/SEA-AD_data_dashboard/multiomics_dashboard/images/neuropath_corr_scatter.svg')

# read in quant neuropath (from sea-ad website)
filepath = '/Users/victoriarachleff/SEA-AD_neuropath_working/adata_MTGneuropath.h5ad'
neuropath_adata = sc.read_h5ad(filepath)

# read in csv form
# st.write('Invisible steam of loading datasets')
filepath = '/Users/victoriarachleff/SEA-AD_neuropath_working/MTG_neuropath.csv'
neuropath_df = pd.read_csv(filepath)
# st.write('...Complete!')

source = neuropath_df

fig = alt.Chart(source).mark_circle(size=100).encode(
    alt.X(alt.repeat("column"), type='quantitative'),
    alt.Y(alt.repeat("row"), type='quantitative'),
    color='donor_pseudotime:Q',
    tooltip=['donor_ID', 'Braak', 'Thal', 'Overall AD neuropathological Change','CERAD score', 'LATE', 'donor_pseudotime'],
).properties(
    width=400,
    height=400
).repeat(
    row=['number of NeuN positive cells per area_Grey matter', 'number of 6e10 positive objects per area_Grey matter', 'number of AT8 positive cells per area_Grey matter'],
    column=['number of AT8 positive cells per area_Grey matter', 'number of 6e10 positive objects per area_Grey matter', 'number of NeuN positive cells per area_Grey matter']
).interactive()

st.altair_chart(fig)










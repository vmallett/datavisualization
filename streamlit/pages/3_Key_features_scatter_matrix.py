'''
Neuropath streamlit website

Created by: Victoria Rachleff
Created on: 4/30/24 
Updated on: 05/24/24
'''
import os
import io
import streamlit as st
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import altair as alt
from streamlit_image_zoom import image_zoom
from PIL import Image
from typing import Union, Optional, Tuple, NewType
from IPython.display import HTML
import sys
sys.path.insert(0, '..')


st.set_page_config(
    page_title="Hello",
    page_icon="👋",
    layout="wide",
)

# read in quant neuropath (from sea-ad website)
filepath = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data", "MTG_neuropath.csv")
neuropath_df = pd.read_csv(filepath)

source = neuropath_df


st.title('Scatter Matrix of  Key Quantitative Neuropathology Features')
st.write('Plot 3: Visualizing correlations between key features from the quantitative neuropathology dataset. This correlation matrix allows one to assess relationships across multiple variables readily and provides hover over funtionality for details on demand! Three key features were selected for visualization in this context: the number of neurons (NeuN positive cells), the number of amyloid beta plaques (6e10 positive objects), and the number of pTau bearing cells (AT8 positive cells).')
st.write('Interaction: Hover mouse over any data point to see details on demand including: Donor ID, Braak Stage, Thal Phase, ADNC, CERAD, LATE, and donor pseudotime.')

# interactive hover over quant neuropath figure


chart = alt.Chart(source)
fig = chart.mark_circle(size=100).encode(
    alt.X(alt.repeat("column"), type='quantitative'),
    alt.Y(alt.repeat("row"), type='quantitative'),
    color='donor_pseudotime:Q',
    tooltip=['donor_ID', 'Braak', 'Thal', 'Overall AD neuropathological Change','CERAD score', 'LATE', 'donor_pseudotime'],
).properties(
    width=350,
    height=350
).repeat(
    row=['number of NeuN positive cells per area_Grey matter', 'number of 6e10 positive objects per area_Grey matter', 'number of AT8 positive cells per area_Grey matter'],
    column=['number of AT8 positive cells per area_Grey matter', 'number of 6e10 positive objects per area_Grey matter', 'number of NeuN positive cells per area_Grey matter']
).interactive()

st.altair_chart(fig)

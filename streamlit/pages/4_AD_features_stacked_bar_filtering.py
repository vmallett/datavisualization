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


st.title('Stacked Bar with Filtering')
st.write('Plot 4: Visualizing one key metric (Number of neurons) is split across Braak Stage and ADNC')
st.write('Interaction: Use mouse to draw a selection box over a given Braak Stage(s) to see the individual breakdown of ADNC categories in non-stacked bar plot.')


st.header('Number of neurons per area')

brush = alt.selection_interval()

chart = alt.Chart(source, width = 1500)
fig = chart.mark_bar(size=100).encode(
    # x=alt.X('x:Q', title=''),
    x='Braak:N',
    y='number of NeuN positive cells per area_Grey matter:Q',
    color='Overall AD neuropathological Change:N'
# ).transform_calculate(
#     x=f'datum[{xcol_param.name}]'
).add_params(brush,
    # xcol_param
)

bars = chart.mark_bar().encode(
    x='count()',
    y='Overall AD neuropathological Change:N',
    color='Overall AD neuropathological Change:N'
).transform_filter(
    brush
)

fig & bars



st.header('Number of pTau-bearing cells per area')


brush = alt.selection_interval()

chart = alt.Chart(source, width = 1500)
fig = chart.mark_bar(size=100).encode(
    # x=alt.X('x:Q', title=''),
    x='Braak:N',
    y='number of AT8 positive cells per area_Grey matter:Q',
    color='Overall AD neuropathological Change:N'
# ).transform_calculate(
#     x=f'datum[{xcol_param.name}]'
).add_params(brush,
    # xcol_param
)

bars = chart.mark_bar().encode(
    x='count()',
    y='Overall AD neuropathological Change:N',
    color='Overall AD neuropathological Change:N'
).transform_filter(
    brush
)

fig & bars


st.header('Number of amyloid beta plaques per area')

brush = alt.selection_interval()

chart = alt.Chart(source, width = 1500)
fig = chart.mark_bar(size=100).encode(
    # x=alt.X('x:Q', title=''),
    x='Braak:N',
    y='number of 6e10 positive objects per area_Grey matter:Q',
    color='Overall AD neuropathological Change:N'
# ).transform_calculate(
#     x=f'datum[{xcol_param.name}]'
).add_params(brush,
    # xcol_param
)

bars = chart.mark_bar().encode(
    x='count()',
    y='Overall AD neuropathological Change:N',
    color='Overall AD neuropathological Change:N'
).transform_filter(
    brush
)

fig & bars
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


st.title('Histograms to visualize distribution of any dataset feature')
st.write('Plot 2: Histogram (non-binned) to visualize counts of all features in the dataset. Note that some features have only unique values (count=1 for all) and others have biological outliers (extreme values) that will limit the utility of the visualization.')
st.write('Interaction: Use the drop down menu to select any feature of interst.')

# Selectbox to choose a feature from the dataset
selected_feature = st.selectbox('Select feature for distribution', source.columns)

# Create a histogram for the selected feature
hist = alt.Chart(source).mark_bar().encode(
    alt.X(selected_feature, bin=False),
    y='count()',
    tooltip=[selected_feature, 'count()']
).properties(
    title=f'Distribution of {selected_feature}'
)

# Display the histogram using Streamlit
st.altair_chart(hist, use_container_width=True)
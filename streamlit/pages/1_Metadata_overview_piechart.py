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

st.title('Pie Charts for Quick Dataset Overview')
st.write('Plot 1: Two pie charts to visualize two categories of metadata features in the dataset. The first pie chart visualizes the breakdown of Alzheimers disease related assessments and diagnoses. The second pie chart visualizes donor-level metadata (age at death, sex) tissue quality metrics (RIN, Brain pH).')
st.write('Interaction: Use the drop down menu to select any feature of interst.')


col1, col_padding, col2 = st.columns([1, 0.05, 1])

with col1:
    # Select the column you want to visualize
    selected_column = st.selectbox('Select column for pie chart', ['Braak', 'Thal', 'Overall AD neuropathological Change', 'CERAD score', 'LATE'])

    # Calculate counts for each category
    counts = source[selected_column].value_counts().reset_index()
    counts.columns = [selected_column, 'count']

    # Create the pie chart
    fig = alt.Chart(counts).mark_arc(size=500).encode(
        theta=alt.Theta(field='count', type='quantitative'),
        color=alt.Color(field=selected_column, type='nominal'),
        tooltip=[selected_column, 'count']
    ).properties(
        title=f'Pie chart of {selected_column}'
    )

    st.altair_chart(fig)

with col2:
    # Select the column you want to visualize
    selected_column = st.selectbox('Select column for pie chart', ['Sex', 'PMI', 'Brain pH', 'RIN', 'Age at Death'])

    # Calculate counts for each category
    counts = source[selected_column].value_counts().reset_index()
    counts.columns = [selected_column, 'count']

    # Create the pie chart
    fig = alt.Chart(counts).mark_arc(size=500).encode(
        theta=alt.Theta(field='count', type='quantitative'),
        color=alt.Color(field=selected_column, type='nominal'),
        tooltip=[selected_column, 'count']
    ).properties(
        title=f'Pie chart of {selected_column}'
    )

    st.altair_chart(fig)

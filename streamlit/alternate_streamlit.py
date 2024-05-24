'''
Alternate streamlit test file 

Created by: Victoria Rachleff
Created on: 4/30/24 
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
import cv2
from typing import Union, Optional, Tuple, NewType
from IPython.display import HTML
import sys

HTML = NewType('HTML', str)

st.set_page_config(layout="wide")

st.title("Quantitative Neuropathology Interactive Data Explorer")

st.write('Data is publicly available data (sea-ad.org)')

st.write('')

col1, col2 = st.columns([1, 0.1, 1]) # column 2 creates padding between the two "real" columns

## anndata image
with col1:

    st.write('tst')


## anndata image
with col1:

    st.write('tst')

# st.image('/Users/victoriarachleff/SEA-AD_data_dashboard/multiomics_dashboard/images/neuropath_corr_scatter.svg')

# read in quant neuropath (from sea-ad website)
# st.write('Invisible steam of loading datasets')
# filepath = '/Users/victoriarachleff/SEA-AD_data_dashboard/datavisualization/data/MTG_neuropath.csv'
filepath = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "MTG_neuropath.csv")
neuropath_df = pd.read_csv(filepath)
# st.write('...Complete!')

source = neuropath_df

# interactive hover over quant neuropath figure
chart = alt.Chart(source)
fig = chart.mark_circle(size=100).encode(
    alt.X(alt.repeat("column"), type='quantitative'),
    alt.Y(alt.repeat("row"), type='quantitative'),
    color='donor_pseudotime:Q',
    tooltip=['donor_ID', 'Braak', 'Thal', 'Overall AD neuropathological Change','CERAD score', 'LATE', 'donor_pseudotime'],
).properties(
    width=300,
    height=300
).repeat(
    row=['number of NeuN positive cells per area_Grey matter', 'number of 6e10 positive objects per area_Grey matter', 'number of AT8 positive cells per area_Grey matter'],
    column=['number of AT8 positive cells per area_Grey matter', 'number of 6e10 positive objects per area_Grey matter', 'number of NeuN positive cells per area_Grey matter']
).interactive()

st.altair_chart(fig)

# New streamlit function: image zoom (working), but need higher qualtiy image
def load_image(image_path: str) -> Image.Image:
    try:
        image = Image.open(image_path)
        return image
    except Exception as e:
        st.error(f"Error loading image: {e}")
        return None

def main():
    image_path = "/Users/victoriarachleff/SEA-AD_data_dashboard/datavisualization/images/tissuewithannotation.png"  # Ensure this path is correct

    # Load the image
    image = load_image(image_path)
    if image is None:
        return

    # Display the image with zoom capabilities
    html_content = image_zoom(image, mode="scroll", size=(800, 600), keep_aspect_ratio=False, zoom_factor=4.0, increment=0.2)
    st.markdown(html_content, unsafe_allow_html=True)

if __name__ == "__main__":
    main()









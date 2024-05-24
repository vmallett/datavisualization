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


def image_zoom(image: Union[Image.Image, np.ndarray],
                mode: Optional[str] = "default",
                size: Optional[Union[int, Tuple[int, int]]] = 512,
                keep_aspect_ratio: Optional[bool] = True,
                keep_resolution: Optional[bool] = False,
                zoom_factor: Optional[Union[float, int]] = 2.0,
                increment: Optional[float] = 0.2,
            ) -> HTML:

    return component

st.set_page_config(layout="wide")

st.title("Quantitative Neuropathology")

st.write('All visualizations created from publicly available data (sea-ad.org)')

# st.image('/Users/victoriarachleff/SEA-AD_data_dashboard/multiomics_dashboard/images/neuropath_corr_scatter.svg')

# read in quant neuropath (from sea-ad website)
# filepath = '/datavisualization/data/adata_MTGneuropath.h5ad'
# neuropath_adata = sc.read_h5ad(filepath)

# read in csv form
# st.write('Invisible steam of loading datasets')
# filepath = '/Users/victoriarachleff/SEA-AD_data_dashboard/datavisualization/data/MTG_neuropath.csv'
filepath = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "MTG_neuropath.csv")
neuropath_df = pd.read_csv(filepath)
# st.write('...Complete!')

source = neuropath_df

chart = alt.Chart(source)
fig = chart.mark_circle(size=100).encode(
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


# Supported Image Formats
# PIL image
from PIL import Image

# Load image using PIL
def load_image(image_path: str) -> Image.Image:
    return Image.open(image_path)

# Function to convert PIL image to OpenCV format
def pil_to_cv2(image: Image.Image) -> np.ndarray:
    return cv2.cvtColor(np.array(image), cv2.COLOR_RGB2BGR)

# Main script
image_path = "/Users/victoriarachleff/SEA-AD_data_dashboard/datavisualization/images/tissuewithannotation.png"  # Ensure this path is correct

# Load the image
image = load_image(image_path)

# Display image with default settings
image_zoom(image)

# Display image with custom settings
image_zoom(image, mode="scroll", size=(800, 600), keep_aspect_ratio=False, zoom_factor=4.0, increment=0.2)











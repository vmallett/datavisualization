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

def load_image(image_path: str) -> Image.Image:
    try:
        image = Image.open(image_path)
        return image
    except Exception as e:
        st.error(f"Error loading image: {e}")
        return None

HTML = NewType('HTML', str)

st.set_page_config(layout="wide")

st.title("Quantitative Neuropathology Interactive Data Explorer")

st.write('All images and data are publicly available data at sea-ag.org')
url = "https://www.sea-ad.org"

st.write('Dataset description: SEA-AD is a large multicenter effort which uses postmortem human brain tissue (N=84 donors) to study Alzheimers disease. The middle temporal gyrus (MTG) was the focus of study for the current dataset. This webpage is dedicated to increasing visability of one aspect of this effort: quantitative neuropathology. To supplement the Neuropathology Image Viewer available on sea-ad.org which allows for high resolution synced zooming of both the original and masked image, this viewer is intended to introduce the dataset in more detail and then provide useful and interactive quantitative visualizations! ')

st.write('')

col1, col_padding, col2 = st.columns([1, 0.05, 1])

## anndata image
with col1:

    st.write('Example image of human middle temporal gyrus (MTG) tissue immunohistochemically stained for 6e10 (amyloid beta plaques, brown) and Iba1 (microglia, blue). Cortical layer annotations are included for layers 1, 2, 3, 4, and 5/6. Move mouse over image and scroll to zoom in!')


    # New streamlit function: image zoom (working), but need higher qualtiy image
    def main():
        image_path = "images/tissuewithannotation.png"  # Ensure this path is correct

        # Load the image
        image = load_image(image_path)
        if image is None:
            return

        # Display the image with zoom capabilities
        html_content = image_zoom(image, mode="scroll", size=(750, 550), keep_aspect_ratio=True, zoom_factor=4.0, increment=0.2)
        st.markdown(html_content, unsafe_allow_html=True)

    if __name__ == "__main__":
        main()



## anndata image
with col2:

    st.write('Example annotation mask from the same slide to visualize the features that are quantified (plaques in green and mircroglia in red. Duplex stains allow for co-localization analysis, which assesses the spatial overlap between plaques and microglia which is relevant for AD pathophysiology. Move mouse over image and scroll to zoom in!)') 

    def main():
        image_path = "images/markupmask.png"  # Ensure this path is correct

        # Load the image
        image = load_image(image_path)
        if image is None:
            return

        # Display the image with zoom capabilities
        html_content = image_zoom(image, mode="scroll", size=(750, 550), keep_aspect_ratio=True, zoom_factor=4.0, increment=0.2)
        st.markdown(html_content, unsafe_allow_html=True)

    if __name__ == "__main__":
        main()


st.title("Interactive Plots")

# read in quant neuropath (from sea-ad website)
# st.write('Invisible steam of loading datasets')
# filepath = '/Users/victoriarachleff/SEA-AD_data_dashboard/datavisualization/data/MTG_neuropath.csv'
filepath = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "MTG_neuropath.csv")
neuropath_df = pd.read_csv(filepath)
# st.write('...Complete!')

source = neuropath_df



# st.title('Parallel Coordinate Plot to Explore Donor Metadata')
# st.write('Plot 1: ')
# st.write('Interaction: ')


# chart = alt.Chart(source, width = 1500)

# fig = chart.transform_window(
#     index='count()'
# ).transform_fold(
#     ['Braak','Thal']
# ).mark_line().encode(
#     x='key:N',
#     y=alt.Y('value:N'),
#     color='donor_ID:N',
#     detail='index:N',
#     opacity=alt.value(0.5)
# )
# st.altair_chart(fig)





st.title('Scatter Matrix of  Key Quantitative Neuropathology Features')
st.write('Plot 1: Visualizing correlations between key features from the quantitative neuropathology dataset. This correlation matrix allows one to assess relationships across multiple variables readily and provides hover over funtionality for details on demand! Three key features were selected for visualization in this context: the number of neurons (NeuN positive cells), the number of amyloid beta plaques (6e10 positive objects), and the number of pTau bearing cells (AT8 positive cells).')
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


dropdown = alt.binding_select(
    options=['number of AT8 positive cells per area_Grey matter', 'number of 6e10 positive objects per area_Grey matter'],
    name='X-axis column '
)
xcol_param = alt.param(
    value='number of AT8 positive cells per area_Grey matter',
    bind=dropdown
)

dropdown = alt.binding_select(
    options=['number of AT8 positive cells per area_Grey matter', 'number of 6e10 positive objects per area_Grey matter'],
    name='X-axis column '
)
ycol_param = alt.param(
    value='number of AT8 positive cells per area_Grey matter',
    bind=dropdown
)

brush = alt.selection_interval()

chart = alt.Chart(source, width = 1500)
fig = chart.mark_circle(size=100).encode(
    x=alt.X('x:Q', title=''),
    y='number of NeuN positive cells per area_Grey matter:Q',
    color='Overall AD neuropathological Change:N'
).transform_calculate(
    x=f'datum[{xcol_param.name}]'
).add_params(
    xcol_param, ycol_param
)

bars = chart.mark_bar().encode(
    x='count()',
    y='Overall AD neuropathological Change:N',
    color='Overall AD neuropathological Change:N'
).transform_filter(
    brush
)

fig & bars

st.altair_chart(fig)


st.title('Model Predictions of Quantitative Neuropathology Markers Over Pseudotime')
st.write('Plot 2: Generated by Janna Hong. Visualizes the observed data and model predictions for various layers and features related to quant neuropath data of Alzheimers Disease patients, showing the relationship over pseudotime to infer the progression and dynamics of neuropathological markers.')
st.write('Interaction: Select a feature of interest (e.g., AT8) and cortical layer (e.g, Layer 3) from the dropdown menu to see specific model predictions.')

import streamlit.components.v1 as components

HtmlFile = open("html/MTG_INFERENCE 1.html", 'r', encoding='utf-8')
source_code = HtmlFile.read() 
# print(source_code)
components.html(source_code,  height = 600, width = 850)


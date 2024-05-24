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

url = "https://www.sea-ad.org"
st.write("check out this [link](%s)" % url)
st.markdown("check out this [link](%s)" % url)


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

st.write('Temporary note for CSE512: This is the prototype for what will be an even more descriptive and interactive website to engage with quantitative neuropathology data. I spent quite a bit of time implementing a new image zoom feature in Streamlit, but until I am able to access high resolution images (waiting on others for this), I realize the usefulness of this tool is limited. Given this is the prototype, I wanted to show that I could load and zoom the images, but plan to 1) sync the zooming and 2) use higher quality input images to see real tissue details on zoom! Additionally, I plan to add a table to neuropathology feature descriptions above the first plot (scatter matrix). I am working on this now and think that this will set the stage much better for looking at this jargony/non-standard data. Finally, a collegue Janna has been working with us on this dataset and is very interseted in data visualization and buiding a website - I included a great interactive visual that she made here (and gave credit in the plot as well). I intend to add one more visualization at the end which is a heatmap correlating all of the neuropathology features. I am still trying to think of the best way to make this plot INTERACTIVE, so until I figure that out, I am waiting to add that visual to the website. This is really an early prototype, but Im hopeful the final product will be something more engaging and useful! Final note, I am not able to remove the text below the images for now - this is something I intend to fix by the final version!')


st.write('All images and data are publicly available data at [link](%s)' % url)

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


st.title('Scatter Matrix of  Key Quantitative Neuropathology Features')
st.write('Plot 1: Visualizing correlations between key features from the quantitative neuropathology dataset. This correlation matrix allows one to assess relationships across multiple variables readily and provides hover over funtionality for details on demand! Three key features were selected for visualization in this context: the number of neurons (NeuN positive cells), the number of amyloid beta plaques (6e10 positive objects), and the number of pTau bearing cells (AT8 positive cells).')
st.write('Interaction: Hover mouse over any data point to see details on demand including: Donor ID, Braak Stage, Thal Phase, ADNC, CERAD, LATE, and donor pseudotime.')

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
    width=350,
    height=350
).repeat(
    row=['number of NeuN positive cells per area_Grey matter', 'number of 6e10 positive objects per area_Grey matter', 'number of AT8 positive cells per area_Grey matter'],
    column=['number of AT8 positive cells per area_Grey matter', 'number of 6e10 positive objects per area_Grey matter', 'number of NeuN positive cells per area_Grey matter']
).interactive()

st.altair_chart(fig)

st.title('Model Predictions of Quantitative Neuropathology Markers Over Pseudotime')
st.write('Plot 2: Generated by Janna Hong. Visualizes the observed data and model predictions for various layers and features related to quant neuropath data of Alzheimers Disease patients, showing the relationship over pseudotime to infer the progression and dynamics of neuropathological markers.')
st.write('Interaction: Select a feature of interest (e.g., AT8) and cortical layer (e.g, Layer 3) from the dropdown menu to see specific model predictions.')

import streamlit.components.v1 as components

HtmlFile = open("html/MTG_INFERENCE 1.html", 'r', encoding='utf-8')
source_code = HtmlFile.read() 
# print(source_code)
components.html(source_code,  height = 600, width = 850)


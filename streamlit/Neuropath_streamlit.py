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

import streamlit as st

st.write("# Welcome to the Neuropath Interactive Viewer! 🧠")

st.sidebar.success("Select a sub-section above.")

st.header("The purpose of this website is to increase visibility and allow for visualization of the SEA-AD quantitative neuropathology dataset with interactive plots to enable users to ask more targeted reserach questions.")

st.markdown(
"""
Dataset description: SEA-AD is a large multicenter effort which uses postmortem human brain tissue (N=84 donors) 
to study Alzheimers disease. The middle temporal gyrus (MTG) was the focus of study for the current dataset. 
This webpage is dedicated to increasing visability of one aspect of this effort: quantitative neuropathology. 
To supplement the Neuropathology Image Viewer available on sea-ad.org which allows for high resolution synced zooming 
of both the original and masked image, this viewer is intended to introduce the dataset in more detail and then provide 
useful and interactive quantitative visualizations! ')

- All images and data are available for download 
- Check the SEA-AD website here:[sea-ad.org](https://www.sea-ad.org)
"""
)


def load_image(image_path: str) -> Image.Image:
    try:
        image = Image.open(image_path)
        return image
    except Exception as e:
        st.error(f"Error loading image: {e}")
        return None

HTML = NewType('HTML', str)

col1, col_padding, col2 = st.columns([1, 0.05, 1])

## anndata image
with col1:
    st.header('Whole slide image with cortical layer annotations')


    st.write('Example image of human middle temporal gyrus (MTG) tissue immunohistochemically stained for 6e10 (amyloid beta plaques, brown) and Iba1 (microglia, blue) with cortical layer annotations. Move mouse over image and scroll to zoom in! High res coming soon...')


    # New streamlit function: image zoom (working), but need higher qualtiy image
    def main():
        image_path = "images/tissuewithannotation.png"  # Ensure this path is correct

        # Load the image
        image = load_image(image_path)
        if image is None:
            return

        # Display the image with zoom capabilities
        image_zoom(image, mode="scroll", size=(500, 300), keep_aspect_ratio=True, zoom_factor=4.0, increment=0.2)

    if __name__ == "__main__":
        main()



## anndata image
with col2:
    st.header('Extracted region of interest with feature masking')


    st.write('Example annotation mask from the same slide to visualize the features that are quantified (plaques in green and mircroglia in red. Move mouse over image and scroll to zoom in!). High res coming soon...') 

    def main():
        image_path = "images/markupmask.png"  # Ensure this path is correct

        # Load the image
        image = load_image(image_path)
        if image is None:
            return

        # Display the image with zoom capabilities
        image_zoom(image, mode="scroll", size=(500, 300), keep_aspect_ratio=True, zoom_factor=4.0, increment=0.2)

    if __name__ == "__main__":
        main()



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



st.title("Interactive Plots")


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


# dropdown = alt.binding_select(
#     options=['number of AT8 positive cells per area_Grey matter', 'number of 6e10 positive objects per area_Grey matter'],
#     name='X-axis column '
# )
# xcol_param = alt.param(
#     value='number of AT8 positive cells per area_Grey matter',
#     bind=dropdown
# )


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

# st.altair_chart(fig)


st.title('Model Predictions of Quantitative Neuropathology Markers Over Pseudotime')
st.write('Plot 5: Generated by Janna Hong using Plotly. Visualizes the observed data and model predictions for various layers and features related to quant neuropath data of Alzheimers Disease patients, showing the relationship over pseudotime to infer the progression and dynamics of neuropathological markers. On the right, the heatmap serves as a tool for visualizing multi-layered data over time, providing insights into the dynamics of the measured features.')
st.write('Interaction: Select a feature of interest (e.g., AT8) and cortical layer (e.g, Layer 3) from the dropdown menu to see specific model predictions. Option to pan, zoom, and select values of interest.')

import streamlit.components.v1 as components


col1, col_padding, col2 = st.columns([1, 0.05, 1])

with col1: 
    HtmlFile = open("html/white_MTG_INFERENCE.html", 'r', encoding='utf-8')
    source_code = HtmlFile.read() 
    # print(source_code)
    components.html(source_code,  height = 600, width = 850)

with col2: 
    HtmlFile = open("html/Pink_MTG_INFERENCE_HEATMAP.html", 'r', encoding='utf-8')
    source_code = HtmlFile.read() 
    # print(source_code)
    components.html(source_code,  height = 600, width = 850)



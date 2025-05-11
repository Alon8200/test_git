
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 16:56:59 2025

@author: YarivB02
"""

import dash
import dash_cytoscape as cyto
from dash import html, dcc
import plotly.express as px
import pandas as pd
import numpy as np
import igraph as ig
from dash.dependencies import Input, Output, State, ALL
import sys
import matplotlib.pyplot as plt
import json
from collections import Counter

print("Mp_Gwas_Zdf_v1.py is starting successfully!")

def log(message):
    print(message)
    sys.stdout.flush()

# File paths for new data
# gwas_file = "C:/Users/YarivB02/OneDrive/OneDrive - Tel-Aviv University/8_demo_data_for_tool/metabolites/gwas_blink_model_info_b_pos_v1.csv"
# metabolites_file = "C:/Users/YarivB02/OneDrive/OneDrive - Tel-Aviv University/8_demo_data_for_tool/metabolites/metabolties_info_b_pos_v1.csv"
# index_file = "C:/Users/YarivB02/OneDrive/OneDrive - Tel-Aviv University/8_demo_data_for_tool/metabolites/metabolites_ids_b_pos_v2.csv"
gwas_file = "C:/Users/YarivBNB02/OneDrive - Tel-Aviv University/8_demo_data_for_tool/metabolites/gwas_blink_model_info_b_pos_v1.csv"
metabolites_file = "C:/Users/YarivBNB02/OneDrive - Tel-Aviv University/8_demo_data_for_tool/metabolites/metabolties_info_b_pos_v1.csv"
index_file = "C:/Users/YarivBNB02/OneDrive - Tel-Aviv University/8_demo_data_for_tool/metabolites/metabolites_ids_b_pos_v2.csv"

# File path for CANOPUS data
# canopus_file = "C:/Users/YarivB02/OneDrive/OneDrive - Tel-Aviv University/8_demo_data_for_tool/sirius/canopus_formula_summary.csv"
canopus_file = "C:/Users/YarivBNB02/OneDrive - Tel-Aviv University/8_demo_data_for_tool/sirius/canopus_formula_summary.csv"

canopus_df = pd.read_csv(canopus_file)
# Filter for ClassyFire classes with probability > 0.7
filtered_classes = canopus_df[canopus_df["ClassyFire#class Probability"] > 0.7]

# Create a mapping of met_index to class
met_class_mapping = filtered_classes.groupby("met_index")["ClassyFire#class"].apply(list).to_dict()
# Alternatively, if your file uses "row_ID":
# met_class_mapping = filtered_classes.groupby("row_ID")["ClassyFire#class"].apply(list).to_dict()

# Generate unique colors for each class
unique_classes = list(set([cls for classes in met_class_mapping.values() for cls in classes]))
color_map = {cls: plt.cm.get_cmap("tab10")(i % 10) for i, cls in enumerate(unique_classes)}

# Convert colors to hex format
color_map_hex = {cls: f"rgb({int(r*255)}, {int(g*255)}, {int(b*255)})" for cls, (r, g, b, _) in color_map.items()}

# Function to assign color for a given met_index (used for the Manhattan plot)
def get_color_for_met(met_index):
    classes = met_class_mapping.get(met_index, [])
    return color_map_hex[classes[0]] if classes else "gray"

# ------------------------
# Process GWAS Data
# ------------------------
gwas_df = pd.read_csv(gwas_file)
gwas_df.rename(columns={
    "row ID": "met_index",
    "pos_all_genome": "xaxis",
    "lod": "yaxis"
}, inplace=True)

gwas_df['met_index'] = gwas_df['met_index'].astype(str)
gwas_df['-log10(P.value)'] = -np.log10(gwas_df['p_value'])
gwas_df["color"] = gwas_df["met_index"].apply(get_color_for_met)

# ------------------------
# Process Metabolite Data & Merge
# ------------------------
met_df = pd.read_csv(metabolites_file)
inf_df = pd.read_csv(index_file)

inf_df['neutral_losses_clean'] = inf_df['neutral_losses'].fillna('').str.replace(':', ',').str.replace(';', ',')
inf_df['NL_list'] = inf_df['neutral_losses_clean'].str.split(',').apply(
    lambda lst: [nl.strip() for nl in lst if nl.strip() != 'NA' and nl.strip() != '']
)


# Count frequencies of each NL across all rows
nl_counts = Counter(nl for sublist in inf_df['NL_list'] for nl in sublist)

met_nl_mapping = dict(zip(inf_df['met_index'], inf_df['NL_list']))



# # Rename columns to match desired structure
# inf_df.rename(columns={"row ID": "met_index"}, inplace=True)
# inf_df.rename(columns={"id_num": "row ID"}, inplace=True)

# Ensure both columns are of the same type before merging
inf_df['row ID'] = inf_df['row ID'].astype(str)
met_df['row ID'] = met_df['row ID'].astype(str)

# Merge to update met_index with id_num
met_df = met_df.merge(inf_df[['met_index', 'row ID']], on='row ID', how='left')
met_df['met_index'] = met_df['met_index'].astype(str)

# ------------------------
# Correlation Analysis
# ------------------------
intensity_columns = [col for col in met_df.columns if "Peak area" in col]
met_intensities = met_df[['met_index'] + intensity_columns].copy()
met_intensities.set_index('met_index', inplace=True)

# Compute correlation matrix
correlation_matrix = met_intensities.T.corr()
threshold = 0.65

# Create igraph Graph
g = ig.Graph()
g.add_vertices(len(met_df))
g.vs['met_index'] = met_df['met_index'].tolist()

edges = []
for i in range(len(correlation_matrix)):
    for j in range(i+1, len(correlation_matrix)):
        if correlation_matrix.iloc[i, j] > threshold:
            edges.append((i, j))
g.add_edges(edges)

# Convert to Dash Cytoscape format
nodes = [{"data": {"id": str(i), "met_index": str(g.vs[i]['met_index'])}} for i in range(len(g.vs))]
edges = [{"data": {"source": str(e.source), "target": str(e.target)}} for e in g.es]

# ------------------------
# Create Manhattan Plot Figure
# ------------------------
x_ticks = gwas_df.groupby('chr')['mid_position'].mean().to_dict()
unique_met = gwas_df['met_index'].nunique()

manhattan_fig = px.scatter(
    gwas_df, x="xaxis", y="-log10(P.value)",
    title=f"GWAS Manhattan Plot (mapped compounds: {unique_met})",
    labels={"xaxis": "Chromosome"},
    hover_data={
        "met_index": True,
        "snp": True,
        "chr": True,
        "chr_arm": True,
        "genotype": True,
        "genes_next_to_snp_num": True
    },
)
manhattan_fig.update_layout(
    dragmode='lasso',  # Enable lasso selection
    xaxis=dict(
        tickmode='array',
        tickvals=list(x_ticks.values()),
        ticktext=list(x_ticks.keys())
    )
)
manhattan_fig.update_traces(marker=dict(size=3, color=gwas_df["color"]))

# ------------------------
# Class Legend Data
# ------------------------
class_counts = filtered_classes["ClassyFire#class"].value_counts().to_dict()
sorted_classes = sorted(unique_classes, key=lambda cls: class_counts.get(cls, 0), reverse=True)

# ------------------------
# Build Dash Layout
# ------------------------
app = dash.Dash(__name__)
app.layout = html.Div([
    html.H1("Metabolite Correlation & Manhattan Plot"),
    
    dcc.Graph(id='manhattan-plot', figure=manhattan_fig),
    
    html.H2("Metabolite Correlation Network"),
    
    cyto.Cytoscape(
        id='cytoscape-network',
        elements=nodes + edges,
        style={'width': '100%', 'height': '500px'},
        boxSelectionEnabled=True,  # Allow multiple selection via box
        layout={'name': 'cose'},
        stylesheet=[
            {'selector': 'node', 'style': {'background-color': 'gray'}},
            {'selector': 'edge', 'style': {'width': '2px', 'line-color': 'gray'}}
        ] + [
            {'selector': f'node[met_index="{met}"]', 'style': {'background-color': get_color_for_met(met)}}
            for met in met_class_mapping.keys()
        ]
    ),
    
    html.Button("Reset Selection", id='reset-button', n_clicks=0),
    
    html.H3("Selected Metabolite: "),
    html.Div(id='output-text', style={'font-size': '20px', 'color': 'blue'}),
    
    # Class Legend as Dropdown
    html.Div([
        html.H4("Class Legend"),
        dcc.Dropdown(
            id='class-dropdown-legend',
            options=[{'label': f"{cls} ({class_counts.get(cls, 0)})", 'value': cls} for cls in sorted_classes],
            multi=True,
            placeholder="Select class(es) from legend"
        )
    ]),
    
    # Neutral Loss Legend as Dropdown (Sorted in Descending Order by Count)
    html.Div([
        html.H4("Neutral Loss Legend"),
        dcc.Dropdown(
            id='nl-dropdown-legend',
            options=[{'label': f"{nl} ({count})", 'value': nl} 
                     for nl, count in sorted(nl_counts.items(), key=lambda x: x[1], reverse=True)],
            multi=True,
            placeholder="Select neutral loss(es) from legend"
        )
    ]),

    # Tooltip Div
    html.Div(
        id='tooltip-div',
        style={
            'position': 'absolute',
            'zIndex': '1000',
            'backgroundColor': 'white',
            'border': '1px solid black',
            'padding': '5px',
            'display': 'none'
        }
    ),
    
    # Data Stores
    dcc.Store(id='selected-met-index'),
    dcc.Store(id='dummy-store')
])




# ------------------------
# Callbacks for Interactivity
# ------------------------
@app.callback(
    Output('selected-met-index', 'data'),
    [Input('manhattan-plot', 'selectedData')],
    [State('selected-met-index', 'data')],
    prevent_initial_call='initial_duplicate'
)
def store_selected_met_index_mp(selected_data, stored_selected):
    if selected_data and 'points' in selected_data:
        selected_ids = [str(point['customdata'][0]) for point in selected_data['points']]
        
        # Ensure stored_selected is a list and update with new selections
        if not stored_selected:
            stored_selected = []
        elif isinstance(stored_selected, str):
            stored_selected = [stored_selected]

        # Add only new selections, avoiding duplicates
        stored_selected = list(set(stored_selected + selected_ids))

        return stored_selected

    return dash.no_update


@app.callback(
    Output('selected-met-index', 'data', allow_duplicate=True),
    Input('cytoscape-network', 'tapNodeData'),
    State('selected-met-index', 'data'),
    prevent_initial_call='initial_duplicate'
)
def store_selected_met_index_mn(network_click, stored_selected):
    if network_click and 'met_index' in network_click:
        selected_id = network_click['met_index']

        # Ensure stored_selected is a list and update with new selections
        if not stored_selected:
            stored_selected = []
        elif isinstance(stored_selected, str):
            stored_selected = [stored_selected]

        if selected_id not in stored_selected:
            stored_selected.append(selected_id)

        return stored_selected

    return dash.no_update


@app.callback(
    [Output('selected-met-index', 'data', allow_duplicate=True),
     Output('manhattan-plot', 'figure'),
     Output('cytoscape-network', 'stylesheet')],
    [Input('reset-button', 'n_clicks'),
     Input('selected-met-index', 'data')],
    [State('manhattan-plot', 'figure')],
    prevent_initial_call='initial_duplicate'
)
def update_selection(n_clicks, selected_met_indices, current_figure):
    ctx = dash.callback_context
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None

    # Reset Selection Triggered
    if triggered_id == 'reset-button' and n_clicks > 0:
        for trace in current_figure['data']:
            trace['marker']['color'] = 'gray'
            trace['marker']['opacity'] = 1  # Ensure all points return to full opacity
        return [], current_figure, [
            {'selector': 'node', 'style': {'background-color': 'gray'}},
            {'selector': 'edge', 'style': {'width': '2px', 'line-color': 'gray'}}
        ]

    # Ensure selected_met_indices is a list (handling single & multiple selections)
    if selected_met_indices:
        if isinstance(selected_met_indices, str):  # Convert single selection to a list
            selected_met_indices = [selected_met_indices]

        # Update Manhattan plot colors
        for trace in current_figure['data']:
            trace['marker']['color'] = [
                'red' if str(met) in selected_met_indices else 'gray'
                for met in gwas_df['met_index']
            ]
            trace['marker']['opacity'] = 1  # **Fix: Prevent dimming effect**
            trace['selectedpoints'] = None  # **Fix: Disable automatic dimming effect**

        # Update Cytoscape node colors
        new_stylesheet = [{'selector': 'node', 'style': {'background-color': 'gray'}}] + [
            {'selector': f'node[met_index="{met}"]', 'style': {'background-color': 'red'}}
            for met in selected_met_indices
        ]

        return selected_met_indices, current_figure, new_stylesheet

    return dash.no_update



@app.callback(
    [Output('manhattan-plot', 'figure', allow_duplicate=True),
     Output('cytoscape-network', 'stylesheet', allow_duplicate=True)],
    [Input('class-dropdown-legend', 'value'),
     Input('nl-dropdown-legend', 'value')],
    [State('manhattan-plot', 'figure')],
    prevent_initial_call=True
)
def update_by_legend(class_selected, nl_selected, current_figure):
    log(f"Class selected: {class_selected}, NL selected: {nl_selected}")

    # Initialize default colors and stylesheet: all markers/nodes gray.
    new_colors = ['gray'] * len(gwas_df['met_index'])
    new_stylesheet = [{'selector': 'node', 'style': {'background-color': 'gray'}}]

    # If a class legend dropdown selection was made...
    if class_selected:
        selected_classes = set(class_selected)
        new_colors = [
            'red' if selected_classes.intersection(met_class_mapping.get(met, [])) else 'gray'
            for met in gwas_df['met_index']
        ]
        for met, classes in met_class_mapping.items():
            if selected_classes.intersection(classes):
                new_stylesheet.append({
                    'selector': f'node[met_index="{met}"]',
                    'style': {'background-color': 'red'}
                })

    # If a neutral loss dropdown selection was made...
    if nl_selected:
        selected_nls = set(nl_selected)
        new_colors = [
            'red' if selected_nls.intersection(met_nl_mapping.get(met, [])) else 'gray'
            for met in gwas_df['met_index']
        ]
        for met, nls in met_nl_mapping.items():
            if selected_nls.intersection(nls):
                new_stylesheet.append({
                    'selector': f'node[met_index="{met}"]',
                    'style': {'background-color': 'red'}
                })

    # Update the marker colors for every trace in the Manhattan plot.
    for trace in current_figure['data']:
        trace['marker']['color'] = new_colors

    return current_figure, new_stylesheet



if __name__ == '__main__':
    app.run_server(debug=True)
    log("Mp_Gwas_Zdf_v1.py is running successfully!")

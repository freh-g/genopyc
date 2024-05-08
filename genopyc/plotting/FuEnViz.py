import numpy as np
import dash
from dash import dcc
from dash import html
from dash.dependencies import Input,Output
import dash_cytoscape as cyto
from gprofiler import GProfiler
import igraph as ig
import networkx as nx
import os
import pandas as pd
from genopyc.mapping.geneId_mapping import *




def FuEnViz(list_of_genes):
    """
    Perform function enrichment analysis and visualize the network.

    Parameters:
    - list_of_genes (list): List of genes for analysis and visualization.

    Returns:
    - None

    Note:
    - This function performs function enrichment analysis on the provided list of genes and visualizes the network.
    - The interactome data is obtained through the HIPPIE dataset.
    - The function enrichment analysis is performed using GProfiler.
    - The network is interactive and visualized using Dash and Cytoscape.
    """
    
    def parse_interactome():
        location = os.path.dirname(os.path.realpath(__file__))
        interactome = pd.read_csv(os.path.join(location, 'data', 'hippie_interactome.sif'), header=None, sep=' ', usecols=[0, 2])
        interactome.columns = ['source', 'target']
        interactome.source = geneId_mapping(interactome.source.astype(int).tolist(), 'entrez', 'symbol')
        interactome.target = geneId_mapping(interactome.target.astype(int).tolist(), 'entrez', 'symbol')
        interactome.dropna(inplace=True)
        return interactome

    def enrichment_analysis(list_of_genes):
        gp = GProfiler(return_dataframe=True)
        enrichment = gp.profile(
            organism='hsapiens',
            query=list_of_genes,
            significance_threshold_method='bonferroni',
            no_iea=True,
            no_evidences=False
        )
        return enrichment

    def build_graph(subgraph, enrichment_dataframe):
        enrichment_dataframe.dropna(subset='p_value', inplace=True)
        ig_subgraph = ig.Graph.from_networkx(subgraph)
        pos_ = dict(zip([v['_nx_name'] for v in ig_subgraph.vs], [coord for coord in ig_subgraph.layout_auto()]))
        app = dash.Dash(__name__)
        cyto_node_data = list(zip(
            pos_.keys(),
            [coord[0] for coord in pos_.values()],
            [coord[1] for coord in pos_.values()]
        ))
        nodes = [
            {
                'data': {'id': str(_id), 'label': str(_id)},
                'position': {'x': 120 * x, 'y': 120 * y}
            }
            for _id, x, y in cyto_node_data
        ]

        edges = [
            {'data': {'source': source, 'target': target}}
            for source, target in subgraph.edges()
        ]

        elements = nodes + edges

        default_stylesheet = [
            {
                'selector': 'node',
                'style': {
                    'background-color': '#F5CEC5',
                    'border-color': 'black',
                    'border-width': '1',
                    'label': 'data(label)',
                    'width': '60',
                    'height': '60'
                }
            },
            {
                'selector': 'edge',
                'style': {
                    'line-color': 'red',
                    'width': '1'
                }
            }
        ]

        app.layout = html.Div([
            html.Header(html.H1(['Function enrichment analysis topology visualization'],
                                style={'textAlign': 'center', 'paddingBottom': '50px', 'border': '0px solid',
                                       'border-bottom': '1px solid black'})),

            html.Main([
                html.Div([
                    html.Label('P-value Slider'),
                    dcc.Slider(
                        id='pvalue_slider',
                        min=round(-np.log10(enrichment_dataframe['p_value'].max())),
                        max=round(-np.log10(enrichment_dataframe['p_value'].min())),
                        value=round(-np.log10(enrichment_dataframe['p_value'].max())),
                        marks=dict(
                            list(zip(
                                set(sorted([round(el) for el in -np.log10(enrichment_dataframe.p_value.tolist())])),
                                [{} for value in
                                 set([round(el) for el in -np.log10(enrichment_dataframe.p_value.tolist())])]))),
                        step=None),
                    html.Div(id='updatemode-output-container', style={'marginTop': 20}),
                    html.Br(style={'lineHeight': '4'}),
                    html.Label('Sources'),
                    dcc.RadioItems(
                        id='sources',
                        labelStyle={'display': 'flex'}
                    ),
                    html.Br(style={'lineHeight': '4'}),
                    html.Label('Function'),
                    dcc.Dropdown(id='function_dropdown'),
                    html.P(id='cytoscape-mouseoverNodeData-output')
                ],
                    style={'width': '20%', 'display': 'inline-block', 'float': 'left', 'paddingTop': '20px',
                           'paddingLeft': '50px'}
                ),
                html.Div([
                    cyto.Cytoscape(
                        id='cytoscape_network',
                        layout={'name': 'preset'},
                        style={'width': '100%', 'height': '800px'},
                        stylesheet=default_stylesheet,
                        elements=elements,
                        autoRefreshLayout=True
                    )
                ],
                    style={'width': '75%', 'float': 'right', 'position': 'relative', 'top': '20px'}
                )
            ])
        ])

        @app.callback(
            Output('updatemode-output-container', 'children'),
            Input('pvalue_slider', 'value')
        )
        def display_value(value):
            return '-log10(P_Value): %s' % value

        @app.callback(
            Output('sources', 'options'),
            Input('pvalue_slider', 'value')
        )
        def set_sources(selected_pvalue):
            return [{'label': i, 'value': i} for i in
                    set(enrichment_dataframe[-np.log10(enrichment_dataframe.p_value) >= selected_pvalue].source.tolist())]

        @app.callback(
            Output('function_dropdown', 'options'),
            Input('pvalue_slider', 'value'),
            Input('sources', 'value')
        )
        def set_functions(p_value, source):
            return [{'label': i, 'value': i} for i in set(
                enrichment_dataframe[
                    (-np.log10(enrichment_dataframe.p_value) >= p_value) & (
                            enrichment_dataframe.source == source)].name.tolist())]

        @app.callback(
            Output('cytoscape_network', 'stylesheet'),
            Input('sources', 'value'),
            Input('function_dropdown', 'value')
        )
        def update_network(fsource, ffunction):
            try:
                filt_enrich = enrichment_dataframe[
                    (enrichment_dataframe.name == ffunction) & (enrichment_dataframe.source == fsource)].intersections.values[0]
                new_stylesheet = [{
                    'selector': "[id='%s']" % ele,
                    'style': {
                        'background-color': 'black',
                        'line-color': 'black'
                    }
                } for ele in filt_enrich]
                return default_stylesheet + new_stylesheet
            except:
                return default_stylesheet

        app.run(debug=True)

    edge_file = parse_interactome()
    hippie_net = nx.from_pandas_edgelist(edge_file)
    subg = hippie_net.subgraph(list_of_genes)
    enrichment = enrichment_analysis(list_of_genes)
    build_graph(subg, enrichment)


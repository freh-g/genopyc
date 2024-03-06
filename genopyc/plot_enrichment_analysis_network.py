from matplotlib.legend_handler import HandlerBase
from matplotlib.text import Text
from matplotlib.legend import Legend
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import matplotlib
from gprofiler import GProfiler



def plot_enrichment_analysis_network(
    
                            list_of_genes, pvalue, colormap='cividis', edgecolor='red', mkcolor='grey', mkfsize=2500, layout='spring',
                            mklinewidths=2, alpha=1, figsize=(20, 15), savefig=False, factor=1, k=100, cbarfontsize=20, labelling=True,
                            legend=False, legend_fontsize=20, legend_titlefontsize=25, legend_location=(1.5,0.8), legend_col=1,
                            legend_labelspacing=3, legend_title='Number of Genes', legend_columnspacing=3, legend_handlelength=3,
                            size_legend_nofelements=3, cbar_orientation='horizontal', cbar_loc=(0,0), fontsize_ann=20,
                            method_of_correction='bonferroni', no_evidences=False, no_iea=True):
    """
    Perform enrichment analysis and visualize the results as a network where nodes are the functions, size of the nodes is proportional to the number of genes involved in the specific function,
    and edges exist between nodes if the functions share genes. The thickness of the edge is proportional to the number of genes involved in the specific function.

    Parameters:
    - list_of_genes (list): List of genes for enrichment analysis and network visualization.
    - pvalue (float): P-value threshold for significance.
    - colormap (str): Colormap for visualization (default: 'cividis').
    - edgecolor (str): Color of the edges (default: 'red').
    - mkcolor (str): Color of the markers (default: 'grey').
    - mkfsize (int): Marker size (default: 10000).
    - layout (str): Layout algorithm for network visualization (default: 'spring').
    - mklinewidths (int): Line width of markers (default: 2).
    - alpha (float): Transparency of markers (default: 1).
    - figsize (tuple): Figure size (default: (40,20)).
    - savefig (bool): Whether to save the figure (default: False).
    - factor (int): Scaling factor for the layout (default: 1).
    - k (int): Optimal distance between nodes for spring layout (default: 10).
    - cbarfontsize (int): Font size of color bar (default: 10).
    - labelling (bool): Whether to label the nodes (default: True).
    - legend (bool): Whether to show legend (default: False).
    - legend_fontsize (int): Font size of legend (default: 20).
    - legend_titlefontsize (int): Font size of legend title (default: 25).
    - legend_location (tuple): Location of legend (default: (0.5,0.0)).
    - legend_col (int): Number of columns in legend (default: 6).
    - legend_labelspacing (float): Spacing between legend labels (default: 1.5).
    - legend_title (str): Title of legend (default: '').
    - legend_columnspacing (float): Spacing between legend columns (default: 1.5).
    - legend_handlelength (int): Length of legend handles (default: 3).
    - size_legend_nofelements (int): Number of elements in size legend (default: 3).
    - cbar_orientation (str): Orientation of color bar (default: 'horizontal').
    - cbar_loc (tuple): Location of color bar (default: (1, 0.5)).
    - method_of_correction (str): Method of multiple testing correction (default: 'bonferroni').
    - no_evidences (bool): Exclude electronic annotations (default: False).
    - no_iea (bool): Exclude Inferred from Electronic Annotation (IEA) evidence (default: True).
    - fontsize_ann: fontsize of the annotations on the nodes of the network

    Returns:
    - Pandas DataFrame of the enriched functions
    """
    gp = GProfiler(return_dataframe=True)
    df = gp.profile(organism='hsapiens',
                    query=list_of_genes,
                    significance_threshold_method=method_of_correction,
                    no_iea=no_iea,
                    no_evidences=no_evidences)

    if df.shape[0] == 0:
        print(f"Couldn't retrieve significantly enriched functions for the query list of genes:\n\n {list_of_genes}")
        return
    def labelling_without_overlapping(x, y, list_of_annotations, ax, fontsize ,verbose=False):
        class Point:
            def __init__(self, x, y):
                self.x = x
                self.y = y

        def doOverlap(ret1, ret2):
            l1 = Point(ret1[0, 0], ret1[1, 1])
            r1 = Point(ret1[1, 0], ret1[0, 1])
            l2 = Point(ret2[0, 0], ret2[1, 1])
            r2 = Point(ret2[1, 0], ret2[0, 1])

            if l1.x >= r2.x or l2.x >= r1.x:
                return False

            if (r1.y >= l2.y or r2.y >= l1.y):
                return False

            return True

        annotations_coord = []
        for i, dot in enumerate(y):
            x_coords = x[i]
            y_coords = y[i]
            annotation = ax.annotate(str(list_of_annotations[i]),
                                     xy=(x[i], y[i]),
                                     xytext=(x_coords, y_coords),
                                     size = fontsize)

            ax.figure.canvas.draw()
            bbox = matplotlib.text.Text.get_window_extent(annotation)
            bbox_data = ax.transData.inverted().transform(bbox)
            factor = 0.2 * (bbox_data[0, 0] - bbox_data[1, 0])
            annotations_coord.append(bbox_data)

            theta = np.radians(np.linspace(1, 360 * 200, 500))
            r = np.linspace(0, max(max(zip(x, y))), len(theta))
            x_2 = r * np.cos(theta) + x_coords
            y_2 = r * np.sin(theta) + y_coords
            n = 0
            keep_cycling = True
            while keep_cycling:
                keep_cycling = False
                for ind, box in enumerate(annotations_coord[0:-1]):
                    if doOverlap(box, bbox_data):
                        annotation.set_x(x_2[n])
                        annotation.set_y(y_2[n])
                        n += 1
                        ax.figure.canvas.draw()
                        bbox = matplotlib.text.Text.get_window_extent(annotation)
                        bbox_data = ax.transData.inverted().transform(bbox)
                        annotations_coord.pop()
                        annotations_coord.append(bbox_data)
                        keep_cycling = True
                        break

    maxpv = max([-np.log10(p) for p in df.p_value.tolist()])
    for i, (s, v) in enumerate(zip(df.source.value_counts().index, df.source.value_counts())):
        data = df[(df.source == s) & (df.p_value < pvalue)].reset_index()
        if data.shape[0] == 0:
            continue
        else:
            nxen = nx.Graph()
            for i, r in data.iterrows():
                nxen.add_node(r['name'], size=r['intersection_size'], pvalue=-np.log10(r['p_value']))

            for i, r in data.iterrows():
                for index, row in data.iloc[i + 1:].reset_index().iterrows():
                    if len(set(r['intersections']).intersection(set(row['intersections']))) > 0:
                        nxen.add_edge(r['name'],
                                      row['name'],
                                      weight=len(set(r['intersections']).intersection(set(row['intersections']))))

            if layout == 'spring':
                pos_ = nx.spring_layout(nxen, k)
            elif layout == 'auto':
                ig_subgraph = ig.Graph.from_networkx(nxen)
                pos_ = dict(zip([v['_nx_name'] for v in ig_subgraph.vs], [coord for coord in ig_subgraph.layout_auto()]))

            
            
            connections = [edge[2]['weight'] for edge in nxen.edges(data=True)]
            if len(connections) != 0:
                norm_connections = [(x - min(connections)) / (max(connections) - min(connections)) for x in connections]
            else:
                norm_connections = connections

            markers = [node[1]['size'] for node in nxen.nodes(data=True)]
            if len(markers) > 1:
                norm_markers = [(x - min(markers)) / (max(markers) - min(markers)) for x in markers]
            else:
                norm_markers = markers

            norm_markers = np.clip(norm_markers, 0.3, 1)
            
             
            fig, ax = plt.subplots(figsize=figsize)
            xses, yses, lab, colors = [], [], [], []
            for node in nxen.nodes(data=True):
                xses.append(pos_[node[0]][0])
                yses.append(pos_[node[0]][1])
                lab.append(node[0])
                colors.append(node[1]['pvalue'])
            number_lab = [str(e[0]) for e in list(enumerate(lab))]
            nodez_for_legend = ax.scatter(xses, yses, s=markers)
            nodez = ax.scatter(xses, yses, s=[mkfsize * size for size in norm_markers],
                               c=colors, cmap=colormap, vmax=maxpv, alpha=alpha, edgecolors=mkcolor,
                               linewidths=mklinewidths, clip_on=False, zorder=1)

            if labelling:
                labelling_without_overlapping(xses, yses, number_lab, ax, fontsize = fontsize_ann)

            for indx, edge in enumerate(nxen.edges(data=True)):
                if edge[2]['weight'] > 0:
                    path_1 = edge[0]
                    path_2 = edge[1]
                    x0, y0 = pos_[path_1]
                    x1, y1 = pos_[path_2]
                    edgez = ax.plot(np.linspace(x0, x1), np.linspace(y0, y1),
                                    color=edgecolor,
                                    linewidth=3 * norm_connections[indx] ** 4,
                                    zorder=0)

            cbar = plt.colorbar(nodez, ax=ax, orientation=cbar_orientation, panchor=cbar_loc)
            cbar.set_label(r'$-log_{10}(p-value)$', fontsize=cbarfontsize + 4)
            cbar.ax.tick_params(labelsize=cbarfontsize)

            if legend:
                class TextHandlerB(HandlerBase):
                    def create_artists(self, legend, text, xdescent, ydescent,
                                        width, height, fontsize, trans):
                        tx = Text(width / 2., height / 2, text, fontsize=fontsize,
                                  ha="center", va="center", fontweight="bold")
                        return [tx]

                Legend.update_default_handler_map({str: TextHandlerB()})

                first_legend = fig.legend(number_lab, lab, bbox_to_anchor=legend_location, loc="lower left",fontsize=legend_fontsize,)
                plt.gca().add_artist(first_legend)

                handles, _ = nodez.legend_elements(prop="sizes", alpha=0.6, num=size_legend_nofelements)
                _, label_markers = nodez_for_legend.legend_elements(prop="sizes", alpha=0.6)

                legend = fig.legend(handles, label_markers, fontsize=legend_fontsize, loc="upper left",
                                    bbox_to_anchor=legend_location, ncol=legend_col, labelspacing=legend_labelspacing,
                                    columnspacing=legend_columnspacing, handlelength=legend_handlelength, frameon=False)

                legend.set_title(legend_title, prop={'size': legend_titlefontsize})

            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.set_title(s, fontsize=35)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            axis = plt.gca()
            axis.set_xlim([factor * x for x in axis.get_xlim()])
            axis.set_ylim([factor * y for y in axis.get_ylim()])
            plt.tight_layout()
            if savefig:
                plt.savefig(str(s) + 'enrichment_analysis.jpeg', dpi=300, bbox_inches='tight')

            plt.show()
    return df

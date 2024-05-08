import networkx as nx
from genopyc.mapping.geneId_mapping import *

def make_eqtl_network(list_of_eqtls_df, tissue=None, gene=None, variant=None, gene_symbol=False, tissue_name=False):
    """
    Generate a network based on eQTL associations.

    Parameters:
    - list_of_eqtls_df (list): List of DataFrames containing eQTL information.
    - tissue (str, optional): Filter the network by tissue. Default is None.
    - gene (str, optional): Filter the network by gene. Default is None.
    - variant (list, optional): List of variants to include in the network. Default is None.
    - gene_symbol (bool, optional): Whether to use gene symbols instead of Ensembl IDs in the network. Default is False.
    - tissue_name (bool, optional): Whether to use tissue names instead of IDs in the network. Default is False.

    Returns:
    - networkx.Graph: A networkx graph representing the eQTL network.

    Note:
    - The function generates a (networkx) network based on eQTL associations from a list of DataFrames.
    - The network can be filtered by tissue, gene, or variant.
    - If 'gene_symbol' is True, gene symbols will be used instead of Ensembl IDs in the network.
    - If 'tissue_name' is True, tissue names will be used instead of IDs in the network.
    """

    Graph = nx.Graph()
    variant_nodes = []
    gene_nodes = []
    tissue_nodes = []
    edgelist_variants_genes = []
    edgelist_genes_tissues = []
    
    if variant:
        filtered_list = []
        for dif in list_of_eqtls_df:
            if dif.iloc[0].variantid in variant:
                filtered_list.append(dif)
        list_of_eqtls_df = filtered_list
            
    
    for eq_df in list_of_eqtls_df:            
        # Decide whether to use tissue name or id in the network (sometimes only the Id is present)
        # Mapping the ids to symbols       
        if gene_symbol:
            gene_identifier_in_the_network = 'gene_symbol'
            gene_ids = eq_df.gene_id.tolist()
            gene_symbols = geneId_mapping(gene_ids, 'ensembl', 'symbol')
            eq_df[gene_identifier_in_the_network] = gene_symbols
            eq_df.dropna(subset='gene_symbol', inplace=True)
        else:
            gene_identifier_in_the_network = 'gene_id'

        if tissue_name:
            tissue_identifier_in_the_network = 'tissue_name'
            eq_df.dropna(subset='tissue_name', inplace=True)  # None Cannot be a node in the network
        else:
            tissue_identifier_in_the_network = 'tissue'
        
        # Apply the filters if needed
        if tissue:
            eq_df = eq_df[eq_df[tissue_identifier_in_the_network] == tissue]
        if gene:
            eq_df = eq_df[eq_df[gene_identifier_in_the_network] == gene]
            
        # Making the edges
        attributes = [{'beta': v, 'color': 'green' if v > 0 else 'red'} for v in eq_df.beta.tolist()]
        annotations = [{'annotation': 'is_expressed_in'} for _ in eq_df[tissue_identifier_in_the_network].tolist()]
        edgelist_variants_genes.append(list(zip(eq_df.variantid.tolist(), eq_df[gene_identifier_in_the_network].tolist(), attributes)))
        edgelist_genes_tissues.append(list(zip(eq_df[gene_identifier_in_the_network].tolist(), eq_df[tissue_identifier_in_the_network].tolist(), annotations)))
        
        # Making the nodes
        variant_nodes.append(list(set(eq_df.variantid.tolist())))
        gene_nodes.append(list(set(eq_df[gene_identifier_in_the_network].tolist())))
        tissue_nodes.append(list(set(eq_df[tissue_identifier_in_the_network].tolist())))
        
    # Adding nodes to the graph
    variant_nodes = list(sum(variant_nodes, []))
    gene_nodes = list(sum(gene_nodes, []))
    tissue_nodes = list(sum(tissue_nodes, []))
    
    Graph.add_nodes_from(variant_nodes, tipo='variant')
    Graph.add_nodes_from(gene_nodes, tipo='gene')
    Graph.add_nodes_from(tissue_nodes, tipo='tissue')
    
    # Merge the data in the edgelists
    edgelist_genes_tissues = list(sum(edgelist_genes_tissues, []))
    edgelist_variants_genes = list(sum(edgelist_variants_genes, []))
        
    Graph.add_edges_from(edgelist_variants_genes)
    Graph.add_edges_from(edgelist_genes_tissues)
    
    return Graph

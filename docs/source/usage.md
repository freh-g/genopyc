# Get variants associated with a trait

**Genopyc allow the user of querying GWAS catalog to query for SNPs associated to a complex trait. The trait has to be expressed as experimental factor ontology ID.**

```
get_associations(efotrait, verbose=False, studyid=False)
```
Retrieve SNPs associated with an EFO trait.

Parameters:
- efotrait (str): EFO id of the trait for which variants are to be retrieved (ex. EFO_0004994)
- verbose (bool, optional): If True, the function returns the steps of variant retrieval
- studyid (bool, optional): If True, the Id of the GWAS where the association was found is retrieved

Returns:
- Pandas DataFrame: Variants associated with the trait

**Once the SNPs are retrieved they can be investigated with the function:**

```
get_variants_info(idlist, chunked=False, chunksize=200)
```

Retrieve information about variants in the Homo sapiens species from the Ensembl variation database.

Parameters:
- idlist (list or str): A list of variant IDs or a single variant ID as a string.
- chunked (bool, optional): Whether to chunk the requests into smaller portions if the number of IDs exceeds 200. Default is False.
- chunksize (int, optional): The size of each chunk when chunking requests. Default is 200.

Returns:
- dict: A dictionary containing information about the variants, where keys are variant IDs and values are variant information.

Note:
This function requires internet connectivity to access the Ensembl REST API.
Ensembl accepts a maximum of 200 IDs per request. If the number of IDs exceeds 200 and chunking is not disabled,
the requests will be split into smaller chunks to retrieve variant information.

**Moreover genopyc allows to retrieve genes or variants associated to a specific trait throug DisGeNET one of the largest repositories of genes and variants associated to human diseases.**

```
get_gdas(query_list, username, password, mode)
```
Retrieve gene or variant-disease associations from the DisGeNET database.

Parameters:
- query_list (list): A list of gene symbols or variant IDs for which associations are requested.
- username (str): Username for accessing the DisGeNET API.
- password (str): Password for accessing the DisGeNET API.
- mode (str): Mode of query, either 'genes' or 'variants'.

Returns:
- pandas DataFrame: A DataFrame containing gene or variant-disease associations.

Note:
This function queries the DisGeNET <sup>1</sup> API to retrieve gene or variant-disease associations.
It requires authentication using a username and password provided by the user.
The mode parameter specifies whether the query is for genes or variants.


# Investigate genomic positions

**Genopyc includes many functions to investigate genomic positions, specifically, given a list of variants the user can retrieve the genomic positions:**

```
get_variants_position(idlist, chunked=False, chunksize=200)
```

Retrieve the positions (chromosome and position) of variant IDs from the Ensembl database for Homo sapiens.

Parameters:
- idlist (list or str): A list of variant IDs or a single variant ID as a string.
- chunked (bool, optional): Whether to chunk the requests into smaller portions. Default is False.
- chunksize (int, optional): The size of each chunk when chunking requests. Default is 200.

Returns:
- list of tuples: A list of tuples containing the variant ID, chromosome, and position.

Note:
This function requires internet connectivity to access the Ensembl REST API.

```
get_genes_position(idlist, chunked=False, chunksize=200)

```

Retrieve the coordinates of genes in the GRCh38 assembly for a list of Ensembl IDs.

Parameters:
- idlist (list or str): A list of Ensembl gene IDs or a single Ensembl gene ID as a string.
- chunked (bool, optional): Whether to chunk the requests into smaller portions if the number of IDs exceeds 200. Default is False.
- chunksize (int, optional): The size of each chunk when chunking requests. Default is 200.

Returns:
- list of tuples: A list of tuples containing the Ensembl gene ID, chromosome number, start position, and end position.

```
get_genes(ch, position, window_size=10000, pop='EUR', features=['gene'], mode='all')
```
Retrieve genes in a window centered around a genomic position and compute the distance between the position and all genes.

Parameters:
- ch (str): Chromosome identifier.
- position (int): Genomic position around which the window is centered.
- window_size (int, optional): Size of the window in base pairs. Default is 10,000.
- pop (str, optional): Population for which to retrieve gene data. Default is 'EUR' (European).
- features (list, optional): List of features to include in the retrieval. Default is ['gene']. Possible values are [band, gene, transcript, cds, exon, repeat, simple, misc, variation, somatic_variation, structural_variation, somatic_structural_variation, constrained, regulatory, motif, other_regulatory, array_probe, man].
- mode (str, optional): Retrieval mode. Options: 'all' (returns all genes and their distances), 
'complete_data' (returns complete response from the API), 'closest_forward' (returns the closest gene 
located forward of the position), 'closest_backward' (returns the closest gene located backward of 
the position), 'closest_overall' (returns the closest gene regardless of direction). Default is 'all'.

Returns:
- dict or str: Depending on the mode parameter:
- If mode is 'all', returns a dictionary where keys are gene names and values are their distances from the position.
- If mode is 'complete_data', returns the complete response from the API.
- If mode is 'closest_forward' or 'closest_backward', returns the name of the closest gene in the specified direction.
- If mode is 'closest_overall', returns the name of the closest gene regardless of direction.
- If no genes are found, returns a message indicating the absence of genes.



```
get_ov_region(snp=None, ch=None, position=None, window=500, features=['gene'], mode='region')
```
Retrieve overlap regions for a specified SNP or genomic position.

Parameters:
- snp (str, optional): The variant ID (SNP) for which overlap regions are requested. Default is None, if an snp is provided as an input also the mode have to be on 'SNP'.
- ch (str or int, optional): Chromosome name or identifier. Required if mode is 'region'.
- position (int, optional): Genomic position (start position) for which overlap regions are requested. Required if mode is 'region'.
- window (int, optional): Window size to extend the genomic region around the SNP or position. Default is 500.
- features (list, optional): List of features to include in the retrieval. Default is ['gene']. Possible values are [band, gene, transcript, cds, exon, repeat, simple, misc, variation, somatic_variation, structural_variation, somatic_structural_variation, constrained, regulatory, motif, other_regulatory, array_probe, man].
- mode (str, optional): Mode of operation, either 'SNP' or 'region'. Default is 'region' i.e. the function is working with genomic coordinates. If SNP as input we need to feed mode = 'SNP'

Returns:
- list of pandas DataFrames: A list of pandas DataFrames containing overlap regions for each specified genomic feature type.

Note:
This function queries the Ensembl REST API to retrieve overlap regions for the specified SNP or genomic position.
It returns a list of pandas DataFrames, where each DataFrame contains overlap regions for a specific genomic feature type.
    

```
get_sequence(ch, start, stop)
```

Retrieve the genomic sequence for a given genomic region.

Parameters:
- ch (str or int): Chromosome name or identifier.
- start (int): Start position of the genomic region.
- stop (int): End position of the genomic region.

Returns:
- str: The genomic sequence for the specified genomic region.

Note:
This function queries the Ensembl REST API to retrieve the genomic sequence for the specified genomic region.
It returns the genomic sequence as a string.

```
get_closest_genes(ch, position, window_size, type_of_gene=False, mode='start',position_id)
```
Find the closest upstream and downstream genes to a given genomic position.

Parameters:
- ch (str or int): Chromosome name or identifier.
- position (int): Genomic position for which closest genes are to be found.
- window_size (int): Window size to extend the genomic region around the position.
- type_of_gene (str, optional): Biotype of genes to consider. Default is False (include all genes).
- position_id (str, optional): User-defined position id.
- mode (str, optional): Mode of calculating distances, currently only distance calculated from start of upstream and downstream gene is implemented'.
    start:  In this case the distance will be calculated from the beginning of both upstream and downstream genes 
    absolute: Upstream distance is calcuated as distance (base pairs) to the end and downstream as distance to the start of the gene

Returns:
- tuple: A tuple containing the position ID, the closest upstream gene ID, and the closest downstream gene ID.

Note:
This function utilizes the `get_ov_region` function to retrieve genes within the specified window around the genomic position.
It calculates distances from the specified position to the start (or end) of genes and finds the closest upstream and downstream genes.

# Linkage Disequilibrium

**Linkage disequilibrium (LD) is an important aspect when it comes to investigate the effect of variants on downstream biological pathways. In fact many times a SNPs significantly correlated with a phenotype through GWAS is associated with many others due to this phenomena. Thus, determining which is the causative SNP is not a trivial task. With genopyc information related to LD in different ways:**

```
get_variants_in_LD(variant, r2, pop='EUR')
```
Retrieve variants in linkage disequilibrium (LD) with the specified variant.

Parameters:
- variant (str): The variant ID for which LD variants are requested.
- r2 (float): The minimum threshold for the linkage disequilibrium coefficient (r2).
- pop (str, optional): The population for which LD information is requested. Default is 'EUR' (European population).

Returns:
- list: A list of variant IDs that are in linkage disequilibrium with the specified variant, based on the given r2 threshold.

Note:
This function queries the Ensembl REST API to retrieve LD information for the specified variant.
It returns a list of variant IDs that are in LD with the specified variant, with an r2 coefficient greater than or equal to the specified threshold.
If no variants are found or an error occurs during retrieval, an empty list is returned.

```
get_ld_matrix(list_of_snps, token, pop='EUR', metric='r2')
```
Get the linkage disequilibrium (LD) matrix for a list of SNPs.

Parameters:
- list_of_snps (list): A list of SNP IDs.
- token (str): The authentication token for accessing the LDlink API.
- pop (str, optional): The population for which LD information is requested. Default is 'EUR' (European population).
- metric (str, optional): The LD metric to be used. Default is 'r2'.

Returns:
- pandas DataFrame: A DataFrame representing the LD matrix. Rows and columns are SNP IDs, and cell values represent the LD metric between corresponding SNPs.

Note:
This function queries the LDlink API to retrieve the LD matrix for the specified list of SNPs, for this reason it needs the user to be registered on LDlink in order to obtain a token.
It returns a pandas DataFrame containing the LD matrix, where rows and columns are SNP IDs, and cell values represent the LD metric.
If any LD metric values are missing (NA), they are replaced with None, and the DataFrame is converted to numeric data types.


```
pairwise_ld(ch, start, end, pop='EUR')
```
Retrieve pairwise linkage disequilibrium (LD) information for variants within a genomic region.

Parameters:
- ch (str): Chromosome name or identifier.
- start (int): Start position of the genomic region.
- end (int): End position of the genomic region.
- pop (str, optional): The population for which LD information is requested. Default is 'EUR' (European population).

Returns:
- pandas DataFrame: A DataFrame containing pairwise LD information between variants within the specified genomic region.
    Columns include 'v1' and 'v2' for variant IDs and 'r2' for the LD coefficient.

Note:
This function queries the Ensembl REST API to retrieve pairwise LD information for variants within the specified genomic region.
It returns a DataFrame containing pairwise LD information, where 'v1' and 'v2' represent variant IDs and 'r2' represents the LD coefficient.


# eQTLS

**An important aspect of the possible imapct of the variants on the phenotype, expecially in relation to complex traits, is their effect on the expression of genes. This data can be gathered through genopyc with the functions:**


```
get_eqtl_df(rsid, p_value=0.005, increase_index=False)
```
Retrieve eQTL (expression quantitative trait loci) associations as a DataFrame for a given variant.

Parameters:
- rsid (str): The variant ID (e.g., rsID) for which eQTL associations are requested.
- p_value (float, optional): Maximum p-value threshold for filtering eQTL associations. Default is 0.005.
- increase_index (bool, optional): Whether to increase the DataFrame index by 1. Default is False.

Returns:
- pandas DataFrame or None: A DataFrame containing eQTL associations for the specified variant if successful, 
    otherwise returns None.

Note:
This function queries the EBI eQTL Catalog API to retrieve eQTL associations for the specified variant.
It returns a pandas DataFrame containing eQTL association information, including variant ID, p-value, beta value, 
target gene ID, tissue, study ID, and tissue name. It optionally filters associations based on the provided p-value threshold.

**A functionality to express associations between variants and genes as network in a specific tissue is available in genopyc with the function**

```
make_eqtl_network(list_of_eqtls_df, tissue=None, gene=None, variant=None, gene_symbol=False, tissue_name=False):
```

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

**In this network there are 3 types of nodes (genes, variants, and tissue) an edge exists if the variant is affecting the transcription of that specific gene in the tissue. Moreover, the edge carries the information of the strength of the association between the variant and the transcription increase or decrease of the gene. The network can be plotted as shown in the [tutorial notebook](https://github.com/freh-g/genopyc/blob/main/tutorials/Genopyc_tutorial_notebook.ipynb).**

# From Variants to Genes

**Genopyc integrates many tools for investigating how the variants could act on genes. One already mentioned way for example is the gathering of eQTL data. To complement the analysis genopyc integrates also Ensembl variant effect predictor (VEP) <sup>2</sup>**:

```
variant_effect_predictor(idlist,input_type = 'rsid', chunked=False,chunksize=200,verbose=False,all_data=False,plot = False,save_plot='')
```
Perform Variant Effect Prediction (VEP) using the Ensembl Variant Effect Predictor for a list of variant IDs.

Parameters:
- idlist (list or str): A list of variant IDs or a single variant ID as a string.
- input_type (str, optional): Type of input IDs. Valid options are 'rsid' (default) or 'hgvs'.
- chunked (bool, optional): Whether to chunk the requests into smaller portions if the number of IDs exceeds 200. Default is False.
- chunksize (int, optional): The size of each chunk when chunking requests. Default is 200.
- verbose (bool, optional): Whether to print verbose information during processing. Default is False.
- all_data (bool, optional): Whether to return all data retrieved from VEP. Default is False.
- plot (bool, optional): Whether to generate a pie chart visualizing the consequences of the variants. Default is False.
- save_plot (str, optional): Filepath to save the generated plot. If not provided, the plot will be displayed instead.

Returns:
- list of DataFrames or dict: A list of DataFrames containing the VEP results (if `all_data` is False), or a dictionary containing all retrieved data (if `all_data` is True). if `plot` is True a pie chart with the effects of the variants is returned.

Note:
This function requires internet connectivity to access the Ensembl REST API.


**Locus to gene (L2G) pipeline <sup>3</sup> from Open Target Genetics:**

```
locus_to_gene(list_of_variants, score=0.1, output='genes')
```
Retrieve genes associated with variants from the Open Targets Genetics API.

Parameters:
- list_of_variants (list): List of variant identifiers.
- score (float, optional): Threshold score for gene association. Genes with scores higher than this threshold will be included. Default is 0.1.
- output (str, optional): Output format. 'genes' returns a list of genes, 'all' returns a DataFrame with variant, gene, and score information. Default is 'genes'.

Returns:
- list or DataFrame: Depending on the output parameter, either a list of genes or a DataFrame with variant, gene, and score information.

Note:
- The function retrieves genes associated with variants from the Open Targets Genetics API.
- The 'score' parameter filters genes based on their association score with the variants.
- 'output' can be set to 'genes' to return a list of genes or 'all' to return detailed information in a DataFrame.


# Mapping

**An utility of genopyc library is the capability of mapping gene Ids between different vocabularies such as ensembl, entrez, Uniprot and gene symbol. It can be easily done with:**

```
geneId_maping(query_list,source,target, return_not_mapped = False)
```

Maps gene identifiers from one Id to another. This function handles ensembleId, genesymbol, UniprotID and EntrezID.

Parameters:
- query_list (list): List of gene identifiers to be mapped.
- source (str): Source Id of the gene identifiers. Possible values are: ensembl, symbol, entrez, uniprot
- target (str): Target Id to which the gene identifiers will be mapped. Possible values are: ensembl, symbol, entrez, uniprot
- return_not_mapped (bool, default: False): print the genes that weren't mapped

Returns:
- list: List of gene identifiers mapped to the target namespace.

Raises:
- AssertionError: If input query_list contains non-integer values when source or target is 'entrez'.

**To do this genopyc relies on ncbi id mapping and and Uniprot gene info files. The datasets can be updated to the most recent version by running:**

```
update_mapping_datasets()
```

**Moreover, genopyc handles variant Id mapping through the function:**

```
variantId_mapping(list_of_variants, source='variantid', target='rsid')
```

Convert variants from one identifier type to another using the Open Targets Genetics API.

Parameters:
- list_of_variants (list): List of variant identifiers to be converted.
- source (str, optional): Identifier type of the variants in list_of_variants. Default is 'variantid'. Possible values are: 'rsid' and 'variantid'
- target (str, optional): Identifier type to which the variants will be converted. Default is 'rsid'. Possible values are: 'rsid' and 'variantid'

Returns:
- list: List of converted variant identifiers.

Note:
- Currently supports conversion between 'variantid' and 'rsid'.
- If a variant cannot be converted, it will be skipped, and a message will be printed.



# Functional Enrichment

**Once variants have been connected to genes, a functional enrichment can be performed in order to understand which pathways could be perturbated from the variants. We exploit GProfiler <sup>4</sup> to carry out the enrichment analysis through the function:**

```
plot_enrichment_analysis_network
```

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
- legend_location (tuple): position of legend (default: (0.5,0.0)).
- legend_col (int): Number of columns in legend (default: 6).
- legend_labelspacing (float): Spacing between legend labels (default: 1.5).
- legend_title (str): Title of legend (default: '').
- legend_columnspacing (float): Spacing between legend columns (default: 1.5).
- legend_handlelength (int): Length of legend handles (default: 3).
- size_legend_nofelements (int): Number of elements in size legend (default: 3).
- cbar_orientation (str): Orientation of color bar (default: 'horizontal').
- cbar_loc (tuple): position of color bar (default: (1, 0.5)).
- method_of_correction (str): Method of multiple testing correction (default: 'bonferroni').
- no_evidences (bool): Exclude electronic annotations (default: False).
- no_iea (bool): Exclude Inferred from Electronic Annotation (IEA) evidence (default: True).
- **kwargs: Additional keyword arguments.

Returns:
- Pandas DataFrame of the enriched functions

**The returned plot is a network where the nodes are the function.**
\
**A new functionality implemented in genopyc is the interactive visualization of the function enrichment through the function**
```
FuEnViz
```


**In this network the nodes are proteins and an edge exists between the proteins if there is an interaction. A dropdown menu listing all the enriched functions is utlized in order to highlight the proteins associated to that specific pathway.**\
**A video showing an example is available in the [GitHub repository](https://github.com/freh-g/genopyc/blob/main/function_enrichment_visualization_tool_tutorial_supplementary.mp4)**


**REFERENCES**

[1] Piñero, J., Ramírez-Anguita, J. M., Saüch-Pitarch, J., Ronzano, F., Centeno, E., Sanz, F., & Furlong, L. I. (2019). The DisGeNET knowledge platform for disease genomics: 2019 update. Nucleic Acids Research, 48(D1), D845–D855. https://doi.org/10.1093/nar/gkz1021

[2] McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1). https://doi.org/10.1186/s13059-016-0974-4

[3] Mountjoy, E., Schmidt, E. M., Carmona, M., Schwartzentruber, J., Peat, G., Miranda, A., Fumis, L., Hayhurst, J., Buniello, A., Karim, M. A., Wright, D., Hercules, A., Papa, E., Fauman, E. B., Barrett, J. C., Todd, J. A., Ochoa, D., Dunham, I., & Ghoussaini, M. (2021). An open approach to systematically prioritize causal variants and genes at all published human GWAS trait-associated loci. Nature Genetics, 53(11), 1527–1533. https://doi.org/10.1038/s41588-021-00945-5

[4] Raudvere, U., Kolberg, L., Kuzmin, I., Arak, T., Adler, P., Peterson, H., & Vilo, J. (2019). g:Profiler: a web server for functional enrichment analysis and conversions of gene lists (2019 update). Nucleic Acids Research, 47(W1), W191–W198. https://doi.org/10.1093/nar/gkz369
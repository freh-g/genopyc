ncbidb = pd.read_csv(os.path.join(location, 'data', 'Homo_sapiens.gene_info.gz'),sep='\t')
tab_uni =  pd.read_csv(os.path.join(location, 'data', 'HUMAN_9606_idmapping.dat.gz'),sep='\t',names=['uniprot','mapper','id'])

def gene_mapping_many(query_list,source,target):
    
    """
    Maps gene identifiers from one namespace to another.

    Args:
        query_list (list): List of gene identifiers to be mapped.
        source (str): Source namespace of the gene identifiers.
        target (str): Target namespace to which the gene identifiers will be mapped.

    Returns:
        list: List of gene identifiers mapped to the target namespace.

    Raises:
        AssertionError: If input query_list contains non-integer values when source or target is 'entrez'.

    Note:
        The function relies on external data sources 'ncbidb' and 'tab_uni'.
    """
    
    

    if source=='ensembl' and target=='symbol':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((x.split('|')[0], y) for x, y in list(zip(ense,ncbidb['Symbol'])) if 'ENS' in x)
        return list(map(dictio.get,query_list))

    elif source=='ensembl' and target=='entrez':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((x.split('|')[0], y) for x, y in list(zip(ense,ncbidb['GeneID'])) if 'ENS' in x)
        return list(map(dictio.get,query_list))

    elif source=='ensembl' and target=='uniprot':
        dictio=dict(zip(list(reversed(tab_uni[tab_uni.mapper=='Ensembl'].id.tolist())),
                list(reversed(tab_uni[tab_uni.mapper=='Ensembl'].uniprot.tolist()))))
        return list(map(dictio.get,query_list))

    elif source=='entrez' and target=='ensembl':
        assert all(isinstance(x, int) for x in query_list), "All entrez_ids in query list should be integers"
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((y, x.split('|')[0]) for x, y in list(zip(ense,ncbidb['GeneID'])) if 'ENS' in x)
        return list(map(dictio.get,query_list))

    elif source == 'entrez' and target == 'symbol':
        assert all(isinstance(x, int) for x in query_list), "All entrez_ids in query list should be integers"
        dictio=dict((x, y) for x, y in list(zip(ncbidb['GeneID'],ncbidb['Symbol'])))
        return list(map(dictio.get,query_list))

    elif source=='entrez' and target=='uniprot':
        assert all(isinstance(x, int) for x in query_list), "All entrez_ids in query list should be integers"
        dictio=dict(zip([int(g) for g in list(reversed(tab_uni[tab_uni.mapper=='GeneID'].id.tolist()))],
                list(reversed(tab_uni[tab_uni.mapper=='GeneID'].uniprot.tolist()))))
        return list(map(dictio.get,query_list))

    elif source=='symbol' and target=='ensembl':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((y, x.split('|')[0]) for x, y in list(zip(ense,ncbidb['Symbol'])) if 'ENS' in x)
        return list(map(dictio.get,query_list))

    elif source == 'symbol' and target == 'entrez':
        dictio=dict((y, x) for x, y in list(zip(ncbidb['GeneID'],ncbidb['Symbol'])))
        return list(map(dictio.get,query_list))

    elif source=='symbol' and target=='uniprot':
        dictio=dict(zip(list(reversed(tab_uni[tab_uni.mapper=='Gene_Name'].id.tolist())),
                list(reversed(tab_uni[tab_uni.mapper=='Gene_Name'].uniprot.tolist()))))
        return list(map(dictio.get,query_list))

    elif source=='uniprot' and target=='ensembl':
        dictio=dict(zip(list(reversed(tab_uni[tab_uni.mapper=='Ensembl'].uniprot.tolist())),
                list(reversed(tab_uni[tab_uni.mapper=='Ensembl'].id.tolist()))))
        return list(map(dictio.get,query_list))

    elif source=='uniprot' and target=='entrez':
        dictio=dict(zip(list(reversed(tab_uni[tab_uni.mapper=='GeneID'].uniprot.tolist())),
                [int(g) for g in list(reversed(tab_uni[tab_uni.mapper=='GeneID'].id.tolist()))]))
        return list(map(dictio.get,query_list))

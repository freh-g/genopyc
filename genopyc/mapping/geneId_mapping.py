import pandas as pd
import os
import re


location = os.path.dirname(os.path.realpath(__file__))
location = location.rstrip('/mapping')
ncbidb = pd.read_csv(os.path.join(location, 'data', 'Homo_sapiens.gene_info.gz'),sep='\t')
tab_uni =  pd.read_csv(os.path.join(location, 'data', 'HUMAN_9606_idmapping.dat.gz'),sep='\t',names=['uniprot','mapper','id'])

def geneId_mapping(query_list,source,target,return_not_mapped = False):
    
    """
    Maps gene identifiers from one Id to another. This function handles ensembleId, genesymbol, UniprotID and EntrezID.

    Args:
        query_list (list): List of gene identifiers to be mapped.
        source (str): Source Id of the gene identifiers.
        target (str): Target Id to which the gene identifiers will be mapped.
        return_not_mapped (bool, default: False): print the genes that weren't mapped

    Returns:
        list: List of gene identifiers mapped to the target namespace.

    Raises:
        AssertionError: If input query_list contains non-integer values when source or target is 'entrez'.

    """
    
    if not isinstance(query_list,list):
        query_list = [query_list]
    if source=='ensembl' and target=='symbol':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((x.split('|')[0], y) for x, y in list(zip(ense,ncbidb['Symbol'])) if 'ENS' in x)
        if return_not_mapped:
            for (original,mapped) in zip(query_list,list(map(dictio.get,query_list))):
                map_str = str(mapped)
                if map_str == "None":
                    print(f"Couldn't Map {original}")
        
        return list(map(dictio.get,query_list))

    elif source=='ensembl' and target=='entrez':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((x.split('|')[0], y) for x, y in list(zip(ense,ncbidb['GeneID'])) if 'ENS' in x)
        if return_not_mapped:
            for (original,mapped) in zip(query_list,list(map(dictio.get,query_list))):
                map_str = str(mapped)
                if map_str == "None":
                    print(f"Couldn't Map {original}")
        
        return list(map(dictio.get,query_list))

    elif source=='ensembl' and target=='uniprot':
        dictio=dict(zip(list(reversed(tab_uni[tab_uni.mapper=='Ensembl'].id.tolist())),
                list(reversed(tab_uni[tab_uni.mapper=='Ensembl'].uniprot.tolist()))))
        if return_not_mapped:
            for (original,mapped) in zip(query_list,list(map(dictio.get,query_list))):
                map_str = str(mapped)
                if map_str == "None":
                    print(f"Couldn't Map {original}")
        return list(map(dictio.get,query_list))

    elif source=='entrez' and target=='ensembl':
        assert all(isinstance(x, int) for x in query_list), "All entrez_ids in query list should be integers"
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((y, x.split('|')[0]) for x, y in list(zip(ense,ncbidb['GeneID'])) if 'ENS' in x)
        if return_not_mapped:
            for (original,mapped) in zip(query_list,list(map(dictio.get,query_list))):
                map_str = str(mapped)
                if map_str == "None":
                    print(f"Couldn't Map {original}")
        return list(map(dictio.get,query_list))

    elif source == 'entrez' and target == 'symbol':
        assert all(isinstance(x, int) for x in query_list), "All entrez_ids in query list should be integers"
        dictio=dict((x, y) for x, y in list(zip(ncbidb['GeneID'],ncbidb['Symbol'])))
        if return_not_mapped:
            for (original,mapped) in zip(query_list,list(map(dictio.get,query_list))):
                map_str = str(mapped)
                if map_str == "None":
                    print(f"Couldn't Map {original}")
        return list(map(dictio.get,query_list))

    elif source=='entrez' and target=='uniprot':
        assert all(isinstance(x, int) for x in query_list), "All entrez_ids in query list should be integers"
        dictio=dict(zip([int(g) for g in list(reversed(tab_uni[tab_uni.mapper=='GeneID'].id.tolist()))],
                list(reversed(tab_uni[tab_uni.mapper=='GeneID'].uniprot.tolist()))))
        if return_not_mapped:
            for (original,mapped) in zip(query_list,list(map(dictio.get,query_list))):
                map_str = str(mapped)
                if map_str == "None":
                    print(f"Couldn't Map {original}")
        return list(map(dictio.get,query_list))

    elif source=='symbol' and target=='ensembl':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((y, x.split('|')[0]) for x, y in list(zip(ense,ncbidb['Symbol'])) if 'ENS' in x)
        if return_not_mapped:
            for (original,mapped) in zip(query_list,list(map(dictio.get,query_list))):
                map_str = str(mapped)
                if map_str == "None":
                    print(f"Couldn't Map {original}")
        return list(map(dictio.get,query_list))

    elif source == 'symbol' and target == 'entrez':
        dictio=dict((y, x) for x, y in list(zip(ncbidb['GeneID'],ncbidb['Symbol'])))
        if return_not_mapped:
            for (original,mapped) in zip(query_list,list(map(dictio.get,query_list))):
                map_str = str(mapped)
                if map_str == "None":
                    print(f"Couldn't Map {original}")
        return list(map(dictio.get,query_list))

    elif source=='symbol' and target=='uniprot':
        dictio=dict(zip(list(reversed(tab_uni[tab_uni.mapper=='Gene_Name'].id.tolist())),
                list(reversed(tab_uni[tab_uni.mapper=='Gene_Name'].uniprot.tolist()))))
        if return_not_mapped:
            for (original,mapped) in zip(query_list,list(map(dictio.get,query_list))):
                map_str = str(mapped)
                if map_str == "None":
                    print(f"Couldn't Map {original}")
        return list(map(dictio.get,query_list))

    elif source=='uniprot' and target=='ensembl':
        dictio=dict(zip(list(reversed(tab_uni[tab_uni.mapper=='Ensembl'].uniprot.tolist())),
                list(reversed(tab_uni[tab_uni.mapper=='Ensembl'].id.tolist()))))
        if return_not_mapped:
            for (original,mapped) in zip(query_list,list(map(dictio.get,query_list))):
                map_str = str(mapped)
                if map_str == "None":
                    print(f"Couldn't Map {original}")
        return list(map(dictio.get,query_list))

    elif source=='uniprot' and target=='entrez':
        dictio=dict(zip(list(reversed(tab_uni[tab_uni.mapper=='GeneID'].uniprot.tolist())),
                [int(g) for g in list(reversed(tab_uni[tab_uni.mapper=='GeneID'].id.tolist()))]))
        if return_not_mapped:
            for (original,mapped) in zip(query_list,list(map(dictio.get,query_list))):
                map_str = str(mapped)
                if map_str == "None":
                    print(f"Couldn't Map {original}")
        return list(map(dictio.get,query_list))

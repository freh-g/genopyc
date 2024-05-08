import requests
import pandas as pd
from genopyc.genomic_features.get_variants_position import *

def get_ov_region(snp=None, ch=None, position=None, window=500, features=['gene'], mode='region'):
    """
    Retrieve overlap regions for a specified SNP or genomic position.

    Parameters:
    - snp (str, optional): The variant ID (SNP) for which overlap regions are requested. Default is None.
    - ch (str or int, optional): Chromosome name or identifier. Required if mode is 'SNP'.
    - position (int, optional): Genomic position (start position) for which overlap regions are requested. Required if mode is 'region'.
    - window (int, optional): Window size to extend the genomic region around the SNP or position. Default is 500.
    - features (list, optional): List of features to include in the retrieval. Default is ['gene']. Possible values are [band, gene, transcript, cds, exon, repeat, simple, misc, variation, somatic_variation, structural_variation, somatic_structural_variation, constrained, regulatory, motif, other_regulatory, array_probe, man].
    - mode (str, optional): Mode of operation, either 'SNP' or 'region'. Default is 'region'.

    Returns:
    - list of pandas DataFrames: A list of pandas DataFrames containing overlap regions for each specified genomic feature type.

    Note:
    This function queries the Ensembl REST API to retrieve overlap regions for the specified SNP or genomic position.
    It returns a list of pandas DataFrames, where each DataFrame contains overlap regions for a specific genomic feature type.
    """

    str_features = ';'.join(['feature=' + x for x in features])

    if mode == 'SNP':
        pos = get_variants_position(snp)
        ch = int(pos[0][1])
        genomic_position = pos[0][2]
        start = genomic_position - window // 2
        stop = genomic_position + window // 2
    else:

        start = position - window // 2
        stop = position + window // 2
    http = "https://rest.ensembl.org/overlap/region/human/%s:%s-%s?%s" % (ch, start, stop, str_features)
    risposta = requests.get(http, headers={"Content-Type": "application/json"})
    if risposta.ok:
        risposta = risposta.json()
        lodfs = []
        list_of_features_retrieved = list(set([e['feature_type'] for e in risposta]))
        for f in features:
            if f not in list_of_features_retrieved:
                print(f"Couldn't retrieve data for feature: {f}")
        for x in list_of_features_retrieved:
            tmp_list = [dict(sorted(e.items())) for e in risposta if e['feature_type'] == x]
            # find the most number of keys dictionary
            length_list = [len(e.items()) for e in tmp_list]
            max_l = max(length_list)
            index_max = length_list.index(max_l)
            cols = tmp_list[index_max].keys()
            e = pd.DataFrame(columns=cols)
            for i, f in enumerate(tmp_list):
                serie = pd.Series(data=f, index=list(f.keys()))
                e.loc[i] = serie

            lodfs.append(e)
        return lodfs
    else:
        print(f'ERROR: Bad Request:\n{risposta.text}')

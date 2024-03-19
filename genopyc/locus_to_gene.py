import requests
import pandas as pd

def OT_L2G(list_of_variants, score=0.1, output='genes'):
    """
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
    """

    query = """{
              genesForVariant(variantId:"%s"){
                overallScore
                gene{id}
              }
            }
            """
    OT_url = 'https://api.genetics.opentargets.org/graphql'

    results = {}
    for variant in list_of_variants:
        r = requests.post(OT_url, json={'query': query % (variant)})
        r = r.json()
        ResultsForVariant = []
        try:
            for data in r['data']['genesForVariant']:
                ResultsForVariant.append((data['gene']['id'], data['overallScore']))
            results[variant] = ResultsForVariant
        except Exception as e:
            print(e,f"Couldn't retrieve data for variant {variant}")
    if output == 'all':
        raw_data = {key: sorted(value, key=lambda x: x[1], reverse=True) for key, value in results.items()}
        cols = ['id', 'gene', 'score']
        data = []
        for k, v in raw_data.items():
            for gene in v:
                data.append((k, gene[0], gene[1]))
        res = pd.DataFrame(data, columns=cols)
        return res
    else:
        return list(set(sum([[value[0] for value in values if value[1] > score] for (key, values) in results.items()], [])))

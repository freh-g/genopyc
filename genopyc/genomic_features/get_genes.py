import requests

def get_genes(ch, position, window_size=10000, pop='EUR', features=['gene'], mode='all'):
    """
    Retrieve genes in a window centered around a genomic position and compute the distance between the position and all genes.

    Parameters:
    - ch (str): Chromosome identifier.
    - position (int): Genomic position around which the window is centered.
    - window_size (int, optional): Size of the window in base pairs. Default is 10,000.
    - pop (str, optional): Population for which to retrieve gene data. Default is 'EUR' (European).
    - features (list, optional): List of features to include in the retrieval. Default is ['gene'].
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
    """

    win_start = position - window_size // 2
    win_end = position + window_size // 2
    str_features = ';'.join(['feature=' + x for x in features])
    http = "https://rest.ensembl.org/overlap/region/human/%s:%s-%s?%s" % (ch, win_start, win_end, str_features)
    response = requests.get(http, headers={"Content-Type": "application/json"})
    if response.ok:
        response = response.json()
        if mode == 'complete_data':
            return response
        elif mode == 'all':
            elements = {}
            for el in response:
                try:
                    elements[el['external_name']] = int(el['start'] - position)
                except:
                    pass
            return elements
        elif mode == 'closest_forward':
            elements = {}
            for el in response:
                try:
                    elements[el['external_name']] = int(el['start'] - position)
                except:
                    pass
            try:
                return min([(k, v) for (k, v) in elements.items() if v > 0], key=lambda x: x[1])
            except:
                return 'no_genes_forward'
        elif mode == 'closest_backward':
            elements = {}
            if el['biotype'] == 'protein_coding':
                try:
                    elements[el['external_name']] = int(el['start'] - position)
                except:
                    pass
            try:
                return max([(k, v) for (k, v) in elements.items() if v < 0], key=lambda x: x[1])
            except:
                return 'no_genes_backward'
        elif mode == 'closest_overall':
            elements = {}
            for el in response:
                try:
                    elements[el['external_name']] = int(el['start'] - position)
                except:
                    pass
            try:
                return min([(k, np.absolute(v)) for (k, v) in elements.items()], key=lambda x: x[1])
            except:
                return 'no_genes'
    else:
        print(f'ERROR: Bad Request:\n{response.text}')

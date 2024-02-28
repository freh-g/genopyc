import get_ov_region

def closest_genes(position_id, chromosome, position, window_size, type_of_gene=False, mode='start'):
    """
    Find the closest upstream and downstream genes to a given genomic position.

    Parameters:
    - position_id (str): Identifier for the genomic position.
    - chromosome (str or int): Chromosome name or identifier.
    - position (int): Genomic position for which closest genes are to be found.
    - window_size (int): Window size to extend the genomic region around the position.
    - type_of_gene (str, optional): Biotype of genes to consider. Default is False (include all genes).
    - mode (str, optional): Mode of calculating distances, either 'start' or 'absolute'. Default is 'start'.
        start:  In this case the distance will be calculated from the beginning of both upstream and downstream genes 
        absolute: Upstream distance is calcuated as distance (base pairs) to the end and downstream as distance to the start of the gene

    Returns:
    - tuple: A tuple containing the position ID, the closest upstream gene ID, and the closest downstream gene ID.

    Note:
    This function utilizes the `get_ov_region` function to retrieve genes within the specified window around the genomic position.
    It calculates distances from the specified position to the start (or end) of genes and finds the closest upstream and downstream genes.
    """
    elements = get_ov_region(chr=chromosome, location=position, features=['gene'], window=window_size)
    all_genes = []
    for i, r in elements[0].iterrows():
        if type_of_gene:
            if r['biotype'] == type_of_gene:
                all_genes.append((r['gene_id'], r['start'] - position, r['end'] - position))
        else:
            all_genes.append((r['gene_id'], r['start'] - position, r['end'] - position))

    positive_r = []
    negative_r = []
    for r in all_genes:
        # Look for the closest upstream gene 
        if r[1] > 0:
            positive_r.append(r)
        elif r[2] < 0:
            negative_r.append(r)

    if mode == 'start':
        if len(negative_r) > 0:
            closest_downstream_gene = max(negative_r, key=lambda x: x[1])[0]
        else:
            print(f'no downstream genes in the selected window for variant {position_id}')
            closest_downstream_gene = ''

        if len(positive_r) > 0:
            closest_upstream_gene = min(positive_r, key=lambda x: x[1])[0]
        else:
            print(f'no upstream genes in the selected window for variant {position_id}')
            closest_upstream_gene = ''

    elif mode == 'absolute':
        """ In this case the distance will be calculated from the beginning of the downstream and the end of the upstream """
        if len(negative_r) > 0:
            closest_downstream_gene = max(negative_r, key=lambda x: x[2])[0]
        else:
            print(f'no downstream genes in the selected window for variant {position_id}')
            closest_downstream_gene = ''

        if len(positive_r) > 0:
            closest_upstream_gene = min(positive_r, key=lambda x: x[2])[0]
        else:
            print(f'no upstream genes in the selected window for variant {position_id}')
            closest_upstream_gene = ''

    return (position_id, closest_upstream_gene, closest_downstream_gene)

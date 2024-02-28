def get_eqtl_variant(rsid):
    """
    Retrieve eQTL (expression quantitative trait loci) associations for a given variant.

    Parameters:
    - rsid (str): The variant ID (e.g., rsID) for which eQTL associations are requested.

    Returns:
    - dict: A dictionary containing eQTL associations for the specified variant.

    Note:
    This function queries the EBI eQTL Catalog API to retrieve eQTL associations for the specified variant.
    It returns a dictionary containing eQTL association information, including target genes and other relevant data, for the specified variant.
    """

    url = 'http://www.ebi.ac.uk/eqtl/api/associations/%s' % (rsid)
    risp = requests.get(url)
    if risp.ok:
        risp = risp.json()
        return risp
    else:
        print(f'ERROR: Bad Request:\n{risp.text}')

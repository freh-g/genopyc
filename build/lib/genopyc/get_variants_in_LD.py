import requests

def get_variants_in_LD(variant, r2, pop='EUR'):
    """
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
    """
    
    http = "https://rest.ensembl.org/ld/human/%s/1000GENOMES:phase_3:%s?r2=%s" % (variant, pop, r2)
    
    try:
        resp = requests.get(http, headers={ "Content-Type": "application/json"})
        if resp.ok:
            variants = resp.json()
            return [x['variation2'] for x in variants if float(x['r2']) >= r2]
        else:
            print(f'ERROR: Bad Request:\n{resp.text}')
    except Exception as e:
        print("Error:", e)
        return []


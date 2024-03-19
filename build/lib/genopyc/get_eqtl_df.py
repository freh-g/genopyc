import pandas as pd
import requests
import os
import pickle
import numpy as np

def get_eqtl_df(rsid, p_value=0.005, increase_index=False):
    """
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
    """

    location = os.path.dirname(os.path.realpath(__file__))
    out = os.path.join(location, 'data')
    with open(out + '/uberon_dict.pickle', 'rb') as f:
        ubdict = pickle.load(f)
    url = 'http://www.ebi.ac.uk/eqtl/api/v1/associations/%s?size=1000' % (rsid)
    response = requests.get(url)
    if response.ok:
        eqtls = response.json()
        try:
            eqtl_df = pd.DataFrame(columns=['variantid', 'p_value', 'log_pval', 'beta', 'alt', 'gene_id', 'tissue', 'study_id'])
            for ass in eqtls['_embedded']['associations'].keys():
                pval = eqtls['_embedded']['associations'][ass]['pvalue']
                nlog_pval = -np.log10(pval)
                beta = eqtls['_embedded']['associations'][ass]['beta']
                alt = eqtls['_embedded']['associations'][ass]['alt']
                geneid = eqtls['_embedded']['associations'][ass]['gene_id']
                tissue = eqtls['_embedded']['associations'][ass]['tissue']
                study = eqtls['_embedded']['associations'][ass]['study_id']
                eqtl_df.loc[ass] = [rsid, pval, nlog_pval, beta, alt, geneid, tissue, study]

            eqtl_df = eqtl_df.loc[eqtl_df.p_value <= p_value]
            eqtl_df.tissue = eqtl_df.tissue.apply(lambda x: x.replace('UBER_', 'UBERON_'))
            eqtl_df['tissue_name'] = list(map(ubdict.get, eqtl_df.tissue.tolist()))

            eqtl_df = eqtl_df.reset_index(drop=True)
            if increase_index:
                eqtl_df.index += 1
        except Exception as er:
            print(er, eqtls)
            return None
        return eqtl_df
    else:
        print(f'ERROR: Bad Request:\n{response.text}')

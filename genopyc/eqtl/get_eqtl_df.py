import pandas as pd
import requests


def get_eqtl_df(rsid,start=0,size=100 ,p_value=0.005, increase_index=False):
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
    url = 'http://www.ebi.ac.uk/eqtl/api/v3/associations?rsid=%s&start=%i&size=%i' % (rsid,start,size)
    response = requests.get(url)
    if (response.ok) & (response.text!='[]'):
        try:
            res_df = pd.DataFrame(response.json())
            dataset_ids = set(res_df.dataset_id.tolist())
            datasets_metadata = [requests.get(f'http://www.ebi.ac.uk/eqtl/api/v2/datasets/{study}').json() for study in dataset_ids]
            mapping_dict = dict(zip([d['dataset_id'] for d in datasets_metadata],[d['tissue_label'] for d in datasets_metadata]))
            list_of_studies = res_df['dataset_id'].tolist()
            res_df['tissue_name'] = list(map(mapping_dict.get,list_of_studies))
            res_df['beta'] = res_df['beta'].astype(float)
            res_df['pvalue'] = res_df['pvalue'].astype(float)
            return res_df
        
        except Exception as er:
            print(er)
            return None
    else:
        print(f'ERROR: Bad Request:\n{response.text}')
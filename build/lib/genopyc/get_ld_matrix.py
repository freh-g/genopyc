import requests
import pandas as pd

def get_ld_matrix(list_of_snps, token, pop='EUR', metric='r2'):
    """
    Get the linkage disequilibrium (LD) matrix for a list of SNPs.

    Parameters:
    - list_of_snps (list): A list of SNP IDs.
    - token (str): The authentication token for accessing the LDlink API.
    - pop (str, optional): The population for which LD information is requested. Default is 'EUR' (European population).
    - metric (str, optional): The LD metric to be used. Default is 'r2'.

    Returns:
    - pandas DataFrame: A DataFrame representing the LD matrix. Rows and columns are SNP IDs, and cell values represent the LD metric between corresponding SNPs.

    Note:
    This function queries the LDlink API to retrieve the LD matrix for the specified list of SNPs, for this reason it needs the user to be registered on LDlink in order to obtain a token.
    It returns a pandas DataFrame containing the LD matrix, where rows and columns are SNP IDs, and cell values represent the LD metric.
    If any LD metric values are missing (NA), they are replaced with None, and the DataFrame is converted to numeric data types.
    """
    snp_string = '\n'.join(list_of_snps)

    headers = {
        'Content-Type': 'application/json',
    }

    params = (
        ('token', token),
    )

    json_data = {
        'snps': snp_string,
        'pop': pop,
        'r2_d': metric,
        'genome_build': 'grch38',
    }

    response = requests.post('https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix', headers=headers, params=params, json=json_data, verify=False)
    if response.ok:
        response = response.text
        if "error" in response:
            print("error : Input variant list does not contain any valid RS numbers or coordinates.")
            return None
        else:
            data_frame = pd.DataFrame([x.split('\t') for x in response.split('\n')])
            new_header = data_frame.iloc[0]
            data_frame = data_frame[1:]  # take the data less the header row
            data_frame.columns = new_header  # set the header row as the df header

            new_rows = data_frame[data_frame.columns[0]]
            data_frame = data_frame[data_frame.columns[1:]].set_index(new_rows)

            data_frame.replace('NA', None, inplace=True)
            data_frame = data_frame.astype(None)

            return data_frame.fillna(0).iloc[:-1]
    else:
        print(f'ERROR: Bad Request:\n{response.text}')
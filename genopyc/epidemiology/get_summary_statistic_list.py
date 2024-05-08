import requests

def get_summary_statistic_list():
    """
    Retrieve a list of summary statistics for all available GWAS studies.

    Returns:
    - list: A list of dictionaries containing summary statistics for all available GWAS studies.

    Note:
    This function queries the EBI GWAS Catalog Summary Statistics API to retrieve summary statistics for all available GWAS studies.
    It returns a list of dictionaries, where each dictionary contains summary statistics information for a specific GWAS study.
    """

    http = 'https://www.ebi.ac.uk/gwas/summary-statistics/api/associations'
    response = requests.get(http)
    if response.ok:
        response = response.json()
        return response()
    else:
        print(f'ERROR: Bad Request:\n{response.text}')

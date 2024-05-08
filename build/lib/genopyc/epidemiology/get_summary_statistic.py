import requests

def get_summary_statistic(study):
    """
    Retrieve summary statistics for a given GWAS study.

    Parameters:
    - study (str): The study identifier for which summary statistics are requested.

    Returns:
    - dict: A dictionary containing summary statistics for the specified GWAS study.

    Note:
    This function queries the EBI GWAS Catalog Summary Statistics API to retrieve summary statistics for the specified study.
    It returns a dictionary containing summary statistics information, including associations, for the specified study.
    """

    http = 'https://www.ebi.ac.uk/gwas/summary-statistics/api/studies/%s/associations' % (study)
    response = requests.get(http)
    if response.ok:
        response = response.json()
        return response()
    else:
        print(f'ERROR: Bad Request:\n{response.text}')

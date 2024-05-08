import requests
import pandas as pd

def pairwise_ld(ch, start, end, pop='EUR'):
    """
    Retrieve pairwise linkage disequilibrium (LD) information for variants within a genomic region.

    Parameters:
    - ch (str): Chromosome name or identifier.
    - start (int): Start position of the genomic region.
    - end (int): End position of the genomic region.
    - pop (str, optional): The population for which LD information is requested. Default is 'EUR' (European population).

    Returns:
    - pandas DataFrame: A DataFrame containing pairwise LD information between variants within the specified genomic region.
      Columns include 'v1' and 'v2' for variant IDs and 'r2' for the LD coefficient.

    Note:
    This function queries the Ensembl REST API to retrieve pairwise LD information for variants within the specified genomic region.
    It returns a DataFrame containing pairwise LD information, where 'v1' and 'v2' represent variant IDs and 'r2' represents the LD coefficient.
    """

    http = f"https://rest.ensembl.org/ld/human/region/{ch}:{str(start)}..{str(end)}/1000GENOMES:phase_3:{pop}"
    response = requests.get(http, headers={"Content-Type": "application/json"})
    if response.ok:
        response = response.json()
        if response != []:
            ld_mat = []

            for i, element in enumerate(response):
                try:
                    v1 = element['variation1']
                    v2 = element['variation2']
                    r2 = element['r2']
                    ld_mat.append((v1, v2, r2))
                except Exception as r:
                    print(r, f'- Error for variant "{element}"')
                ld_mat_df = pd.DataFrame(ld_mat, columns=['v1', 'v2', 'r2'])

            return ld_mat_df.sort_values(by='v1')
        else:
            print(f"No data relative to the genomic location {ch}:{start}-{end}")
    else:
        print(f'ERROR: Bad Request:\n{response.text}')
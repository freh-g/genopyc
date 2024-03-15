import requests

def get_sequence(ch, start, stop):
    """
    Retrieve the genomic sequence for a given genomic region.

    Parameters:
    - ch (str or int): Chromosome name or identifier.
    - start (int): Start position of the genomic region.
    - stop (int): End position of the genomic region.

    Returns:
    - str: The genomic sequence for the specified genomic region.

    Note:
    This function queries the Ensembl REST API to retrieve the genomic sequence for the specified genomic region.
    It returns the genomic sequence as a string.
    """
    http = "https://rest.ensembl.org/sequence/region/human/%s:%s..%s?" % (ch, start, stop)
    risposta = requests.get(http, headers={"Content-Type": "application/json"})
    if risposta.ok:
        risposta = risposta.json()
        return risposta['seq']
    else:
        print(f'ERROR: Bad Request:\n{risposta.text}')
        return None

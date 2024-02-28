import requests

def grch_liftover(chromosome, start, end, source, target):
    """
    Lift over genomic coordinates from one assembly version to another.

    Parameters:
    - chromosome (str or int): Chromosome name or identifier.
    - start (int): Start position of the genomic region.
    - end (int): End position of the genomic region.
    - source (str): Source genome assembly version.
    - target (str): Target genome assembly version.

    Returns:
    - tuple or None: A tuple containing the chromosome, start, and end positions in the target assembly if successful, 
      otherwise returns None.

    Note:
    This function queries the Ensembl REST API to perform coordinate liftover from the source genome assembly to the target assembly.
    It returns a tuple containing the chromosome, start, and end positions in the target assembly if liftover is successful,
    otherwise returns None.
    """
    url = "https://rest.ensembl.org/map/human/%s/%s:%i..%i/%s?" % (source, chromosome, start, end, target)
    r = requests.get(url, headers={"Content-Type": "application/json"}).json()
    try:
        return (chromosome, r['mappings'][0]['mapped']['start'], r['mappings'][0]['mapped']['end'])
    except:
        return None

import requests

def get_phenotypes(chromosome, start, stop, feature_type='Genes', only_phenotypes=1):
    """
    Retrieve phenotype annotations annotations of a given genomic region.

    Parameters:
    - chromosome (str): Chromosome name or identifier.
    - start (int): Start position of the genomic region.
    - stop (int): End position of the genomic region.
    - feature_type (str, optional): Type of genomic feature to include in annotations. Default is 'Genes'.
    - only_phenotypes (int, optional): Whether to include only phenotypes. Default is 1 (True).

    Returns:
    - dict: A dictionary containing annotations of the specified genomic region.

    Note:
    This function queries the Ensembl REST API to retrieve annotations for the specified genomic region.
    It returns a dictionary containing annotations, including phenotypes, for the specified genomic region.
    """

    http = "https://rest.ensembl.org/phenotype/region/homo_sapiens/%s:%s-%s?only_phenotypes=%s;feature_type=%s" % (
    chromosome, start, stop, only_phenotypes, feature_type)
    annot = requests.get(http, headers={"Content-Type": "application/json"})
    return annot.json()

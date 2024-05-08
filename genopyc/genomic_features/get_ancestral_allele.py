from genopyc.genomic_features.get_variants_info import *

def get_ancestral_allele(variants):
    """
    Retrieve the ancestral allele for the specified variants.
    Parameters:
    - variants (list or str): A list of variant IDs or a single variant ID as a string.

    Returns:
    - dict: A dictionary where keys are variant IDs and values are their corresponding ancestral alleles.

    Note:
    This function relies on the get_variants_info function to retrieve variant information.
    It extracts the ancestral alleles from the retrieved variant mappings based on the specified mode. It calculates the most frequently occurring allele among mappings.
    """
    
    variant_info = get_variants_info(variants)
    variant_mappings = dict(zip(variant_info.keys(),[info['mappings'] for info in variant_info.values()]))
    
    def most_frequent(List):
        return max(set(List), key = List.count)

    variant_ancestral_alleles = {}
    for rsid, mappings in variant_mappings.items():
        if len(mappings) > 1:
            list_of_alleles = []
            for element in mappings:
                ref_allele = element['allele_string'].split('/')[0]
                list_of_alleles.append(ref_allele)        
            try:
                variant_ancestral_alleles[rsid] = most_frequent(list_of_alleles)
            except Exception as error:
                print(f'Error for rsid {rsid}: ' + repr(error))
        else:
            try:
                ref_allele = mappings[0]['allele_string'].split('/')[0]
                variant_ancestral_alleles[rsid] = mappings[0]['ancestral_allele']
            except Exception as error:
                print(f'Error for rsid {rsid}: ' + repr(error) + ' MAPPINGS = ' + str(mappings))

    variant_ancestral_alleles = {k:v for k,v in variant_ancestral_alleles.items() if v != None}
    return variant_ancestral_alleles


import requests

def variantId_mapping(list_of_variants, source='variantid', target='rsid'):
    """
    Convert variants from one identifier type to another using the Open Targets Genetics API.

    Parameters:
    - list_of_variants (list): List of variant identifiers to be converted.
    - source (str, optional): Identifier type of the variants in list_of_variants. Default is 'variantid'.
    - target (str, optional): Identifier type to which the variants will be converted. Default is 'rsid'.

    Returns:
    - list: List of converted variant identifiers.

    Note:
    - Currently supports conversion between 'variantid' and 'rsid'.
    - If a variant cannot be converted, it will be skipped, and a message will be printed.
    """

    OT_url = 'https://api.genetics.opentargets.org/graphql'
    MappingDict = {}
    if (source == 'variantid') & (target == 'rsid'):
        query = """{
               variantInfo(variantId:"%s"){
                   rsId
                   }
              }"""
        for variant in list_of_variants:
            r = requests.post(OT_url, json={'query': query % (variant)})
            json_response = r.json()
            MappingDict[variant] = json_response['data']['variantInfo']['rsId']
        return list(map(MappingDict.get, list_of_variants))
    elif (source == 'rsid') & (target == 'variantid'):
        query = """{
              search(queryString:"%s"){
                  variants{
                      id
                      }
                   }
              }
              """
        for variant in list_of_variants:
            try:
                r = requests.post(OT_url, json={'query': query % (variant)})
                json_response = r.json()
                MappingDict[variant] = json_response['data']['search']['variants'][0]['id']
            except Exception as e:
                print(e, f"Couldn't Convert Variant {variant}")
        return list(map(MappingDict.get, list_of_variants))


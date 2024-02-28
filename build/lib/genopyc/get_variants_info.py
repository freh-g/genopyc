import requests

def get_variants_info(idlist, chunked=False, chunksize=200):
    """
    Retrieve information about variants in the Homo sapiens species from the Ensembl variation database.

    Parameters:
    - idlist (list or str): A list of variant IDs or a single variant ID as a string.
    - chunked (bool, optional): Whether to chunk the requests into smaller portions if the number of IDs exceeds 200. Default is False.
    - chunksize (int, optional): The size of each chunk when chunking requests. Default is 200.

    Returns:
    - dict: A dictionary containing information about the variants, where keys are variant IDs and values are variant information.

    Note:
    This function requires internet connectivity to access the Ensembl REST API.
    Ensembl accepts a maximum of 200 IDs per request. If the number of IDs exceeds 200 and chunking is not disabled,
    the requests will be split into smaller chunks to retrieve variant information.
    """
    
    if not isinstance(idlist, list):
        idlist = [idlist]
    http = "https://rest.ensembl.org/variation/homo_sapiens"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    chunked_idlist = []
    if chunked or len(idlist) > 200:
        print('total number of chunks: %s' % (int(len(idlist) / chunksize) + 1))
        for i in range(0, len(idlist), chunksize):
            chunked_idlist.append(idlist[i:i + chunksize])
        results = {}
        for i, chunk in enumerate(chunked_idlist):
            response = requests.post(http, headers=headers,
                                     data="{" + '"ids" : {}'.format(str(chunk).replace("'", '"')) + "}")
            if response.ok:
                results.update(response.json())
                print('chunk %s processed' % (i))
            else:
                print(f'ERROR: Bad Request:\n{response.text}')
        return results

    else:
        response = requests.post(http, headers=headers,
                                 data="{" + '"ids" : {}'.format(str(idlist).replace("'", '"')) + "}")
        if response.ok:
            return response.json()
        else:
            print(f'ERROR: Bad Request:\n{response.text}')

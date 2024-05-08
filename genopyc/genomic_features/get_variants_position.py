import requests

def get_variants_position(idlist, chunked=False, chunksize=200):
    """
    Retrieve the positions (chromosome and position) of variant IDs from the Ensembl database for Homo sapiens.

    Parameters:
    - idlist (list or str): A list of variant IDs or a single variant ID as a string.
    - chunked (bool, optional): Whether to chunk the requests into smaller portions. Default is False.
    - chunksize (int, optional): The size of each chunk when chunking requests. Default is 200.

    Returns:
    - list of tuples: A list of tuples containing the variant ID, chromosome, and position.
    
    Note:
    This function requires internet connectivity to access the Ensembl REST API.
    """

    if not isinstance(idlist, list):
        idlist = [idlist]

    http = "https://rest.ensembl.org/variation/homo_sapiens"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    if chunked or len(idlist) > 200:
        chunked_idlist = []
        print('total number of chunks: %s' % (int(len(idlist) / chunksize) + 1))
        for i in range(0, len(idlist), chunksize):
            chunked_idlist.append(idlist[i:i + chunksize])
        results = []
        for i, chunk in enumerate(chunked_idlist):
            response = requests.post(http, headers=headers,
                                     data="{" + '"ids" : {}'.format(str(chunk).replace("'", '"')) + "}")
            if response.ok:
                response = response.json()
                for key, value in response.items():
                    try:
                        chr = value['mappings'][0]['location'].split(':')[0]
                        pos = value['mappings'][0]['start']
                        results.append((key, chr, pos))

                    except:
                        print(f"Couldn't Retrieve Position for variant {key}")
                        pass
                print(f"chunk {i} processed")
            else:
                print(f'ERROR: Bad Request:\n{response.text}')
        return results
    else:
        response = requests.post(http, headers=headers,
                                 data="{" + '"ids" : {}'.format(str(idlist).replace("'", '"')) + "}")
        if response.ok:
            response = response.json()
            results = []
            for key, value in response.items():
                try:
                    chr = value['mappings'][0]['location'].split(':')[0]
                    pos = value['mappings'][0]['start']
                    results.append((key, chr, pos))

                except:
                    print(f"Couldn't Retrieve Position for variant {key}")
                    pass
            return results
        else:
            print(f'ERROR: Bad Request:\n{response.text}')

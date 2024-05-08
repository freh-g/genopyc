import requests
def get_genes_position(idlist, chunked=False, chunksize=200):
    """
    Retrieve the coordinates of genes in the GRCh38 assembly for a list of Ensembl IDs.

    Parameters:
    - idlist (list or str): A list of Ensembl gene IDs or a single Ensembl gene ID as a string.
    - chunked (bool, optional): Whether to chunk the requests into smaller portions if the number of IDs exceeds 200. Default is False.
    - chunksize (int, optional): The size of each chunk when chunking requests. Default is 200.

    Returns:
    - list of tuples: A list of tuples containing the Ensembl gene ID, chromosome number, start position, and end position.
    """
    if not isinstance(idlist, list):
        idlist = [idlist]

    http = "https://rest.ensembl.org/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    if chunked or len(idlist) > 200:
        chunked_idlist = []
        print(f'Total number of chunks: {int(len(idlist) / chunksize) + 1}')
        for i in range(0, len(idlist), chunksize):
            chunked_idlist.append(idlist[i:i + chunksize])
        results = []
        for i, chunk in enumerate(chunked_idlist):
            response = requests.post(http, headers=headers, data="{" + '"ids" : {}'.format(str(chunk).replace("'", '"')) + "}")
            if response.ok:
                response = response.json()
                ListOfTuples = []
                for k, v in response.items():
                    try:
                        ListOfTuples.append((k, int(v['seq_region_name']), v['start'], v['end']))
                    except:
                        print(f"Couldn't retrieve position for gene {k}")
                        continue
                results.append(ListOfTuples)
                print(f'Chunk {i + 1} done')
            else:
                print(f'ERROR: Bad Request:\n{response.text}')
        return results
    else:
        ListOfTuples = []
        for gene in idlist:
            response = requests.get(f"{http}/{gene}", headers=headers)
            if response.ok:
                response = response.json()
                try:
                    ListOfTuples.append((gene, int(response['seq_region_name']), response['start'], response['end']))
                except:
                    print(f"Couldn't retrieve position for gene {gene}")
                    continue
            else:
                print(f'ERROR: Bad Request:\n{response.text}')       
        return ListOfTuples

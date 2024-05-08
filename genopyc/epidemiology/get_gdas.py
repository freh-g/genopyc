import pandas as pd
import requests

def get_gdas(query_list, username, password, mode):
    """
    Retrieve gene or variant-disease associations from the DisGeNET database.

    Parameters:
    - query_list (list): A list of gene symbols or variant IDs for which associations are requested.
    - username (str): Username for accessing the DisGeNET API.
    - password (str): Password for accessing the DisGeNET API.
    - mode (str): Mode of query, either 'genes' or 'variants'.

    Returns:
    - pandas DataFrame: A DataFrame containing gene or variant-disease associations.

    Note:
    This function queries the DisGeNET API to retrieve gene or variant-disease associations.
    It requires authentication using a username and password provided by the user.
    The mode parameter specifies whether the query is for genes or variants.
    """
    
    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    auth_params = {"email": username, "password": password}
    api_host = 'https://www.disgenet.org/api'
    req = requests.Session()
    url = api_host + '/auth/'
    response = req.post(url, data=auth_params)
    token = response.json()['token']
    req.headers.update({"Authorization": "Bearer %s" % token}) 
       

    if mode == "genes":
        if len(query_list) > 100:
            list_dfs = []
            chunks_query = list(chunks(query_list, 100))
            for c in chunks_query:
                query_str = "%2C".join(c)
                resp = req.get(api_host + '/gda/gene/{}'.format(query_str)).json()
                df = pd.DataFrame(resp)
                list_dfs.append(df)
            
            df = pd.concat(list_dfs)    
        else:
            query_str = "%2C".join(query_list)
            resp = req.get(api_host + '/gda/gene/{}'.format(query_str)).json()
            df = pd.DataFrame(resp)
    elif mode == "variants":
        if len(query_list) > 100:
            list_dfs = []
            chunks_query = list(chunks(query_list, 100))
            for c in chunks_query:
                query_str = "%2C".join(c)
                resp = req.get(api_host + '/vda/gene/{}'.format(query_str)).json()
                df = pd.DataFrame(resp)
                list_dfs.append(df)
            
            df = pd.concat(list_dfs)    
        else:
            query_str = "%2C".join(query_list)
            resp = req.get(api_host + '/vda/variant/{}'.format(query_str)).json()
            df = pd.DataFrame(resp)
    return df

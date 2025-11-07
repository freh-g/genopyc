# import requests
# import pandas as pd

# def OT_L2G(list_of_variants, score=0.1, output='genes'):
#     """
#     Retrieve genes associated with variants from the Open Targets Genetics API.

#     Parameters:
#     - list_of_variants (list): List of variant identifiers.
#     - score (float, optional): Threshold score for gene association. Genes with scores higher than this threshold will be included. Default is 0.1.
#     - output (str, optional): Output format. 'genes' returns a list of genes, 'all' returns a DataFrame with variant, gene, and score information. Default is 'genes'.

#     Returns:
#     - list or DataFrame: Depending on the output parameter, either a list of genes or a DataFrame with variant, gene, and score information.

#     Note:
#     - The function retrieves genes associated with variants from the Open Targets Genetics API.
#     - The 'score' parameter filters genes based on their association score with the variants.
#     - 'output' can be set to 'genes' to return a list of genes or 'all' to return detailed information in a DataFrame.
#     """

#     query = """{
#               genesForVariant(variantId:"%s"){
#                 overallScore
#                 gene{id}
#               }
#             }
#             """
#     OT_url = 'https://api.genetics.opentargets.org/graphql'

#     results = {}
#     for variant in list_of_variants:
#         r = requests.post(OT_url, json={'query': query % (variant)})
#         r = r.json()
#         ResultsForVariant = []
#         try:
#             for data in r['data']['genesForVariant']:
#                 ResultsForVariant.append((data['gene']['id'], data['overallScore']))
#             results[variant] = ResultsForVariant
#         except Exception as e:
#             print(e,f"Couldn't retrieve data for variant {variant}")
#     if output == 'all':
#         raw_data = {key: sorted(value, key=lambda x: x[1], reverse=True) for key, value in results.items()}
#         cols = ['id', 'gene', 'score']
#         data = []
#         for k, v in raw_data.items():
#             for gene in v:
#                 data.append((k, gene[0], gene[1]))
#         res = pd.DataFrame(data, columns=cols)
#         return res
#     else:
#         return list(set(sum([[value[0] for value in values if value[1] > score] for (key, values) in results.items()], [])))




# import requests
# import pandas sa pd
# import json

# def def OT_L2G(list_of_variants,size=100):

#     list_of_dataframes = []
#     # --- GraphQL endpoint ---
#     GRAPHQL_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"

#     # --- Minimal query (only L2G predictions) ---
#     query = """
#     query GWASCredibleSetsL2GOnly($variantId: String!, $size: Int!, $index: Int!) {
#     variant(variantId: $variantId) {
#         gwasCredibleSets: credibleSets(
#         studyTypes: [gwas],
#         page: { size: $size, index: $index }
#         ) {
#         rows {
#             l2GPredictions(page: { size: 100, index: 0 }) {
#             rows {
#                 target {
#                 id
#                 approvedSymbol
#                 }
#                 score
#             }
#             }
#         }
#         }
#     }
#     }
#     """
#     list_of_res = []
#     for variant in list_of_variants:

#         # --- Variables ---
#         variables = {
#             "variantId": variant,  # change to your variant
#             "size": size,
#             "index": 0
#         }

#         # --- Headers (add token if needed) ---
#         headers = {"Content-Type": "application/json"}

#         # --- Send request ---
#         data_json = requests.post(
#             GRAPHQL_ENDPOINT,
#             json={"query": query, "variables": variables},
#             headers=headers
#         ).json()
        
#         rows = []
#         for i, cs in enumerate(data_json['data']['variant']['gwasCredibleSets']['rows']):
#             for pred in cs['l2GPredictions']['rows']:
#                 rows.append({
#                     'row_index': i,
#                     'gene': pred['target']['approvedSymbol'],
#                     'score': pred['score']
#                 })

#         df = pd.DataFrame(rows)
#         df = df.groupby('gene').agg({'score':'mean'}) 
#         list_of_dataframes.append(df)
#         list_of_dataframes['variant'] = variant
#     final_df = pd.concat(list_of_dataframes, axis=0)
#     return final_df
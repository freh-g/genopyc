import pandas as pd
import requests

def OT_L2G(lov,size,verbose=False):
    GRAPHQL_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
    query = """
        query GWASCredibleSetsL2GOnly($variantId: String!, $size: Int!, $index: Int!) {
        variant(variantId: $variantId) {
            gwasCredibleSets: credibleSets(
            studyTypes: [gwas],
            page: { size: $size, index: $index }
            ) {
            rows {
                l2GPredictions(page: { size: 100, index: 0 }) {
                rows {
                    target {
                    id
                    approvedSymbol
                    }
                    score
                }
                }
            }
            }
        }
        }
        """
    results = []
    for variant in lov:
        try:
            # --- Variables ---
            variables = {
                "variantId": variant,  # change to your variant
                "size": size,
                "index": 0
            }

            # --- Headers (add token if needed) ---
            headers = {"Content-Type": "application/json"}
            # headers["Authorization"] = "Bearer YOUR_TOKEN"

            # --- Send request ---
            data_json = requests.post(
                GRAPHQL_ENDPOINT,
                json={"query": query, "variables": variables},
                headers=headers
            ).json()
            rows = []
            # Create the dataframe
            for i, cs in enumerate(data_json['data']['variant']['gwasCredibleSets']['rows']):
                for pred in cs['l2GPredictions']['rows']:
                    rows.append({
                        'row_index': i,
                        'gene': pred['target']['approvedSymbol'],
                        'score': pred['score']
                    })

            df = pd.DataFrame(rows)
            df = df.groupby('gene').agg({'score':'mean'})
            df['variant']=variant
            df = df.sort_values(by='score', ascending=False)
            results.append(df)
        except Exception as e:
            if verbose:
                print(f"Couldn't retrieve data for {variant}: {e}")
            continue
    re = pd.concat(results).reset_index()
    return re
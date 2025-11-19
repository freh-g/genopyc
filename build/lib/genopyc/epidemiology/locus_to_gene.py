import pandas as pd
import requests

def OT_L2G(lov, size, verbose=False):
    """
    Retrieve Open Targets L2G (Locus-to-Gene) scores for a list of variants.

    This function queries the Open Targets GraphQL API and fetches L2G predictions 
    from GWAS credible sets for each variant in the input list. For every variant, 
    all predicted target genes are collected, the mean score per gene is computed, 
    and results are returned as a combined dataframe.

    Parameters
    ----------
    lov : list of str
        List of variant IDs to query (e.g., ['1_55516888_T_C', ...]).
    size : int
        Page size for querying credible sets.
    verbose : bool, optional
        If True, prints warnings when a variant fails to retrieve.
        Default is False.

    Returns
    -------
    pandas.DataFrame
        A dataframe with columns:
        - 'gene' : gene symbol
        - 'score' : mean L2G score across credible sets
        - 'variant' : original variant ID

        Rows are sorted by descending score within each variant.

    Notes
    -----
    - This function uses the Open Targets API v4 GraphQL endpoint.
    - API rate limits may apply.
    - If the API response structure changes, the parser may need updates.

    Examples
    --------
    >>> OT_L2G(["1_55516888_T_C"], size=20)
    gene   score    variant
    BRCA1  0.72     1_55516888_T_C
    ...
    """
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
            variables = {
                "variantId": variant,
                "size": size,
                "index": 0
            }

            headers = {"Content-Type": "application/json"}

            data_json = requests.post(
                GRAPHQL_ENDPOINT,
                json={"query": query, "variables": variables},
                headers=headers
            ).json()

            rows = []
            for i, cs in enumerate(data_json['data']['variant']['gwasCredibleSets']['rows']):
                for pred in cs['l2GPredictions']['rows']:
                    rows.append({
                        'row_index': i,
                        'gene': pred['target']['approvedSymbol'],
                        'score': pred['score']
                    })

            df = pd.DataFrame(rows)
            df = df.groupby('gene').agg({'score':'mean'})
            df['variant'] = variant
            df = df.sort_values(by='score', ascending=False)
            results.append(df)

        except Exception as e:
            if verbose:
                print(f"Couldn't retrieve data for {variant}: {e}")
            continue

    re = pd.concat(results).reset_index()
    return re

import pandas as pd
import requests
import numpy as np

def get_associations(efotrait, verbose=False, studyid=False):
    """
    Retrieve SNPs associated with an EFO trait.

    Parameters:
    - efotrait (str): EFO id of the trait for which variants are to be retrieved (ex. EFO_0004994)
    - verbose (bool, optional): If True, the function returns the steps of variant retrieval
    - studyid (bool, optional): If True, the Id of the GWAS where the association was found is retrieved

    Returns:
    - Pandas DataFrame: Variants associated with the trait

    """
    df = pd.DataFrame(columns=['variantid', 'p-value', 'risk_allele', 'RAF', 'beta', 'CI', 'mapped_gene'])
    http = 'https://www.ebi.ac.uk/gwas/rest/api/efoTraits/%s/associations' % (efotrait)
    if verbose:
        print(f"Querying associations for {efotrait}...\n")
    resp = requests.get(http)
    if resp.ok:
        associ = resp.json()
        if verbose:
            print('Building the dataframe...')
        for i, element in enumerate(associ['_embedded']['associations']):
            try:
                variantid = ''.join(element['loci'][0]['strongestRiskAlleles'][0]['riskAlleleName'].split('-')[0:1])
                df.at[i, 'variantid'] = variantid
                df.at[i, 'risk_allele'] = element['loci'][0]['strongestRiskAlleles'][0]['riskAlleleName'].split('-')[-1]
                df.at[i, 'mapped_gene'] = ' '.join([str(elem).strip() for elem in [e['geneName'] for e in element['loci'][0]['authorReportedGenes']]])
                df.at[i, 'p-value'] = float(element['pvalueMantissa']) * 10 ** int(element['pvalueExponent'])

                try:
                    df.at[i, 'RAF'] = float(element['riskFrequency'])
                except:
                    df.at[i, 'RAF'] = np.nan

                df.at[i, 'beta'] = [float(element['betaNum']) if type(element['betaNum']) == float else np.nan][0]
                df.at[i, 'SE'] = element['standardError']
                df.at[i, 'CI'] = element['range']
                try:
                    study_link = element['_links']['study']['href']
                    df.at[i, 'study_url'] = study_link
                except Exception as e:
                    df.at[i, 'study_url'] = np.nan
                    print(e, f"Couldn't retrieve studylink for variant {variantid}")
                if studyid:
                    try:
                        jsonres = requests.get(study_link).json()
                        studyid = jsonres['accessionId']
                        df.at[i, 'studyid'] = studyid
                    except Exception as e:
                        df.at[i, 'studyid'] = np.nan
                        print(e, f"Couldn't retrieve studyID for variant {variantid}")
            except Exception as e:
                print(f'Error {e} for element {element}')
                pass
        df.fillna(np.nan, method=None, inplace=True)
        df['p-value'] = df['p-value'].map("{:.1e}".format)
    else:
        print(f'ERROR: Bad Request:\n{resp.text}')
        return None
    return df

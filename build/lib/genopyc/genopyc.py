import requests
import pandas as pd
import biomapy as bp
import pickle
import os
import ast
from collections import Counter
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib
import requests
from gprofiler import GProfiler
import igraph as ig 
import warnings
warnings.filterwarnings('ignore')
from matplotlib.legend_handler import HandlerBase
from matplotlib.text import Text
from matplotlib.legend import Legend
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input,Output
import dash_cytoscape as cyto
import pandas as pd
import requests
import os
import re
import wget


location = os.path.dirname(os.path.realpath(__file__))
ncbidb = pd.read_csv(os.path.join(location, 'data', 'Homo_sapiens.gene_info.gz'),sep='\t')
tab_uni =  pd.read_csv(os.path.join(location, 'data', 'HUMAN_9606_idmapping.dat.gz'),sep='\t',names=['uniprot','mapper','id'])


def update_mapping_datasets():
    
    """
    Updates mapping datasets required for gene mapping, Downloads NCBI and UniProt mapping files.

    """
    
    location = os.path.dirname(os.path.realpath(__file__))
    out = os.path.join(location, 'datasets')
    ncbidb='https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
    tab_uni='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz'	 
    print('Downloading NCBI mapping file\n')
    wget.download(ncbidb,out)
    print('\nDownloading UniProt mapping file ')
    wget.download(tab_uni,out)
    print('\n!!DATABASES UPDATED!!')



def gene_mapping_many(query_list,source,target):
    
    """
    Maps gene identifiers from one namespace to another.

    Args:
        query_list (list): List of gene identifiers to be mapped.
        source (str): Source namespace of the gene identifiers.
        target (str): Target namespace to which the gene identifiers will be mapped.

    Returns:
        list: List of gene identifiers mapped to the target namespace.

    Raises:
        AssertionError: If input query_list contains non-integer values when source or target is 'entrez'.

    Note:
        The function relies on external data sources 'ncbidb' and 'tab_uni'.
    """
    
    

    if source=='ensembl' and target=='symbol':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((x.split('|')[0], y) for x, y in list(zip(ense,ncbidb['Symbol'])) if 'ENS' in x)
        return list(map(dictio.get,query_list))

    elif source=='ensembl' and target=='entrez':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((x.split('|')[0], y) for x, y in list(zip(ense,ncbidb['GeneID'])) if 'ENS' in x)
        return list(map(dictio.get,query_list))

    elif source=='ensembl' and target=='uniprot':
        dictio=dict(zip(list(reversed(tab_uni[tab_uni.mapper=='Ensembl'].id.tolist())),
                list(reversed(tab_uni[tab_uni.mapper=='Ensembl'].uniprot.tolist()))))
        return list(map(dictio.get,query_list))

    elif source=='entrez' and target=='ensembl':
        assert all(isinstance(x, int) for x in query_list), "All entrez_ids in query list should be integers"
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((y, x.split('|')[0]) for x, y in list(zip(ense,ncbidb['GeneID'])) if 'ENS' in x)
        return list(map(dictio.get,query_list))

    elif source == 'entrez' and target == 'symbol':
        assert all(isinstance(x, int) for x in query_list), "All entrez_ids in query list should be integers"
        dictio=dict((x, y) for x, y in list(zip(ncbidb['GeneID'],ncbidb['Symbol'])))
        return list(map(dictio.get,query_list))

    elif source=='entrez' and target=='uniprot':
        assert all(isinstance(x, int) for x in query_list), "All entrez_ids in query list should be integers"
        dictio=dict(zip([int(g) for g in list(reversed(tab_uni[tab_uni.mapper=='GeneID'].id.tolist()))],
                list(reversed(tab_uni[tab_uni.mapper=='GeneID'].uniprot.tolist()))))
        return list(map(dictio.get,query_list))

    elif source=='symbol' and target=='ensembl':
        ense=ncbidb['dbXrefs'].apply(lambda x : re.sub(r'.+?(?=ENS)', '', x))
        dictio=dict((y, x.split('|')[0]) for x, y in list(zip(ense,ncbidb['Symbol'])) if 'ENS' in x)
        return list(map(dictio.get,query_list))

    elif source == 'symbol' and target == 'entrez':
        dictio=dict((y, x) for x, y in list(zip(ncbidb['GeneID'],ncbidb['Symbol'])))
        return list(map(dictio.get,query_list))

    elif source=='symbol' and target=='uniprot':
        dictio=dict(zip(list(reversed(tab_uni[tab_uni.mapper=='Gene_Name'].id.tolist())),
                list(reversed(tab_uni[tab_uni.mapper=='Gene_Name'].uniprot.tolist()))))
        return list(map(dictio.get,query_list))

    elif source=='uniprot' and target=='ensembl':
        dictio=dict(zip(list(reversed(tab_uni[tab_uni.mapper=='Ensembl'].uniprot.tolist())),
                list(reversed(tab_uni[tab_uni.mapper=='Ensembl'].id.tolist()))))
        return list(map(dictio.get,query_list))

    elif source=='uniprot' and target=='entrez':
        dictio=dict(zip(list(reversed(tab_uni[tab_uni.mapper=='GeneID'].uniprot.tolist())),
                [int(g) for g in list(reversed(tab_uni[tab_uni.mapper=='GeneID'].id.tolist()))]))
        return list(map(dictio.get,query_list))


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


def VEP(idlist,input_type = 'rsid', chunked=False,chunksize=200,verbose=False,all_data=False,plot = False,save_plot=''):
    
    """
    Perform Variant Effect Prediction (VEP) using the Ensembl Variant Effect Predictor for a list of variant IDs.

    Parameters:
    - idlist (list or str): A list of variant IDs or a single variant ID as a string.
    - input_type (str, optional): Type of input IDs. Valid options are 'rsid' (default) or 'hgvs'.
    - chunked (bool, optional): Whether to chunk the requests into smaller portions if the number of IDs exceeds 200. Default is False.
    - chunksize (int, optional): The size of each chunk when chunking requests. Default is 200.
    - verbose (bool, optional): Whether to print verbose information during processing. Default is False.
    - all_data (bool, optional): Whether to return all data retrieved from VEP. Default is False.
    - plot (bool, optional): Whether to generate a pie chart visualizing the consequences of the variants. Default is False.
    - save_plot (str, optional): Filepath to save the generated plot. If not provided, the plot will be displayed instead.

    Returns:
    - list of DataFrames or dict: A list of DataFrames containing the VEP results (if `all_data` is False), or a dictionary containing all retrieved data (if `all_data` is True).

    Note:
    This function requires internet connectivity to access the Ensembl REST API.
    """
    
    if input_type == 'hgvs':
        http="https://rest.ensembl.org/vep/human/hgvs"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        chunked_idlist=[]
        if chunked | len(idlist) > 200:
            print('total number of chunks: %s' %(int(len(idlist)/chunksize)+1))
            for i in range(0,len(idlist),chunksize):
                chunked_idlist.append(idlist[i:i+chunksize])
            results=[]
            for i,chunk in enumerate(chunked_idlist):
                if response.ok:
                    response = requests.post(http, headers=headers, 
                                            data="{" + '"hgvs_notations" : {}'.format(str(chunk).replace("'",'"'))+"}")
                    results.append(response.json())
                    print('chunk %s processed' % (i))
                else:
                    print(f'ERROR: Bad Request:\n{response.text}')
        
            req_results=sum(results,[])
        else:
            results = requests.post(http, headers=headers, 
                                        data="{" + '"hgvs_notations" : {}'.format(str(idlist).replace("'",'"'))+"}")
            req_results=results.json()
    else:
        http="https://rest.ensembl.org/vep/human/id/"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        chunked_idlist=[]
        if chunked | len(idlist) > 200:
            print('total number of chunks: %s' %(int(len(idlist)/chunksize)+1))
            for i in range(0,len(idlist),chunksize):
                chunked_idlist.append(idlist[i:i+chunksize])
            results=[]
            for i,chunk in enumerate(chunked_idlist):
                response = requests.post(http, headers=headers, 
                                        data="{" + '"ids" : {}'.format(str(chunk).replace("'",'"'))+"}")
                if response.ok:
                    results.append(response.json())
                    print('chunk %s processed' % (i))
                else:
                    print(f'ERROR: Bad Request:\n{response.text}')
        
            req_results=sum(results,[])
        else:
            results = requests.post(http, headers=headers, 
                                        data="{" + '"ids" : {}'.format(str(idlist).replace("'",'"'))+"}")
            req_results=results.json()
        
    
    final_transcript_consequences=[]
    final_intergenic_consequences=[]
    final_regulatory_feature_consequences=[]
    final_motif_consequences=[]

    for dict_of_resu in req_results:
        variant_id=dict_of_resu['input']
        #Check Transcript Consequences
        try:
            transcript_consequences=dict_of_resu['transcript_consequences']
            for tc in transcript_consequences :
                consequence_terms=tc['consequence_terms']
                gene_id=tc['gene_id']
                biotype=tc['biotype']
                final_transcript_consequences.append((variant_id,str(consequence_terms),biotype,gene_id))

            n_of_tc=len(transcript_consequences)
        except Exception as error:
            n_of_tc=0



        #Check Intergenic Consequences
        try:
            intergenic_consequences=dict_of_resu['intergenic_consequences']
            for ic in intergenic_consequences :
                consequence_terms=ic['consequence_terms']
                impact=ic['impact']
                final_intergenic_consequences.append((variant_id,str(consequence_terms),impact))
                n_of_ic=len(intergenic_consequences)


        except Exception as error:
            n_of_ic=0




        #Check regulatory_feature_consequences
        try:
            regulatory_feature_consequences=dict_of_resu['regulatory_feature_consequences']

            for rfc in regulatory_feature_consequences :
                consequence_terms=rfc['consequence_terms']
                impact=rfc['impact']
                biotype=rfc['biotype']
                final_regulatory_feature_consequences.append((variant_id,str(consequence_terms),biotype,impact))
                n_of_rfc=len(regulatory_feature_consequences)


        except Exception as error:
            n_of_rfc=0
        
        #Check motif feature consequences
        try:
            motif_consequences=dict_of_resu['motif_feature_consequences']

            for mfc in motif_consequences :
                consequence_terms=mfc['consequence_terms']
                tf=mfc['transcription_factors']
                impact=mfc['motif_score_change']
                final_motif_consequences.append((variant_id,str(consequence_terms),str(tf),impact))
                n_of_mfc=len(motif_consequences)


        except Exception as error:
            n_of_mfc=0



        if verbose:
            print(f"{variant_id} has {n_of_tc} transcript consquence, {n_of_ic} intergenic consquence, {n_of_rfc} regulatory feature consequence, {n_of_mfc} motif feature consequences")

    if all_data:
        return req_results

    final_transcript_consequences_df = pd.DataFrame(list(set(final_transcript_consequences)),columns=['variantid','effects','biotype','ENSID'])
    final_regulatory_feature_consequences_df = pd.DataFrame(list(set(final_regulatory_feature_consequences)),columns=['variantid','effects','biotype','impact'])
    final_intergenic_consequences_df = pd.DataFrame(list(set(final_intergenic_consequences)),columns=['variantid','effects','impact'])
    final_motif_consequences_df = pd.DataFrame(list(set(final_motif_consequences)),columns=['variantid','effects','TF','impact_score'])
    lodfs = [final_intergenic_consequences_df,final_transcript_consequences_df,final_regulatory_feature_consequences_df,final_motif_consequences_df]
    if plot:
        rsid_consequences=sum([ast.literal_eval(tup[1]) for tup in final_transcript_consequences]+[ast.literal_eval(tup[1]) for tup in final_regulatory_feature_consequences]+[ast.literal_eval(tup[1]) for tup in final_intergenic_consequences]+[ast.literal_eval(tup[1]) for tup in final_motif_consequences],[])
        
        data=dict(Counter(rsid_consequences))
        counts=list(data.values())
        percentages = [(val/sum(counts))*100 for val in counts]
        names=data.keys()


        legend_labels=['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(names, percentages)]

        colors=plt.cm.tab20c.colors

        fig,ax=plt.subplots(figsize=(20,10))
        patches,text=ax.pie(data.values(),
            colors = colors,
            textprops={'fontsize': 12})


        ax.legend(patches, legend_labels, bbox_to_anchor=(1,1.050), loc="upper left",fontsize=12)
        
        fig.suptitle("VEP RESULTS",fontsize=15)
        plt.tight_layout()
        if save_plot:
            plt.savefig(save_plot,dpi=350)
        else:
            plt.show()
    return lodfs

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


def get_variants_in_LD(variant, r2, pop='EUR'):
    """
    Retrieve variants in linkage disequilibrium (LD) with the specified variant.

    Parameters:
    - variant (str): The variant ID for which LD variants are requested.
    - r2 (float): The minimum threshold for the linkage disequilibrium coefficient (r2).
    - pop (str, optional): The population for which LD information is requested. Default is 'EUR' (European population).

    Returns:
    - list: A list of variant IDs that are in linkage disequilibrium with the specified variant, based on the given r2 threshold.

    Note:
    This function queries the Ensembl REST API to retrieve LD information for the specified variant.
    It returns a list of variant IDs that are in LD with the specified variant, with an r2 coefficient greater than or equal to the specified threshold.
    If no variants are found or an error occurs during retrieval, an empty list is returned.
    """
    
    http = "https://rest.ensembl.org/ld/human/%s/1000GENOMES:phase_3:%s?r2=%s" % (variant, pop, r2)
    
    try:
        resp = requests.get(http, headers={ "Content-Type": "application/json"})
        if resp.ok:
            variants = resp.json()
            return [x['variation2'] for x in variants if float(x['r2']) >= r2]
        else:
            print(f'ERROR: Bad Request:\n{resp.text}')
    except Exception as e:
        print("Error:", e)
        return []


def get_ld_matrix(list_of_snps, token, pop='EUR', metric='r2'):
    """
    Get the linkage disequilibrium (LD) matrix for a list of SNPs.

    Parameters:
    - list_of_snps (list): A list of SNP IDs.
    - token (str): The authentication token for accessing the LDlink API.
    - pop (str, optional): The population for which LD information is requested. Default is 'EUR' (European population).
    - metric (str, optional): The LD metric to be used. Default is 'r2'.

    Returns:
    - pandas DataFrame: A DataFrame representing the LD matrix. Rows and columns are SNP IDs, and cell values represent the LD metric between corresponding SNPs.

    Note:
    This function queries the LDlink API to retrieve the LD matrix for the specified list of SNPs, for this reason it needs the user to be registered on LDlink in order to obtain a token.
    It returns a pandas DataFrame containing the LD matrix, where rows and columns are SNP IDs, and cell values represent the LD metric.
    If any LD metric values are missing (NA), they are replaced with None, and the DataFrame is converted to numeric data types.
    """
    snp_string = '\n'.join(list_of_snps)

    headers = {
        'Content-Type': 'application/json',
    }

    params = (
        ('token', token),
    )

    json_data = {
        'snps': snp_string,
        'pop': pop,
        'r2_d': metric,
        'genome_build': 'grch38',
    }

    response = requests.post('https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix', headers=headers, params=params, json=json_data, verify=False)
    if response.ok:
        response = response.text
        data_frame = pd.DataFrame([x.split('\t') for x in response.split('\n')])
        new_header = data_frame.iloc[0]
        data_frame = data_frame[1:]  # take the data less the header row
        data_frame.columns = new_header  # set the header row as the df header

        new_rows = data_frame[data_frame.columns[0]]
        data_frame = data_frame[data_frame.columns[1:]].set_index(new_rows)

        data_frame.replace('NA', None, inplace=True)
        data_frame = data_frame.astype(None)

        return data_frame.fillna(0).iloc[:-1]
    else:
        print(f'ERROR: Bad Request:\n{response.text}')


def pairwise_ld(ch, start, end, pop='EUR'):
    """
    Retrieve pairwise linkage disequilibrium (LD) information for variants within a genomic region.

    Parameters:
    - ch (str): Chromosome name or identifier.
    - start (int): Start position of the genomic region.
    - end (int): End position of the genomic region.
    - pop (str, optional): The population for which LD information is requested. Default is 'EUR' (European population).

    Returns:
    - pandas DataFrame: A DataFrame containing pairwise LD information between variants within the specified genomic region.
      Columns include 'v1' and 'v2' for variant IDs and 'r2' for the LD coefficient.

    Note:
    This function queries the Ensembl REST API to retrieve pairwise LD information for variants within the specified genomic region.
    It returns a DataFrame containing pairwise LD information, where 'v1' and 'v2' represent variant IDs and 'r2' represents the LD coefficient.
    """

    http = f"https://rest.ensembl.org/ld/human/region/{ch}:{str(start)}..{str(end)}/1000GENOMES:phase_3:{pop}"
    response = requests.get(http, headers={"Content-Type": "application/json"})
    if response.ok:
        response = response.json()
        ld_mat = []

        for i, element in enumerate(response):
            try:
                v1 = element['variation1']
                v2 = element['variation2']
                r2 = element['r2']
                ld_mat.append((v1, v2, r2))
            except Exception as r:
                print(r, f'- Error for variant "{element}"')
                ld_mat_df = pd.DataFrame(ld_mat, columns=['v1', 'v2', 'r2'])

        return ld_mat_df.sort_values(by='v1')
    else:
        print(f'ERROR: Bad Request:\n{response.text}')
        


def get_summary_statistic(study):
    """
    Retrieve summary statistics for a given GWAS study.

    Parameters:
    - study (str): The study identifier for which summary statistics are requested.

    Returns:
    - dict: A dictionary containing summary statistics for the specified GWAS study.

    Note:
    This function queries the EBI GWAS Catalog Summary Statistics API to retrieve summary statistics for the specified study.
    It returns a dictionary containing summary statistics information, including associations, for the specified study.
    """

    http = 'https://www.ebi.ac.uk/gwas/summary-statistics/api/studies/%s/associations' % (study)
    response = requests.get(http)
    if response.ok:
        response = response.json()
        return response()
    else:
        print(f'ERROR: Bad Request:\n{response.text}')


def get_summary_statistic_list():
    """
    Retrieve a list of summary statistics for all available GWAS studies.

    Returns:
    - list: A list of dictionaries containing summary statistics for all available GWAS studies.

    Note:
    This function queries the EBI GWAS Catalog Summary Statistics API to retrieve summary statistics for all available GWAS studies.
    It returns a list of dictionaries, where each dictionary contains summary statistics information for a specific GWAS study.
    """

    http = 'https://www.ebi.ac.uk/gwas/summary-statistics/api/associations'
    response = requests.get(http)
    if response.ok:
        response = response.json()
        return response()
    else:
        print(f'ERROR: Bad Request:\n{response.text}')


def get_phenotypes(chromosome, start, stop, feature_type='Genes', only_phenotypes=1):
    """
    Retrieve phenotype annotations annotations of a given genomic region.

    Parameters:
    - chromosome (str): Chromosome name or identifier.
    - start (int): Start position of the genomic region.
    - stop (int): End position of the genomic region.
    - feature_type (str, optional): Type of genomic feature to include in annotations. Default is 'Genes'.
    - only_phenotypes (int, optional): Whether to include only phenotypes. Default is 1 (True).

    Returns:
    - dict: A dictionary containing annotations of the specified genomic region.

    Note:
    This function queries the Ensembl REST API to retrieve annotations for the specified genomic region.
    It returns a dictionary containing annotations, including phenotypes, for the specified genomic region.
    """

    http = "https://rest.ensembl.org/phenotype/region/homo_sapiens/%s:%s-%s?only_phenotypes=%s;feature_type=%s" % (
    chromosome, start, stop, only_phenotypes, feature_type)
    annot = requests.get(http, headers={"Content-Type": "application/json"})
    return annot.json()


def get_ov_region(snp=None, chr=None, location=None, window=500, features=list, mode='region'):
    """
    Retrieve overlap regions for a specified SNP or genomic location.

    Parameters:
    - snp (str, optional): The variant ID (SNP) for which overlap regions are requested. Default is None.
    - chr (str or int, optional): Chromosome name or identifier. Required if mode is 'SNP'.
    - location (int, optional): Genomic location (start position) for which overlap regions are requested. Required if mode is 'region'.
    - window (int, optional): Window size to extend the genomic region around the SNP or location. Default is 500.
    - features (list, optional): A list of genomic feature types to include in the overlap regions. Default is an empty list.
    - mode (str, optional): Mode of operation, either 'SNP' or 'region'. Default is 'region'.

    Returns:
    - list of pandas DataFrames: A list of pandas DataFrames containing overlap regions for each specified genomic feature type.

    Note:
    This function queries the Ensembl REST API to retrieve overlap regions for the specified SNP or genomic location.
    It returns a list of pandas DataFrames, where each DataFrame contains overlap regions for a specific genomic feature type.
    """

    str_features = ';'.join(['feature=' + x for x in features])

    if mode == 'SNP':
        pos = get_variants_position(snp)
        chr = int(pos[0][1])
        genomic_location = pos[0][2]
        start = genomic_location - window // 2
        stop = genomic_location + window // 2
    else:

        start = location - window // 2
        stop = location + window // 2
    http = "https://rest.ensembl.org/overlap/region/human/%s:%s-%s?%s" % (chr, start, stop, str_features)
    risposta = requests.get(http, headers={"Content-Type": "application/json"})
    if risposta.ok:
        risposta = risposta.json()
        lodfs = []
        list_of_features_retrieved = list(set([e['feature_type'] for e in risposta]))
        for f in features:
            if f not in list_of_features_retrieved:
                print(f"Couldn't retrieve data for feature: {f}")
        for x in list_of_features_retrieved:
            tmp_list = [dict(sorted(e.items())) for e in risposta if e['feature_type'] == x]
            # find the most number of keys dictionary
            length_list = [len(e.items()) for e in tmp_list]
            max_l = max(length_list)
            index_max = length_list.index(max_l)
            cols = tmp_list[index_max].keys()
            e = pd.DataFrame(columns=cols)
            for i, f in enumerate(tmp_list):
                serie = pd.Series(data=f, index=list(f.keys()))
                e.loc[i] = serie

            lodfs.append(e)
        return lodfs
    else:
        print(f'ERROR: Bad Request:\n{response.text}')

        
        


def closest_genes(position_id, chromosome, position, window_size, type_of_gene=False, mode='start'):
    """
    Find the closest upstream and downstream genes to a given genomic position.

    Parameters:
    - position_id (str): Identifier for the genomic position.
    - chromosome (str or int): Chromosome name or identifier.
    - position (int): Genomic position for which closest genes are to be found.
    - window_size (int): Window size to extend the genomic region around the position.
    - type_of_gene (str, optional): Biotype of genes to consider. Default is False (include all genes).
    - mode (str, optional): Mode of calculating distances, either 'start' or 'absolute'. Default is 'start'.
        start:  In this case the distance will be calculated from the beginning of both upstream and downstream genes 
        absolute: Upstream distance is calcuated as distance (base pairs) to the end and downstream as distance to the start of the gene

    Returns:
    - tuple: A tuple containing the position ID, the closest upstream gene ID, and the closest downstream gene ID.

    Note:
    This function utilizes the `get_ov_region` function to retrieve genes within the specified window around the genomic position.
    It calculates distances from the specified position to the start (or end) of genes and finds the closest upstream and downstream genes.
    """
    elements = get_ov_region(chr=chromosome, location=position, features=['gene'], window=window_size)
    all_genes = []
    for i, r in elements[0].iterrows():
        if type_of_gene:
            if r['biotype'] == type_of_gene:
                all_genes.append((r['gene_id'], r['start'] - position, r['end'] - position))
        else:
            all_genes.append((r['gene_id'], r['start'] - position, r['end'] - position))

    positive_r = []
    negative_r = []
    for r in all_genes:
        # Look for the closest upstream gene 
        if r[1] > 0:
            positive_r.append(r)
        elif r[2] < 0:
            negative_r.append(r)

    if mode == 'start':
        if len(negative_r) > 0:
            closest_downstream_gene = max(negative_r, key=lambda x: x[1])[0]
        else:
            print(f'no downstream genes in the selected window for variant {position_id}')
            closest_downstream_gene = ''

        if len(positive_r) > 0:
            closest_upstream_gene = min(positive_r, key=lambda x: x[1])[0]
        else:
            print(f'no upstream genes in the selected window for variant {position_id}')
            closest_upstream_gene = ''

    elif mode == 'absolute':
        """ In this case the distance will be calculated from the beginning of the downstream and the end of the upstream """
        if len(negative_r) > 0:
            closest_downstream_gene = max(negative_r, key=lambda x: x[2])[0]
        else:
            print(f'no downstream genes in the selected window for variant {position_id}')
            closest_downstream_gene = ''

        if len(positive_r) > 0:
            closest_upstream_gene = min(positive_r, key=lambda x: x[2])[0]
        else:
            print(f'no upstream genes in the selected window for variant {position_id}')
            closest_upstream_gene = ''

    return (position_id, closest_upstream_gene, closest_downstream_gene)


def get_sequence(chromosome, start, stop):
    """
    Retrieve the genomic sequence for a given genomic region.

    Parameters:
    - chromosome (str or int): Chromosome name or identifier.
    - start (int): Start position of the genomic region.
    - stop (int): End position of the genomic region.

    Returns:
    - str: The genomic sequence for the specified genomic region.

    Note:
    This function queries the Ensembl REST API to retrieve the genomic sequence for the specified genomic region.
    It returns the genomic sequence as a string.
    """
    http = "https://rest.ensembl.org/sequence/region/human/%s:%s..%s?" % (chromosome, start, stop)
    risposta = requests.get(http, headers={"Content-Type": "application/json"}).json()
    return risposta['seq']


def grch_liftover(chromosome, start, end, source, target):
    """
    Lift over genomic coordinates from one assembly version to another.

    Parameters:
    - chromosome (str or int): Chromosome name or identifier.
    - start (int): Start position of the genomic region.
    - end (int): End position of the genomic region.
    - source (str): Source genome assembly version.
    - target (str): Target genome assembly version.

    Returns:
    - tuple or None: A tuple containing the chromosome, start, and end positions in the target assembly if successful, 
      otherwise returns None.

    Note:
    This function queries the Ensembl REST API to perform coordinate liftover from the source genome assembly to the target assembly.
    It returns a tuple containing the chromosome, start, and end positions in the target assembly if liftover is successful,
    otherwise returns None.
    """
    url = "https://rest.ensembl.org/map/human/%s/%s:%i..%i/%s?" % (source, chromosome, start, end, target)
    r = requests.get(url, headers={"Content-Type": "application/json"}).json()
    try:
        return (chromosome, r['mappings'][0]['mapped']['start'], r['mappings'][0]['mapped']['end'])
    except:
        return None


def get_eqtl_variant(rsid):
    """
    Retrieve eQTL (expression quantitative trait loci) associations for a given variant.

    Parameters:
    - rsid (str): The variant ID (e.g., rsID) for which eQTL associations are requested.

    Returns:
    - dict: A dictionary containing eQTL associations for the specified variant.

    Note:
    This function queries the EBI eQTL Catalog API to retrieve eQTL associations for the specified variant.
    It returns a dictionary containing eQTL association information, including target genes and other relevant data, for the specified variant.
    """

    url = 'http://www.ebi.ac.uk/eqtl/api/associations/%s' % (rsid)
    risp = requests.get(url)
    if risp.ok:
        risp = risp.json()
        return risp
    else:
        print(f'ERROR: Bad Request:\n{risp.text}')



def get_eqtl_df(rsid, p_value=0.005, increase_index=False):
    """
    Retrieve eQTL (expression quantitative trait loci) associations as a DataFrame for a given variant.

    Parameters:
    - rsid (str): The variant ID (e.g., rsID) for which eQTL associations are requested.
    - p_value (float, optional): Maximum p-value threshold for filtering eQTL associations. Default is 0.005.
    - increase_index (bool, optional): Whether to increase the DataFrame index by 1. Default is False.

    Returns:
    - pandas DataFrame or None: A DataFrame containing eQTL associations for the specified variant if successful, 
      otherwise returns None.

    Note:
    This function queries the EBI eQTL Catalog API to retrieve eQTL associations for the specified variant.
    It returns a pandas DataFrame containing eQTL association information, including variant ID, p-value, beta value, 
    target gene ID, tissue, study ID, and tissue name. It optionally filters associations based on the provided p-value threshold.
    """

    location = os.path.dirname(os.path.realpath(__file__))
    out = os.path.join(location, 'data')
    with open(out + '/uberon_dict.pickle', 'rb') as f:
        ubdict = pickle.load(f)
    url = 'http://www.ebi.ac.uk/eqtl/api/associations/%s?size=1000' % (rsid)
    response = requests.get(url)
    if response.ok:
        eqtls = response.json()
        try:
            eqtl_df = pd.DataFrame(columns=['variantid', 'p_value', 'log_pval', 'beta', 'alt', 'gene_id', 'tissue', 'study_id'])
            for ass in eqtls['_embedded']['associations'].keys():
                pval = eqtls['_embedded']['associations'][ass]['pvalue']
                nlog_pval = -np.log10(pval)
                beta = eqtls['_embedded']['associations'][ass]['beta']
                alt = eqtls['_embedded']['associations'][ass]['alt']
                geneid = eqtls['_embedded']['associations'][ass]['gene_id']
                tissue = eqtls['_embedded']['associations'][ass]['tissue']
                study = eqtls['_embedded']['associations'][ass]['study_id']
                eqtl_df.loc[ass] = [rsid, pval, nlog_pval, beta, alt, geneid, tissue, study]

            eqtl_df = eqtl_df.loc[eqtl_df.p_value <= p_value]
            eqtl_df.tissue = eqtl_df.tissue.apply(lambda x: x.replace('UBER_', 'UBERON_'))
            eqtl_df['tissue_name'] = list(map(ubdict.get, eqtl_df.tissue.tolist()))

            eqtl_df = eqtl_df.reset_index(drop=True)
            if increase_index:
                eqtl_df.index += 1
        except Exception as er:
            print(er, eqtls)
            return None
        return eqtl_df
    else:
        print(f'ERROR: Bad Request:\n{response.text}')


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

def get_genes(cr, location, window_size=10000, pop='EUR', features=['gene'], mode='all'):
    """
    Retrieve genes in a window centered around a genomic position and compute the distance between the position and all genes.

    Parameters:
    - cr (str): Chromosome identifier.
    - location (int): Genomic position around which the window is centered.
    - window_size (int, optional): Size of the window in base pairs. Default is 10,000.
    - pop (str, optional): Population for which to retrieve gene data. Default is 'EUR' (European).
    - features (list, optional): List of features to include in the retrieval. Default is ['gene'].
    - mode (str, optional): Retrieval mode. Options: 'all' (returns all genes and their distances), 
      'complete_data' (returns complete response from the API), 'closest_forward' (returns the closest gene 
      located forward of the position), 'closest_backward' (returns the closest gene located backward of 
      the position), 'closest_overall' (returns the closest gene regardless of direction). Default is 'all'.

    Returns:
    - dict or str: Depending on the mode parameter:
        - If mode is 'all', returns a dictionary where keys are gene names and values are their distances from the position.
        - If mode is 'complete_data', returns the complete response from the API.
        - If mode is 'closest_forward' or 'closest_backward', returns the name of the closest gene in the specified direction.
        - If mode is 'closest_overall', returns the name of the closest gene regardless of direction.
        - If no genes are found, returns a message indicating the absence of genes.
    """

    win_start = location - window_size // 2
    win_end = location + window_size // 2
    str_features = ';'.join(['feature=' + x for x in features])
    http = "https://rest.ensembl.org/overlap/region/human/%s:%s-%s?%s" % (cr, win_start, win_end, str_features)
    response = requests.get(http, headers={"Content-Type": "application/json"}).json()
    if response.ok:
        response = response.json()
        if mode == 'complete_data':
            return response
        elif mode == 'all':
            elements = {}
            for el in response:
                if el['biotype'] == 'protein_coding':
                    try:
                        elements[el['external_name']] = int(el['start'] - start)
                    except:
                        pass
            return elements
        elif mode == 'closest_forward':
            elements = {}
            for el in response:
                if el['biotype'] == 'protein_coding':
                    try:
                        elements[el['external_name']] = int(el['start'] - start)
                    except:
                        pass
            try:
                return min([(k, v) for (k, v) in elements.items() if v > 0], key=lambda x: x[1])
            except:
                return 'no_genes_forward'
        elif mode == 'closest_backward':
            elements = {}
            for el in response:
                if el['biotype'] == 'protein_coding':
                    try:
                        elements[el['external_name']] = int(el['start'] - start)
                    except:
                        pass
            try:
                return max([(k, v) for (k, v) in elements.items() if v < 0], key=lambda x: x[1])
            except:
                return 'no_genes_backward'
        elif mode == 'closest_overall':
            elements = {}
            for el in response:
                if el['biotype'] == 'protein_coding':
                    try:
                        elements[el['external_name']] = int(el['start'] - start)
                    except:
                        pass
            try:
                return min([(k, np.absolute(v)) for (k, v) in elements.items()], key=lambda x: x[1])
            except:
                return 'no_genes'
    else:
        print(f'ERROR: Bad Request:\n{response.text}')



def convert_variants(list_of_variants, source='variantid', target='rsid'):
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


def OT_L2G(list_of_variants, score=0.1, output='genes'):
    """
    Retrieve genes associated with variants from the Open Targets Genetics API.

    Parameters:
    - list_of_variants (list): List of variant identifiers.
    - score (float, optional): Threshold score for gene association. Genes with scores higher than this threshold will be included. Default is 0.1.
    - output (str, optional): Output format. 'genes' returns a list of genes, 'all' returns a DataFrame with variant, gene, and score information. Default is 'genes'.

    Returns:
    - list or DataFrame: Depending on the output parameter, either a list of genes or a DataFrame with variant, gene, and score information.

    Note:
    - The function retrieves genes associated with variants from the Open Targets Genetics API.
    - The 'score' parameter filters genes based on their association score with the variants.
    - 'output' can be set to 'genes' to return a list of genes or 'all' to return detailed information in a DataFrame.
    """

    query = """{
              genesForVariant(variantId:"%s"){
                overallScore
                gene{id}
              }
            }
            """
    OT_url = 'https://api.genetics.opentargets.org/graphql'

    results = {}
    for variant in list_of_variants:
        r = requests.post(OT_url, json={'query': query % (variant)})
        r = r.json()
        ResultsForVariant = []
        for data in r['data']['genesForVariant']:
            ResultsForVariant.append((data['gene']['id'], data['overallScore']))
        results[variant] = ResultsForVariant
    if output == 'all':
        raw_data = {key: sorted(value, key=lambda x: x[1], reverse=True) for key, value in results.items()}
        cols = ['id', 'gene', 'score']
        data = []
        for k, v in raw_data.items():
            for gene in v:
                data.append((k, gene[0], gene[1]))
        res = pd.DataFrame(data, columns=cols)
        return res
    else:
        return list(set(sum([[value[0] for value in values if value[1] > score] for (key, values) in results.items()], [])))

def make_eqtl_network(list_of_eqtls_df, tissue=None, gene=None, variant=None, gene_symbol=False, tissue_name=False):
    """
    Generate a network based on eQTL associations.

    Parameters:
    - list_of_eqtls_df (list): List of DataFrames containing eQTL information.
    - tissue (str, optional): Filter the network by tissue. Default is None.
    - gene (str, optional): Filter the network by gene. Default is None.
    - variant (list, optional): List of variants to include in the network. Default is None.
    - gene_symbol (bool, optional): Whether to use gene symbols instead of Ensembl IDs in the network. Default is False.
    - tissue_name (bool, optional): Whether to use tissue names instead of IDs in the network. Default is False.

    Returns:
    - networkx.Graph: A networkx graph representing the eQTL network.

    Note:
    - The function generates a (networkx) network based on eQTL associations from a list of DataFrames.
    - The network can be filtered by tissue, gene, or variant.
    - If 'gene_symbol' is True, gene symbols will be used instead of Ensembl IDs in the network.
    - If 'tissue_name' is True, tissue names will be used instead of IDs in the network.
    """

    Graph = nx.Graph()
    variant_nodes = []
    gene_nodes = []
    tissue_nodes = []
    edgelist_variants_genes = []
    edgelist_genes_tissues = []
    
    if variant:
        filtered_list = []
        for dif in list_of_eqtls_df:
            if dif.iloc[0].variantid in variant:
                filtered_list.append(dif)
        list_of_eqtls_df = filtered_list
            
    
    for eq_df in list_of_eqtls_df:            
        # Decide whether to use tissue name or id in the network (sometimes only the Id is present)
        # Mapping the ids to symbols       
        if gene_symbol:
            gene_identifier_in_the_network = 'gene_symbol'
            gene_ids = eq_df.gene_id.tolist()
            gene_symbols = gene_mapping_many(gene_ids, 'ensembl', 'symbol')
            eq_df[gene_identifier_in_the_network] = gene_symbols
            eq_df.dropna(subset='gene_symbol', inplace=True)
        else:
            gene_identifier_in_the_network = 'gene_id'

        if tissue_name:
            tissue_identifier_in_the_network = 'tissue_name'
            eq_df.dropna(subset='tissue_name', inplace=True)  # None Cannot be a node in the network
        else:
            tissue_identifier_in_the_network = 'tissue'
        
        # Apply the filters if needed
        if tissue:
            eq_df = eq_df[eq_df[tissue_identifier_in_the_network] == tissue]
        if gene:
            eq_df = eq_df[eq_df[gene_identifier_in_the_network] == gene]
            
        # Making the edges
        attributes = [{'beta': v, 'color': 'green' if v > 0 else 'red'} for v in eq_df.beta.tolist()]
        annotations = [{'annotation': 'is_expressed_in'} for _ in eq_df[tissue_identifier_in_the_network].tolist()]
        edgelist_variants_genes.append(list(zip(eq_df.variantid.tolist(), eq_df[gene_identifier_in_the_network].tolist(), attributes)))
        edgelist_genes_tissues.append(list(zip(eq_df[gene_identifier_in_the_network].tolist(), eq_df[tissue_identifier_in_the_network].tolist(), annotations)))
        
        # Making the nodes
        variant_nodes.append(list(set(eq_df.variantid.tolist())))
        gene_nodes.append(list(set(eq_df[gene_identifier_in_the_network].tolist())))
        tissue_nodes.append(list(set(eq_df[tissue_identifier_in_the_network].tolist())))
        
    # Adding nodes to the graph
    variant_nodes = list(sum(variant_nodes, []))
    gene_nodes = list(sum(gene_nodes, []))
    tissue_nodes = list(sum(tissue_nodes, []))
    
    Graph.add_nodes_from(variant_nodes, tipo='variant')
    Graph.add_nodes_from(gene_nodes, tipo='gene')
    Graph.add_nodes_from(tissue_nodes, tipo='tissue')
    
    # Merge the data in the edgelists
    edgelist_genes_tissues = list(sum(edgelist_genes_tissues, []))
    edgelist_variants_genes = list(sum(edgelist_variants_genes, []))
        
    Graph.add_edges_from(edgelist_variants_genes)
    Graph.add_edges_from(edgelist_genes_tissues)
    
    return Graph

# Function enrichment
def function_enrichment_visualization(list_of_genes):
    """
    Perform function enrichment analysis and visualize the network.

    Parameters:
    - list_of_genes (list): List of genes for analysis and visualization.

    Returns:
    - None

    Note:
    - This function performs function enrichment analysis on the provided list of genes and visualizes the network.
    - The interactome data is obtained through the HIPPIE dataset.
    - The function enrichment analysis is performed using GProfiler.
    - The network is interactive and visualized using Dash and Cytoscape.
    """
    
    def parse_interactome():
        location = os.path.dirname(os.path.realpath(__file__))
        interactome = pd.read_csv(os.path.join(location, 'data', 'hippie_interactome.sif'), header=None, sep=' ', usecols=[0, 2])
        interactome.columns = ['source', 'target']
        interactome.source = bp.gene_mapping_many(interactome.source.astype(int).tolist(), 'entrez', 'symbol')
        interactome.target = bp.gene_mapping_many(interactome.target.astype(int).tolist(), 'entrez', 'symbol')
        interactome.dropna(inplace=True)
        return interactome

    def enrichment_analysis(list_of_genes):
        gp = GProfiler(return_dataframe=True)
        enrichment = gp.profile(
            organism='hsapiens',
            query=list_of_genes,
            significance_threshold_method='bonferroni',
            no_iea=True,
            no_evidences=False
        )
        return enrichment

    def build_graph(subgraph, enrichment_dataframe):
        enrichment_dataframe.dropna(subset='p_value', inplace=True)
        ig_subgraph = ig.Graph.from_networkx(subgraph)
        pos_ = dict(zip([v['_nx_name'] for v in ig_subgraph.vs], [coord for coord in ig_subgraph.layout_auto()]))
        app = dash.Dash(__name__)
        cyto_node_data = list(zip(
            pos_.keys(),
            [coord[0] for coord in pos_.values()],
            [coord[1] for coord in pos_.values()]
        ))
        nodes = [
            {
                'data': {'id': str(_id), 'label': str(_id)},
                'position': {'x': 120 * x, 'y': 120 * y}
            }
            for _id, x, y in cyto_node_data
        ]

        edges = [
            {'data': {'source': source, 'target': target}}
            for source, target in subgraph.edges()
        ]

        elements = nodes + edges

        default_stylesheet = [
            {
                'selector': 'node',
                'style': {
                    'background-color': '#F5CEC5',
                    'border-color': 'black',
                    'border-width': '1',
                    'label': 'data(label)',
                    'width': '60',
                    'height': '60'
                }
            },
            {
                'selector': 'edge',
                'style': {
                    'line-color': 'red',
                    'width': '1'
                }
            }
        ]

        app.layout = html.Div([
            html.Header(html.H1(['Function enrichment analysis topology visualization'],
                                style={'textAlign': 'center', 'paddingBottom': '50px', 'border': '0px solid',
                                       'border-bottom': '1px solid black'})),

            html.Main([
                html.Div([
                    html.Label('P-value Slider'),
                    dcc.Slider(
                        id='pvalue_slider',
                        min=round(-np.log10(enrichment_dataframe['p_value'].max())),
                        max=round(-np.log10(enrichment_dataframe['p_value'].min())),
                        value=round(-np.log10(enrichment_dataframe['p_value'].max())),
                        marks=dict(
                            list(zip(
                                set(sorted([round(el) for el in -np.log10(enrichment_dataframe.p_value.tolist())])),
                                [{} for value in
                                 set([round(el) for el in -np.log10(enrichment_dataframe.p_value.tolist())])]))),
                        step=None),
                    html.Div(id='updatemode-output-container', style={'marginTop': 20}),
                    html.Br(style={'lineHeight': '4'}),
                    html.Label('Sources'),
                    dcc.RadioItems(
                        id='sources',
                        labelStyle={'display': 'flex'}
                    ),
                    html.Br(style={'lineHeight': '4'}),
                    html.Label('Function'),
                    dcc.Dropdown(id='function_dropdown'),
                    html.P(id='cytoscape-mouseoverNodeData-output')
                ],
                    style={'width': '20%', 'display': 'inline-block', 'float': 'left', 'paddingTop': '20px',
                           'paddingLeft': '50px'}
                ),
                html.Div([
                    cyto.Cytoscape(
                        id='cytoscape_network',
                        layout={'name': 'preset'},
                        style={'width': '100%', 'height': '800px'},
                        stylesheet=default_stylesheet,
                        elements=elements,
                        autoRefreshLayout=True
                    )
                ],
                    style={'width': '75%', 'float': 'right', 'position': 'relative', 'top': '20px'}
                )
            ])
        ])

        @app.callback(
            Output('updatemode-output-container', 'children'),
            Input('pvalue_slider', 'value')
        )
        def display_value(value):
            return '-log10(P_Value): %s' % value

        @app.callback(
            Output('sources', 'options'),
            Input('pvalue_slider', 'value')
        )
        def set_sources(selected_pvalue):
            return [{'label': i, 'value': i} for i in
                    set(enrichment_dataframe[-np.log10(enrichment_dataframe.p_value) >= selected_pvalue].source.tolist())]

        @app.callback(
            Output('function_dropdown', 'options'),
            Input('pvalue_slider', 'value'),
            Input('sources', 'value')
        )
        def set_functions(p_value, source):
            return [{'label': i, 'value': i} for i in set(
                enrichment_dataframe[
                    (-np.log10(enrichment_dataframe.p_value) >= p_value) & (
                            enrichment_dataframe.source == source)].name.tolist())]

        @app.callback(
            Output('cytoscape_network', 'stylesheet'),
            Input('sources', 'value'),
            Input('function_dropdown', 'value')
        )
        def update_network(fsource, ffunction):
            try:
                filt_enrich = enrichment_dataframe[
                    (enrichment_dataframe.name == ffunction) & (enrichment_dataframe.source == fsource)].intersections.values[0]
                new_stylesheet = [{
                    'selector': "[id='%s']" % ele,
                    'style': {
                        'background-color': 'black',
                        'line-color': 'black'
                    }
                } for ele in filt_enrich]
                return default_stylesheet + new_stylesheet
            except:
                return default_stylesheet

        app.run()

    edge_file = parse_interactome()
    hippie_net = nx.from_pandas_edgelist(edge_file)
    subg = hippie_net.subgraph(list_of_genes)
    enrichment = enrichment_analysis(list_of_genes)
    build_graph(subg, enrichment)


def plot_enrichment_analysis_network(list_of_genes, pvalue, colormap='cividis', edgecolor='red', mkcolor='grey', mkfsize=10000, layout='spring',
                                     mklinewidths=2, alpha=1, figsize=(40, 20), savefig=False, factor=1, k=10, cbarfontsize=10, labelling=True,
                                     legend=False, legend_fontsize=20, legend_titlefontsize=25, legend_location=(0.5, 0.0), legend_col=6,
                                     legend_labelspacing=1.5, legend_title='', legend_columnspacing=1.5, legend_handlelength=3,
                                     size_legend_nofelements=3, cbar_orientation='horizontal', cbar_loc=(1, 0.5),
                                     method_of_correction='bonferroni', no_evidences=False, no_iea=True, **kwargs):
    """
    Perform enrichment analysis and visualize the results as a network where nodes are the functions, size of the nodes is proportional to the number of genes involved in the specific function,
    and edges exist between nodes if the functions share genes. The thickness of the edge is proportional to the number of genes involved in the specific function.

    Parameters:
    - list_of_genes (list): List of genes for enrichment analysis and network visualization.
    - pvalue (float): P-value threshold for significance.
    - colormap (str): Colormap for visualization (default: 'cividis').
    - edgecolor (str): Color of the edges (default: 'red').
    - mkcolor (str): Color of the markers (default: 'grey').
    - mkfsize (int): Marker size (default: 10000).
    - layout (str): Layout algorithm for network visualization (default: 'spring').
    - mklinewidths (int): Line width of markers (default: 2).
    - alpha (float): Transparency of markers (default: 1).
    - figsize (tuple): Figure size (default: (40,20)).
    - savefig (bool): Whether to save the figure (default: False).
    - factor (int): Scaling factor for the layout (default: 1).
    - k (int): Optimal distance between nodes for spring layout (default: 10).
    - cbarfontsize (int): Font size of color bar (default: 10).
    - labelling (bool): Whether to label the nodes (default: True).
    - legend (bool): Whether to show legend (default: False).
    - legend_fontsize (int): Font size of legend (default: 20).
    - legend_titlefontsize (int): Font size of legend title (default: 25).
    - legend_location (tuple): Location of legend (default: (0.5,0.0)).
    - legend_col (int): Number of columns in legend (default: 6).
    - legend_labelspacing (float): Spacing between legend labels (default: 1.5).
    - legend_title (str): Title of legend (default: '').
    - legend_columnspacing (float): Spacing between legend columns (default: 1.5).
    - legend_handlelength (int): Length of legend handles (default: 3).
    - size_legend_nofelements (int): Number of elements in size legend (default: 3).
    - cbar_orientation (str): Orientation of color bar (default: 'horizontal').
    - cbar_loc (tuple): Location of color bar (default: (1, 0.5)).
    - method_of_correction (str): Method of multiple testing correction (default: 'bonferroni').
    - no_evidences (bool): Exclude electronic annotations (default: False).
    - no_iea (bool): Exclude Inferred from Electronic Annotation (IEA) evidence (default: True).
    - **kwargs: Additional keyword arguments.

    Returns:
    - None
    """
    gp = GProfiler(return_dataframe=True)
    df = gp.profile(organism='hsapiens',
                    query=list_of_genes,
                    significance_threshold_method=method_of_correction,
                    no_iea=no_iea,
                    no_evidences=no_evidences)

    if df.shape[0] == 0:
        print(f"Couldn't retrieve significantly enriched functions for the query list of genes:\n\n {list_of_genes}")
        return

    def labelling_without_overlapping(x, y, list_of_annotations, ax, verbose=False, **kwargs):
        class Point:
            def __init__(self, x, y):
                self.x = x
                self.y = y

        def doOverlap(ret1, ret2):
            l1 = Point(ret1[0, 0], ret1[1, 1])
            r1 = Point(ret1[1, 0], ret1[0, 1])
            l2 = Point(ret2[0, 0], ret2[1, 1])
            r2 = Point(ret2[1, 0], ret2[0, 1])

            if l1.x >= r2.x or l2.x >= r1.x:
                return False

            if (r1.y >= l2.y or r2.y >= l1.y):
                return False

            return True

        annotations_coord = []
        for i, dot in enumerate(y):
            x_coords = x[i]
            y_coords = y[i]
            annotation = ax.annotate(str(list_of_annotations[i]),
                                     xy=(x[i], y[i]),
                                     xytext=(x_coords, y_coords),
                                     **kwargs)

            ax.figure.canvas.draw()
            bbox = matplotlib.text.Text.get_window_extent(annotation)
            bbox_data = ax.transData.inverted().transform(bbox)
            factor = 0.2 * (bbox_data[0, 0] - bbox_data[1, 0])
            annotations_coord.append(bbox_data)

            theta = np.radians(np.linspace(1, 360 * 200, 500))
            r = np.linspace(0, max(max(zip(x, y))), len(theta))
            x_2 = r * np.cos(theta) + x_coords
            y_2 = r * np.sin(theta) + y_coords
            n = 0
            keep_cycling = True
            while keep_cycling:
                keep_cycling = False
                for ind, box in enumerate(annotations_coord[0:-1]):
                    if doOverlap(box, bbox_data):
                        annotation.set_x(x_2[n])
                        annotation.set_y(y_2[n])
                        n += 1
                        ax.figure.canvas.draw()
                        bbox = matplotlib.text.Text.get_window_extent(annotation)
                        bbox_data = ax.transData.inverted().transform(bbox)
                        annotations_coord.pop()
                        annotations_coord.append(bbox_data)
                        keep_cycling = True
                        break

    maxpv = max([-np.log10(p) for p in df.p_value.tolist()])
    for i, (s, v) in enumerate(zip(df.source.value_counts().index, df.source.value_counts())):
        data = df[(df.source == s) & (df.p_value < pvalue)].reset_index()
        if data.shape[0] == 0:
            continue
        else:
            nxen = nx.Graph()
            for i, r in data.iterrows():
                nxen.add_node(r['name'], size=r['intersection_size'], pvalue=-np.log10(r['p_value']))

            for i, r in data.iterrows():
                for index, row in data.iloc[i + 1:].reset_index().iterrows():
                    if len(set(r['intersections']).intersection(set(row['intersections']))) > 0:
                        nxen.add_edge(r['name'],
                                      row['name'],
                                      weight=len(set(r['intersections']).intersection(set(row['intersections']))))

            if layout == 'spring':
                pos_ = nx.spring_layout(nxen, k)
            elif layout == 'auto':
                ig_subgraph = ig.Graph.from_networkx(nxen)
                pos_ = dict(zip([v['_nx_name'] for v in ig_subgraph.vs], [coord for coord in ig_subgraph.layout_auto()]))

            connections = [edge[2]['weight'] for edge in nxen.edges(data=True)]
            if len(connections) != 0:
                norm_connections = [(x - min(connections)) / (max(connections) - min(connections)) for x in connections]
            else:
                connections = norm_connections

            markers = [node[1]['size'] for node in nxen.nodes(data=True)]
            if len(markers) != 0:
                norm_markers = [(x - min(markers)) / (max(markers) - min(markers)) for x in markers]
            else:
                markers = norm_markers

            norm_markers = np.clip(norm_markers, 0.3, 1)
            
            fig, ax = plt.subplots(figsize=figsize)
            xses, yses, lab, colors = [], [], [], []
            for node in nxen.nodes(data=True):
                xses.append(pos_[node[0]][0])
                yses.append(pos_[node[0]][1])
                lab.append(node[0])
                colors.append(node[1]['pvalue'])
            number_lab = [str(e[0]) for e in list(enumerate(lab))]
            nodez_for_legend = ax.scatter(xses, yses, s=markers)
            nodez = ax.scatter(xses, yses, s=[mkfsize * size for size in norm_markers],
                               c=colors, cmap=colormap, vmax=maxpv, alpha=alpha, edgecolors=mkcolor,
                               linewidths=mklinewidths, clip_on=False, zorder=1)

            if labelling:
                labelling_without_overlapping(xses, yses, number_lab, ax, **kwargs)

            for indx, edge in enumerate(nxen.edges(data=True)):
                if edge[2]['weight'] > 0:
                    path_1 = edge[0]
                    path_2 = edge[1]
                    x0, y0 = pos_[path_1]
                    x1, y1 = pos_[path_2]
                    edgez = ax.plot(np.linspace(x0, x1), np.linspace(y0, y1),
                                    color=edgecolor,
                                    linewidth=3 * norm_connections[indx] ** 4,
                                    zorder=0)

            cbar = plt.colorbar(nodez, ax=ax, orientation=cbar_orientation, panchor=cbar_loc)
            cbar.set_label(r'$-log_{10}(p-value)$', fontsize=cbarfontsize + 4)
            cbar.ax.tick_params(labelsize=cbarfontsize)

            if legend:
                class TextHandlerB(HandlerBase):
                    def create_artists(self, legend, text, xdescent, ydescent,
                                        width, height, fontsize, trans):
                        tx = Text(width / 2., height / 2, text, fontsize=fontsize,
                                  ha="center", va="center", fontweight="bold")
                        return [tx]

                Legend.update_default_handler_map({str: TextHandlerB()})

                first_legend = fig.legend(number_lab, lab, bbox_to_anchor=(1, 0.5), loc="lower left",fontsize=legend_fontsize,)
                plt.gca().add_artist(first_legend)

                handles, _ = nodez.legend_elements(prop="sizes", alpha=0.6, num=size_legend_nofelements)
                _, label_markers = nodez_for_legend.legend_elements(prop="sizes", alpha=0.6)

                legend = fig.legend(handles, label_markers, fontsize=legend_fontsize, loc="upper left",
                                    bbox_to_anchor=(1, 0.5), ncol=legend_col, labelspacing=legend_labelspacing,
                                    columnspacing=legend_columnspacing, handlelength=legend_handlelength, frameon=False)

                legend.set_title(legend_title, prop={'size': legend_titlefontsize})

            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.set_title(s, fontsize=35)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            axis = plt.gca()
            axis.set_xlim([factor * x for x in axis.get_xlim()])
            axis.set_ylim([factor * y for y in axis.get_ylim()])
            plt.tight_layout()
            if savefig:
                plt.savefig(str(s) + 'enrichment_analysis.jpeg', dpi=300, bbox_inches='tight')

            plt.show()
  
    


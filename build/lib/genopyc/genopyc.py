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
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input,Output
import dash_cytoscape as cyto



def get_associations(efotrait,verbose=False, studyid = False):
    loe=[]
    """Retrieve snps associated to an EFO trait"""
    df=pd.DataFrame(columns=['variantid','p-value','risk_allele','RAF','beta','CI','mapped_gene'])
    http= 'https://www.ebi.ac.uk/gwas/rest/api/efoTraits/%s/associations' %(efotrait)
    if verbose:
        print('querying associations... \n')
    resp=requests.get(http)
    if resp.ok:
        associ=resp.json()
        if verbose:
            print('building the dataframe...')
        for i,element in enumerate(associ['_embedded']['associations']):
            try:
                variantid = ''.join(element['loci'][0]['strongestRiskAlleles'][0]['riskAlleleName'].split('-')[0:1])
                df.at[i,'variantid']= variantid
                df.at[i,'risk_allele']=element['loci'][0]['strongestRiskAlleles'][0]['riskAlleleName'].split('-')[-1]
                df.at[i,'mapped_gene']= ' '.join([str(elem) for elem in [e['geneName'] for e in element['loci'][0]['authorReportedGenes']]])
                df.at[i,'p-value']=float(element['pvalueMantissa'])*10**int(element['pvalueExponent'])
                
                try: 
                    df.at[i,'RAF']=float(element['riskFrequency'])
                except:
                    df.at[i,'RAF']=np.nan
                    
                df.at[i,'beta']=[float(element['betaNum']) if type(element['betaNum'])==float else np.nan][0]
                df.at[i,'SE']=element['standardError']
                df.at[i,'CI']=element['range']
                try:
                    study_link = element['_links']['study']['href']

                    df.at[i,'study_url'] = study_link
                except Exception as e:
                    
                    df.at[i,'study_url'] = np.nan
                    print(e, f"Couldn't retreive studylink for variant {variantid}")
                if studyid:
                    try:
                        jsonres = requests.get(study_link).json()
                        studyid = jsonres['accessionId']
                        df.at[i,'studyid'] = studyid
                    
                    except Exception as e:
                    
                        df.at[i,'studyid'] = np.nan
                        print(e, f"Couldn't retreive studyID for variant {variantid}")
                            
                        
            except Exception as e:
                print(f'error {e} for element {element}')
                pass
        df.fillna(np.nan, method = None,inplace = True)
        df['p-value'] = df['p-value'].map("{:.1e}".format)
    else:
         print(f'ERROR: Bad Resquest: \n {resp}')
    return df



#retrieve the coordinates of many genes
def get_genes_position(idlist,chunked=False,chunksize=200):
    if not type(idlist) == list:
        idlist = [idlist]

    """ This function accept a list of ensembl Ids and return the coordinates in the GHRC38 if the list is longer than 200 it needs to be chunked because ensembl accepts maximum a request of 200 Ids per time"""

    http="https://rest.ensembl.org/lookup/id"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    if chunked | len(idlist) > 200 :
        chunked_idlist=[]
        print('total number of chunks: %s' %(int(len(idlist)/chunksize)+1))
        for i in range(0,len(idlist),chunksize):
            chunked_idlist.append(idlist[i:i+chunksize])
        results=[]
        for i,chunk in enumerate(chunked_idlist):
            response = requests.post(http, headers=headers, 
                                     data="{" + '"ids" : {}'.format(str(chunk).replace("'",'"'))+"}").json()

            ListOfTuples=[]
            for k,v in response.items():
                try:
                    ListOfTuples.append((k,int(v['seq_region_name']),v['start'],v['end']))

                except:

                    print(f"Couldn't retrieve position for gene {k}") 
                    continue

            results.append(ListOfTuples)
            print('chunk %s processed' % (i))
        return sum(results,[])
    else:
        response = requests.post(http, headers=headers, 
                                     data="{" + '"ids" : {}'.format(str(idlist).replace("'",'"'))+"}").json()

        ListOfTuples=[]

        for k,v in response.items():
            try:
                ListOfTuples.append((k,int(v['seq_region_name']),v['start'],v['end']))
            except:

                print(f"Couldn't retrieve position for gene {k}") 
                pass
        return ListOfTuples

def get_variants_position(idlist,chunked=False,chunksize=200):
    if not type(idlist) == list:
        idlist = [idlist]
        
    http="https://rest.ensembl.org/variation/homo_sapiens"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    if chunked | len(idlist) > 200:
        chunked_idlist=[]
        print('total number of chunks: %s' %(int(len(idlist)/chunksize)+1))
        for i in range(0,len(idlist),chunksize):
            chunked_idlist.append(idlist[i:i+chunksize])
        results=[]
        for i,chunk in enumerate(chunked_idlist):
            response = requests.post(http, headers=headers, 
                                     data="{" + '"ids" : {}'.format(str(chunk).replace("'",'"'))+"}").json()
            for key,value in response.items():
                try:
                    chr=value['mappings'][0]['location'].split(':')[0]
                    pos=value['mappings'][0]['start']
                    results.append((key,chr,pos))
                
                except:
                    print(f"Couldn't Retrieve Position for variant {key}")
                    pass
            print(f"chunk {i} processed")        
        return results        
            
    else:
        response = requests.post(http, headers=headers, 
                                     data="{" + '"ids" : {}'.format(str(idlist).replace("'",'"'))+"}").json()

        results=[]
        for key,value in response.items():
            try:
                chr=value['mappings'][0]['location'].split(':')[0]
                pos=value['mappings'][0]['start']
                results.append((key,chr,pos))
            
            except:
                print(f"Couldn't Retrieve Position for variant {key}")
                pass        
        return results

        

        

def VEP(idlist,input_type = 'rsid', chunked=False,chunksize=200,verbose=False,all_data=False,plot = False,save_plot=False):
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
                response = requests.post(http, headers=headers, 
                                        data="{" + '"hgvs_notations" : {}'.format(str(chunk).replace("'",'"'))+"}")
                results.append(response.json())
                print('chunk %s processed' % (i))
        
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
                results.append(response.json())
                print('chunk %s processed' % (i))
        
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
            plt.savefig('Vep_results.png',dpi=350)
        plt.show()
    return lodfs
    
## return single variant info ##
def get_variants_info(idlist, chunked=False,chunksize=200):
    if not type(idlist) == list:
        idlist = [idlist]
    http="https://rest.ensembl.org/variation/homo_sapiens"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    chunked_idlist=[]
    if chunked | len(idlist)>200:
        print('total number of chunks: %s' %(int(len(idlist)/chunksize)+1))
        for i in range(0,len(idlist),chunksize):
            chunked_idlist.append(idlist[i:i+chunksize])
        results={}
        for i,chunk in enumerate(chunked_idlist):
            response = requests.post(http, headers=headers, 
                                     data="{" + '"ids" : {}'.format(str(chunk).replace("'",'"'))+"}")
            results.update(response.json())
            print('chunk %s processed' % (i))
        return results

    else:
        response = requests.post(http, headers=headers, 
                                     data="{" + '"ids" : {}'.format(str(idlist).replace("'",'"'))+"}")
        return response.json()

        
def get_ancestral_allele(variants, mode = 'most_frequent'):
    variant_info = get_variants_info(variants)
    variant_mappings = dict(zip(variant_info.keys(),[info['mappings'] for info in variant_info.values()]))
    if mode == 'most_frequent':
        def most_frequent(List):
            return max(set(List), key = List.count)


        variant_ancestral_alleles={}
        for rsid,mappings in variant_mappings.items():
            if len(mappings)>1:
                list_of_alleles=[]
                for element in mappings:
                    ref_allele=element['allele_string'].split('/')[0]
                    list_of_alleles.append(ref_allele)        
                try:
                    variant_ancestral_alleles[rsid]=most_frequent(list_of_alleles)
                except Exception as error:
                    print(f'Error for rsid {rsid}: ' + repr(error))
            else:
                try:
                    ref_allele=mappings[0]['allele_string'].split('/')[0]
                    variant_ancestral_alleles[rsid]=mappings[0]['ancestral_allele']
                except Exception as error:
                    print(f'Error for rsid {rsid}: ' + repr(error) + ' MAPPINGS = ' + str(mappings))

        variant_ancestral_alleles={k:v for k,v in variant_ancestral_alleles.items() if v !=None}
        return variant_ancestral_alleles        
#Retrieves variant in LD with a given variant
def get_variants_in_LD(variant,r2,pop='EUR'):
    http= "https://rest.ensembl.org/ld/human/%s/1000GENOMES:phase_3:%s?r2=%s" %(variant,pop,r2)
    try:
        variants=requests.get(http,headers={ "Content-Type" : "application/json"}).json()

        return [x['variation2'] for x in variants if float(x['r2'])>=r2]
    except:
        pass

    
    
    
def GetLDMatrix(ListOfSnps,token,pop='EUR',metric='r2'):
    SnpString='\n'.join(ListOfSnps)

    headers = {
        'Content-Type': 'application/json',
    }

    params = (
        ('token', token),
    )

    json_data = {
        'snps': SnpString,
        'pop': pop,
        'r2_d': metric,
        'genome_build': 'grch38',
    }

    response = requests.post('https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix', headers=headers, params=params, json=json_data, verify=False).text
    dataf=pd.DataFrame([x.split('\t') for x in response.split('\n')])
    new_header=dataf.iloc[0]
    dataf= dataf[1:] #take the data less the header row
    dataf.columns = new_header #set the header row as the df header
    
    new_rows=dataf[dataf.columns[0]]
    dataf=dataf[dataf.columns[1:]].set_index(new_rows)
    
    dataf.replace('NA',None,inplace = True)
    dataf =dataf.astype(None)
    
    return dataf.fillna(0).iloc[:-1]
    
    
    
def PairwiseLD(ch,start,end,pop = 'EUR'):
    http= f"https://rest.ensembl.org/ld/human/region/{ch}:{str(start)}..{str(end)}/1000GENOMES:phase_3:{pop}"
    response = requests.get(http,headers={ "Content-Type" : "application/json"}).json()
    ldmat = []
    
    
    for i,element in enumerate(response):
        try:
            v1 = element['variation1']
            v2 = element['variation2']
            r2 = element['r2']
            ldmat.append((v1,v2,r2))
        except Exception as r:
            print(r,f'- Error for variant "{element}"')
        ldmatdf = pd.DataFrame(ldmat,columns=['v1','v2','r2'])
        
        


        
    
    
    return ldmatdf.sort_values(by = 'v1')
   
    

#Retrieve summary statistic of a given study    
def get_summary_statistic(study):
    http= 'https://www.ebi.ac.uk/gwas/summary-statistics/api/studies/%s/associations' %(study)
    ss=requests.get(http).json()
    return ss


#Retrieves the list of studies with summary statistics available for a given trait
def get_summary_statistic_list():
    http= 'https://www.ebi.ac.uk/gwas/summary-statistics/api/associations'
    ss=requests.get(http).json()
    return ss

#Retrieves annotations of a given genomic region
def get_phenotypes(chromosome,start,stop,feature_type='Genes',only_phenotypes=1):
    http="https://rest.ensembl.org/phenotype/region/homo_sapiens/%s:%s-%s?only_phenotypes=%s;feature_type=%s"%(chromosome,start,stop,only_phenotypes,feature_type)
    annot=requests.get(http, headers={ "Content-Type" : "application/json"})
    return annot.json()

#Retrieves overlapping elements of a given region 
#Retrieves overlapping elements of a given region
def get_ov_region(snp= None, chr=None,location=None, window=500, features=list, mode='region'):
    str_features=';'.join(['feature='+x for x in features])

    if mode == 'SNP':
        pos = get_variants_position(snp)
        chr = int(pos[0][1])
        genomic_location = pos[0][2]
        start=genomic_location-window//2
        stop=genomic_location+window//2
    else:

        start=location-window//2
        stop=location+window//2
    http="https://rest.ensembl.org/overlap/region/human/%s:%s-%s?%s"%(chr,start,stop,str_features)
    risposta=requests.get(http,headers={ "Content-Type" : "application/json"}).json()
    lodfs = []
    list_of_features_retrieved = list(set([e['feature_type'] for e in risposta]))
    for x in list_of_features_retrieved:
        tmp_list = [dict(sorted(e.items())) for e in risposta if e['feature_type'] == x]
        #find the most number of keys dictionary
        length_list = [len(e.items()) for e in tmp_list]
        max_l = max(length_list)
        index_max = length_list.index(max_l)
        cols = tmp_list[index_max].keys()
        e = pd.DataFrame(columns=cols)
        for i,f in enumerate(tmp_list):
            serie = pd.Series(data = f, index=list(f.keys()))
            e.loc[i] = serie

        lodfs.append(e)
    return lodfs

def ClosestGenes(positionid,chromosome,position,window_size,type_of_gene = False):
    """ Retrieve the closeset upstream and downstrem genes given a point genomic location """
    elements =  get_ov_region(chr=chromosome,location = position, features=['gene'],window=window_size)
    all_genes = []
    for i,r in elements[0].iterrows():
        if type_of_gene:
            if r['biotype'] == type_of_gene:
                all_genes.append((r['gene_id'],r['start'] - position,r['end'] - position))
        else:
            
            all_genes.append((r['gene_id'],r['start'] - position,r['end'] - position))

            
    
    positive_r = []
    negative_r = []
    for r in all_genes:
        #Look for the closest upstream gene 
        if r[1] > 0:
            positive_r.append(r)
        elif r[2] < 0:
            negative_r.append(r)
    if len(negative_r)>0:
        closest_downstream_gene = max(negative_r,key = lambda x: x[2])[0]
    else:
        print(f'no downstream genes in the selected window for variant {positionid}')
        closest_downstream_gene = ''

    if len(positive_r)>0:
        closest_upstream_gene = min(positive_r,key = lambda x: x[2])[0]
    else:
        print(f'no upstream genes in the selected window for variant {positionid}')
        closest_upstream_gene = ''

    return (positionid,closest_upstream_gene,closest_downstream_gene)



#Retrieves the nucleotide sequence of a given position 
def get_sequence(chromosome,start,stop):
    http="https://rest.ensembl.org/sequence/region/human/%s:%s..%s?" %(chromosome,start,stop)
    risposta=requests.get(http,headers={ "Content-Type" : "application/json"}).json()
    return risposta['seq']

## lift from grch37 to 38
def grch_liftover(chromosome,start,end,source,target):
    url="https://rest.ensembl.org/map/human/%s/%s:%i..%i/%s?"%(source,chromosome,start,end,target)
    r = requests.get(url, headers={ "Content-Type" : "application/json"}).json()
    try:
        return (chromosome,r['mappings'][0]['mapped']['start'],r['mappings'][0]['mapped']['end'])
    except:
        return None

## function to get variant list of eqtls
def get_eqtl_variant(rsid):
    url='http://www.ebi.ac.uk/eqtl/api/associations/%s'%(rsid)
    risp=requests.get(url).json()
    return risp


def get_eqtl_df(rsid,p_value=0.005,increase_index=False):
    location = os.path.dirname(os.path.realpath(__file__))
    out = os.path.join(location, 'data')
    with open(out+'/uberon_dict.pickle','rb') as f:
        ubdict = pickle.load(f)
    url='http://www.ebi.ac.uk/eqtl/api/associations/%s?size=1000'%(rsid)
    response=requests.get(url)
    eqtls=response.json()
    try:
        eqtl_df=pd.DataFrame(columns=['variantid','p_value','log_pval','beta','alt','gene_id','tissue','study_id'])
        for ass in eqtls['_embedded']['associations'].keys():
            pval=eqtls['_embedded']['associations'][ass]['pvalue']
            nlog_pval=-np.log10(pval)
            beta=eqtls['_embedded']['associations'][ass]['beta']
            alt=eqtls['_embedded']['associations'][ass]['alt']
            geneid=eqtls['_embedded']['associations'][ass]['gene_id']
            tissue=eqtls['_embedded']['associations'][ass]['tissue']
            study=eqtls['_embedded']['associations'][ass]['study_id']
            eqtl_df.loc[ass]=[rsid,pval,nlog_pval,beta,alt,geneid,tissue,study]
            
        eqtl_df=eqtl_df.loc[eqtl_df.p_value<=p_value]
        eqtl_df.tissue = eqtl_df.tissue.apply(lambda x: x.replace('UBER_','UBERON_'))
        eqtl_df['tissue_name'] = list(map(ubdict.get,eqtl_df.tissue.tolist()))
        
        eqtl_df=eqtl_df.reset_index(drop=True)
        if increase_index:
            eqtl_df.index+=1
    except Exception as er:
        print(er,eqtls)
        return None
    return eqtl_df






def get_genes(cr,location,window_size=10000,pop='EUR',features=['gene'],mode='all'):
    
    """
    Retrieve genes in a window centered in a genomic position and compute the distance between the position and all the genes
    """
    winstart=location-window_size//2
    winend=location+window_size//2
    str_features=';'.join(['feature='+x for x in features])
    http="https://rest.ensembl.org/overlap/region/human/%s:%s-%s?%s"%(cr,winstart,winend,str_features)
    risposta=requests.get(http,headers={ "Content-Type" : "application/json"}).json()
    if mode=='complete_data':
        return risposta
    elif mode=='all':
        elements={}
        for el in risposta:
            if el['biotype']=='protein_coding':
                try:
                    elements[el['external_name']]=int(el['start']-start)
                except:
                    pass
        return elements
    elif mode=='closest_forward':
        elements={}
        for el in risposta:
            if el['biotype']=='protein_coding':
                try:
                    elements[el['external_name']]=int(el['start']-start)
                except:
                    pass
        try:
            return min([(k,v) for (k,v) in elements.items() if v>0], key=lambda x:x[1])
        except:
            return 'no_genes_forward'
    elif mode=='closest_backward':
        elements={}
        for el in risposta:
            if el['biotype']=='protein_coding':
                try:
                    elements[el['external_name']]=int(el['start']-start)
                except:
                    pass
        try:
            return max([(k,v) for (k,v) in elements.items() if v<0], key=lambda x:x[1])
        except:
            return 'no_genes_backward'
    elif mode=='closest_overall':
        elements={}
        for el in risposta:
            if el['biotype']=='protein_coding':
                try:
                    elements[el['external_name']]=int(el['start']-start)
                except:
                    pass
        try:
            return min([(k,np.absolute(v)) for (k,v) in elements.items()], key=lambda x:x[1])
        except:
            return 'no_genes'


        
def convert_variants(ListOfVariants,source='variantid',target='rsid'):
    
    OT_url='https://api.genetics.opentargets.org/graphql'
    MappingDict={}
    if (source=='variantid') & (target=='rsid'):
        query="""{
               variantInfo(variantId:"%s"){
                   rsId
                   }
              }"""
        for variant in ListOfVariants:
            r = requests.post(OT_url, json={'query': query % (variant)})
            JsonResponse=r.json()
            MappingDict[variant]=JsonResponse['data']['variantInfo']['rsId']
        return list(map(MappingDict.get,ListOfVariants))

    elif (source=='rsid') & (target=='variantid'):
        query="""{
              search(queryString:"%s"){
                  variants{
                      id
                      }
                   }
              }
              """
        for variant in ListOfVariants:
            try:
                r = requests.post(OT_url, json={'query': query % (variant)})
                JsonResponse=r.json()
                MappingDict[variant]=JsonResponse['data']['search']['variants'][0]['id']
            except Exception as e:
                print(e, f"Couldn't Convert Variant {variant}")
        return list(map(MappingDict.get,ListOfVariants))


    
def OT_L2G(ListofVariants,score=0.1,output='genes'):
    query="""{
              genesForVariant(variantId:"%s"){
                overallScore
                gene{id}
              }
            }
            """
    OT_url='https://api.genetics.opentargets.org/graphql'
    
    results={}
    for variant in ListofVariants:
        r = requests.post(OT_url, json={'query': query % (variant)})
        r = r.json()
        ResultsForVariant=[]
        for data in r['data']['genesForVariant']:
            ResultsForVariant.append((data['gene']['id'],data['overallScore']))
        results[variant]=ResultsForVariant
    if output=='all':
        raw_data = {key:sorted(value, key=lambda x: x[1],reverse=True) for key,value in results.items()}
        cols = ['id','gene','score']
        data = []
        for k,v in raw_data.items():
            for gene in v:
                data.append((k,gene[0],gene[1]))
        res = pd.DataFrame(data,columns=cols)
        return res
    else:
        return list(set(sum([[value[0] for value in values if value[1]>score] for (key,values) in results.items()],[])))
    



#Function enrichment
def FuEnViz(ListOfGenes):
    def ParseInteractome():
        location = os.path.dirname(os.path.realpath(__file__))
        Interactome=pd.read_csv(os.path.join(location, 'data', 'hippie_interactome.sif'),header=None, sep=' ',usecols=[0,2])
        Interactome.columns=['source','target']
        # if ((args.sep=='space') & (args.id=='symbol')):
        Interactome.source=bp.gene_mapping_many(Interactome.source.astype(int).tolist(),'entrez','symbol')
        Interactome.target=bp.gene_mapping_many(Interactome.target.astype(int).tolist(),'entrez','symbol')
        Interactome.dropna(inplace=True)
        return Interactome



    def EnrichmentAnalisys(ListOfGenes):
        ## PERFORM ENRICHMENT ANALISYS##
        gp = GProfiler(return_dataframe=True)
        Enrichment=gp.profile(organism='hsapiens',
                    query=ListOfGenes,
                    significance_threshold_method='bonferroni',
                    no_iea=True,
                    no_evidences=False)
        return Enrichment

    def BuildGraph(SubGraph,EnrichmentDataframe):
        ## implement the possibility to make different layouts and to add attributes to the edges and to the nodes
        EnrichmentDataframe.dropna(subset = 'p_value',inplace = True)
        ig_subgraph=ig.Graph.from_networkx(SubGraph)
        pos_= dict(zip([v['_nx_name'] for v in ig_subgraph.vs],[coord for coord in ig_subgraph.layout_auto()]))
        app=dash.Dash(__name__)
        cyto_node_data=list(
                            zip(
                                    pos_.keys(),
                                    [coord[0] for coord in pos_.values()],
                                    [coord[1] for coord in pos_.values()]
                                    
                                    
                                )
                        )
        nodes = [
        {
            'data': {'id': str(_id), 'label':str(_id)},
            'position': {'x': 120*x, 'y': 120*y}
        }
        for _id, x, y, in cyto_node_data
        ]



        edges = [
            {'data': {'source': source, 'target': target}}
            for source, target in SubGraph.edges()
            ]

        elements = nodes + edges




        default_stylesheet = [
            {
                'selector': 'node',
                'style': {
                    'background-color': '#F5CEC5',
                    'border-color':'black',
                    'border-width':'1',
                    'label': 'data(label)',
                    'width':'60',
                    'height':'60'
                }
            },
            {
                'selector': 'edge',
                'style': {
                    'line-color': 'red',
                    'width':'1'
                }
            }
        ]

        app.layout=html.Div([
                            html.Header(html.H1(['Function enrichment analysis topology visualization'],
                            style={'textAlign':'center','paddingBottom':'50px','border':'0px solid','border-bottom':'1px solid black'})),

                            html.Main([html.Div([html.Label('P-value Slider'),
                                                dcc.Slider(id='pvalue_slider',
                                                        min=round(-np.log10(EnrichmentDataframe['p_value'].max())),
                                                        max=round(-np.log10(EnrichmentDataframe['p_value'].min())),
                                                        value=round(-np.log10(EnrichmentDataframe['p_value'].max())),
                                                        marks=dict(list(zip(set(sorted([round(el) for el in -np.log10(EnrichmentDataframe.p_value.tolist())])),
                                                        [{} for value in set([round(el) for el in -np.log10(EnrichmentDataframe.p_value.tolist())])]))),
                                                            step=None),
                                                    
                                                    html.Div(id='updatemode-output-container', style={'marginTop': 20}),
                                                    
                                                    html.Br(style={'lineHeight':'4'}),
                                                    html.Label('Sources'),
                                                    dcc.RadioItems(id='sources',
                                                                    labelStyle={'display': 'flex'}
                                                                    ),
                                                    html.Br(style={'lineHeight':'4'}), 
                                                    html.Label('Function'),
                                                    dcc.Dropdown(id='function_dropdown'),
                                                    html.P(id='cytoscape-mouseoverNodeData-output')
                                                    


                                                    
                                                    
                                                ],
                                                style={'width':'20%','display':'inline-block','float':'left','paddingTop':'20px','paddingLeft':'50px'}
                                            ),
                                    
                                    
                                    
                                    html.Div([cyto.Cytoscape(id='cytoscape_network',
                                                            layout={'name': 'preset'},
                                                            style={'width': '100%', 'height': '800px'},
                                                            stylesheet=default_stylesheet,
                                                            elements=elements,
                                                            autoRefreshLayout=True
                                                            )
                                                
                                                ],
                                    style={'width':'75%','float':'right','position':'relative','top':'20px'}
                                            
                                            )])
                            
                            ])





        @app.callback(Output('updatemode-output-container', 'children'),
                    Input('pvalue_slider', 'value'))
        def display_value(value):
            return '-log10(P_Value): %s' %value  

        @app.callback(
            Output('sources', 'options'),
            Input('pvalue_slider', 'value'))
        def set_sources(selected_pvalue):
            return [{'label': i, 'value': i} for i in set(EnrichmentDataframe[-np.log10(EnrichmentDataframe.p_value)>=selected_pvalue].source.tolist())]


        @app.callback(Output('function_dropdown', 'options'),
                    Input('pvalue_slider', 'value'),
                    Input('sources', 'value'))
        def set_functions(p_value,source):
            return [{'label': i, 'value': i} for i in set(EnrichmentDataframe[(-np.log10(EnrichmentDataframe.p_value)>=p_value)&(EnrichmentDataframe.source==source)].name.tolist())]


        @app.callback(Output('cytoscape_network', 'stylesheet'),
                    Input('sources', 'value'),
                    Input('function_dropdown', 'value'))



        def update_network(fsource,ffunction):
            """Filter the functions in the dataset"""
            try:
                filt_enrich=EnrichmentDataframe[(EnrichmentDataframe.name==ffunction)&(EnrichmentDataframe.source==fsource)].intersections.values[0]

                new_stylesheet=[{
                                    'selector':"[id='%s']"%ele ,
                                    'style': {
                                        'background-color': 'black',
                                        'line-color': 'black'
                                    }
                                    } for ele in filt_enrich]
            
                return default_stylesheet+new_stylesheet
            except:
                return default_stylesheet
        

        
        
        app.run_server()
    
    EdgeFile=ParseInteractome()

    HippieNet = nx.from_pandas_edgelist(EdgeFile)
    SubG = HippieNet.subgraph(ListOfGenes)
    Enrichment = EnrichmentAnalisys(ListOfGenes)
    BuildGraph(SubG,Enrichment)


#use this function#


def plot_enrichment_analisys_network(list_of_genes,pvalue,colormap='cividis',edgecolor='red',mkcolor='grey',mkfsize=10000,layout='spring',
                                     mklinewidths=2,alpha=1,figsize=(40,20),savefig=False,factor=1,k=10,
                                     cbarfontsize=10,labelling=True,legend=False, legend_fontsize = 20, legend_titlefontsize = 25,
                                     legend_location = (0.5,0.0), legend_col = 6, legend_labelspacing = 1.5, legend_title = '',
                                     legend_columnspacing=1.5, legend_handlelength = 3, size_legend_nofelements=3, cbar_orientation= 'horizontal',
                                     cbar_loc=(1, 0.5),**kwargs):
    def labelling_without_overlapping(x,y,list_of_annotations,ax,verbose=False,**kwargs):
    
        class Point:
            def __init__(self, x, y):
                self.x = x
                self.y = y
    
    
        def doOverlap(ret1,ret2):
            l1 = Point(ret1[0,0],ret1[1,1])
            r1 = Point(ret1[1,0],ret1[0,1])
            l2 = Point(ret2[0,0],ret2[1,1])
            r2 = Point(ret2[1,0],ret2[0,1])

            # If one rectangle is on left side of other
            if l1.x >= r2.x or l2.x >= r1.x:
                return False

            # If one rectangle is above other
            if(r1.y >= l2.y or r2.y >= l1.y):
                return False

            return True

        annotations_coord=[]
        for i, dot in enumerate(y):
            x_coords=x[i]
            y_coords=y[i]
            annotation=ax.annotate(str(list_of_annotations[i]),
                                    xy=(x[i],y[i]),
                                    xytext=(x_coords,y_coords),
                                        **kwargs)

            ax.figure.canvas.draw()
            bbox=matplotlib.text.Text.get_window_extent(annotation)
            bbox_data = ax.transData.inverted().transform(bbox)
            factor=0.2*(bbox_data[0,0]-bbox_data[1,0])
            annotations_coord.append(bbox_data)
            ##BUILD THE SPIRAL##
            theta=np.radians(np.linspace(1,360*200,500))
            r=np.linspace(0,max(max(zip(x,y))),len(theta))
            x_2 = r*np.cos(theta)+x_coords#move the spiral onto the data point
            y_2 = r*np.sin(theta)+y_coords
            n=0
            keep_cycling=True
            while keep_cycling:
                keep_cycling=False
                if verbose==True:
                    print('start checking box %s'% i)
                for ind, box in enumerate (annotations_coord[0:-1]):
                    if verbose:
                        print('checking %s and %s' % (i,ind))
                    if doOverlap(box,bbox_data):
                        if verbose:
                            print('%s and %s overlap' % (i,ind))
                        annotation.set_x(x_2[n])
                        annotation.set_y(y_2[n])
                        n+=1
                        ax.figure.canvas.draw()
                        bbox=matplotlib.text.Text.get_window_extent(annotation)
                        bbox_data = ax.transData.inverted().transform(bbox)
                        annotations_coord.pop()
                        annotations_coord.append(bbox_data)
                        if verbose:
                            print('new coords (x=%i,y=%i)'%(x_coords,y_coords))
                            print('new bbox data',bbox_data)
                            print('annotation coordinates',box)
                            print('restart iteration')
                        keep_cycling=True
                        break


    gp = GProfiler(return_dataframe=True)
    df=gp.profile(organism='hsapiens',
                    query=list_of_genes,
                    significance_threshold_method='bonferroni',
                    no_iea=True,
                    no_evidences=False)
    maxpv=max([-np.log10(p) for p in df.p_value.tolist()])
    for i, (s,v) in enumerate(zip(df.source.value_counts().index,df.source.value_counts())):
        data=df[(df.source==s)&(-np.log10(df.p_value)>pvalue)].reset_index()
        if data.shape[0]==0:
            continue
        else:
            
            nxen=nx.Graph()
            #add nodes
            for i,r in data.iterrows():
                 nxen.add_node(r['name'],size=r['intersection_size'],pvalue=-np.log10(r['p_value']))

            #add edges
            for i,r in data.iterrows():
                for index,row in data.iloc[i+1:].reset_index().iterrows():
                    if len(set(r['intersections']).intersection(set(row['intersections'])))>0:
                        nxen.add_edge(r['name'],
                                    row['name'], 
                                    weight= len(set(r['intersections']).intersection(set(row['intersections']))))
            # Get positions for the nodes in G
            if layout=='spring':
                pos_ = nx.spring_layout(nxen,k)
            
            elif layout=='auto':
                ig_subgraph=ig.Graph.from_networkx(nxen)
                pos_= dict(zip([v['_nx_name'] for v in ig_subgraph.vs],[coord for coord in ig_subgraph.layout_auto()]))
                    
                    


            

            #Normalize connections
            connections=[]
            for edge in nxen.edges(data=True):
                connections.append(edge[2]['weight'])
            if len(connections)!=0:
                
                if ((max(connections)-min(connections)==0) | (len(connections)==0)):
                    norm_connections=[x/100 for x in connections]
                else:
                    norm_connections=[(x-min(connections))/(max(connections)-min(connections)) for x in connections]
            else:
                connections=norm_connections


            #Normalize sizes
            markers=[]
            for node in nxen.nodes(data=True):
                markers.append(node[1]['size'])
            if len(markers)!=0:
                
                if ((max(markers)-min(markers)==0) | (len(markers)==0)):
                    norm_markers=[x/100 for x in markers]
                else:
                    norm_markers=[(x-min(markers))/(max(markers)-min(markers)) for x in markers]
            else:
                markers=norm_markers
            
               
            norm_markers=np.clip(norm_markers,0.3, 1)
            
            
            fig,ax=plt.subplots(figsize=figsize)
            
            ##Plot the nodes
            xses,yses=[],[]
            lab=[]
            colors=[]
            for node in nxen.nodes(data=True):
                xses.append(pos_[node[0]][0])
                yses.append(pos_[node[0]][1])
                lab.append(node[0])
                colors.append(node[1]['pvalue'])
            
            nodez_for_legend = ax.scatter(xses,yses,s=markers)
            nodez=ax.scatter(xses,yses,s=[mkfsize*size for size in norm_markers],
                           c=colors,cmap=colormap,vmax=maxpv,alpha=alpha,edgecolors=mkcolor,
                             linewidths=mklinewidths,clip_on=False,zorder=1)

            ##Mark the labels
            if labelling:
                labelling_without_overlapping(xses,yses,lab,ax,**kwargs)
            


            ##Plot the edges
            for indx, edge in enumerate(nxen.edges(data=True)):
                if edge[2]['weight'] > 0:
                    path_1 = edge[0]#prepare the data to insert in make edge
                    path_2 = edge[1]
                    x0, y0 = pos_[path_1]
                    x1, y1 = pos_[path_2]
                    edgez=ax.plot(np.linspace(x0,x1),np.linspace(y0,y1),
                            color=edgecolor,
                            linewidth = 3*norm_connections[indx]**4,
                                 zorder=0)


            cbar=plt.colorbar(nodez,ax=ax,orientation=cbar_orientation,panchor=cbar_loc)
            cbar.set_label(r'$-log_{10}(p-value)$',fontsize=cbarfontsize+4)
            cbar.ax.tick_params(labelsize=cbarfontsize)
            
            
            if legend:
                handles, _ = nodez.legend_elements(prop="sizes", alpha=0.6, num = size_legend_nofelements)
                _, label_markers = nodez_for_legend.legend_elements(prop="sizes", alpha=0.6)


                legend = ax.legend(handles, label_markers,fontsize = legend_fontsize,
                                    bbox_to_anchor=legend_location,ncol=legend_col,labelspacing = legend_labelspacing,
                                   columnspacing=legend_columnspacing,handlelength=legend_handlelength,frameon = False)

                legend.set_title(legend_title,prop={'size':legend_titlefontsize})

            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.set_title(s,fontsize=35)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            axis = plt.gca()
            # maybe smaller factors work as well, but 1.1 works fine for this minimal example
            axis.set_xlim([factor*x for x in axis.get_xlim()])
            axis.set_ylim([factor*y for y in axis.get_ylim()])
            plt.tight_layout()
            if savefig:
                plt.savefig(str(s)+'enrichment_analysis.jpeg', dpi=300,bbox_inches='tight')
            
            
            
            plt.show()

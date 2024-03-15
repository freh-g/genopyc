import pandas as pd
import matplotlib.pyplot as plt
import requests
from collections import Counter
import ast

def variant_effect_predictor(idlist,input_type = 'rsid', chunked=False,chunksize=200,verbose=False,all_data=False,plot = False,save_plot=''):
    
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

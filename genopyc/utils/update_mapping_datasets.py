import wget
import os

def update_mapping_datasets():
    
    """
    Updates mapping datasets required for gene mapping, Downloads NCBI and UniProt mapping files.

    """
    
    location = os.path.dirname(os.path.realpath(__file__))
    location = location.rstrip('/utils')
    out = os.path.join(location, 'data')
    files = os.listdir(out)
    for f in files:
        if '.gz' in f:
            os.remove(out+'/'+f)
    ncbidb='https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
    tab_uni='https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz'	 
    print('Downloading NCBI mapping file\n')
    wget.download(ncbidb,out)
    print('\nDownloading UniProt mapping file ')
    wget.download(tab_uni,out)
    print('\n!!DATABASES UPDATED!!')

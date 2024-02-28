import wget


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

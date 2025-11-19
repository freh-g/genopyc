import requests
import re
from genopyc.genomic_features.get_variants_info import get_variants_info

def variantId_mapping(list_of_variants, source='variantid', target='rsid'):
    """
    Convert genetic variant identifiers between multiple formats.

    This function allows conversion of genetic variant identifiers across
    several commonly used formats, using parsing logic, Ensembl REST API calls,
    and variant metadata from `genopyc.genomic_features.get_variants_info`.

    Supported identifier types:
        - **variantid** : "CHR_POS_REF_ALT" format (e.g., "1_55516888_T_C")
        - **rsid / rsids** : dbSNP identifiers (e.g., "rs12345")
        - **hgvs** : HGVS genomic notation (e.g., "NC_000001.11:g.55516888T>C")

    Supported conversion directions:
        - variantid → rsid
        - rsid → variantid
        - variantid → hgvs
        - hgvs → variantid
        - hgvs → rsids
        - rsids → hgvs

    Parameters
    ----------
    list_of_variants : list of str
        List of variant identifiers to convert. The format depends on `source`.
    source : {'variantid', 'rsid', 'rsids', 'hgvs'}, default 'variantid'
        The format of the input variants.
    target : {'variantid', 'rsid', 'rsids', 'hgvs'}, default 'rsid'
        The desired output format.

    Returns
    -------
    dict
        Standardized output mapping each input identifier to its converted value(s).
        Depending on the conversion:
            - variantid → rsid: dict {variantid → rsid or None}
            - rsid → variantid: dict {rsid → list of variantid(s)}
            - hgvs → rsids: dict {hgvs → rsid(s)}
            - rsids → hgvs: dict {rsid → list of hgvs strings}
            - variantid → hgvs: dict {variantid → hgvs string}
            - hgvs → variantid: dict {hgvs → variantid string}

        Failed conversions return `None` or an explanatory string.

    Notes
    -----
    - Uses Ensembl REST API `/variant_recoder/homo_sapiens` to resolve rsIDs.
    - Uses `get_variants_info` to extract genomic positions when converting rsid → variantid or rsids → hgvs.
    - Handles multi-allelic variants by creating separate entries for each ALT allele.
    - Conversion to HGVS relies on a predefined GRCh38 RefSeq mapping for common chromosomes.
    - Non-standard chromosomes or malformed inputs are skipped or return descriptive errors.

    Examples
    --------
    >>> variantId_mapping(["1_55516888_T_C"], source="variantid", target="rsid")
    {'1_55516888_T_C': 'rs123456'}

    >>> variantId_mapping(["rs123456"], source="rsid", target="variantid")
    {'rs123456': ['1_55516888_T_C']}

    >>> variantId_mapping(["NC_000001.11:g.55516888T>C"], source="hgvs", target="variantid")
    {'NC_000001.11:g.55516888T>C': '1_55516888_T_C'}
    """

    def extract_positions(variant_dict):
        """
        Extracts chr, pos, ref, alt for standard chromosomes,
        supporting multi-allelic SNPs (e.g., T/A/C/G).
        
        Returns:
            List of tuples: [(chr, pos, ref, alt), ...]
        """
        results = []
        standard_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
        
        for mapping in variant_dict.get('mappings', []):
            chr_ = mapping.get('seq_region_name')
            if chr_ not in standard_chromosomes:
                continue  # skip non-standard chromosomes
            
            pos = mapping.get('start')
            allele_string = mapping.get('allele_string')
            
            if allele_string and pos:
                alleles = allele_string.split('/')
                if len(alleles) >= 2:
                    ref = alleles[0]
                    alts = alleles[1:]
                    for alt in alts:
                        results.append((str(chr_), str(pos), str(ref), str(alt)))
                    results = ['_'.join(x) for x in results]

                else:
                    # caso raro, solo un allele
                    results.append((chr_, pos, alleles[0], None))
                    results = ['_'.join(x) for x in results]
        
        return results


    def variantid2hgvs(lov):
        # Dictionary mapping chromosomes to RefSeq genomic accessions (GRCh38)
        chr_to_refseq = {
            '1': 'NC_000001.11', '2': 'NC_000002.12', '3': 'NC_000003.12', '6': 'NC_000006.12',
            '8': 'NC_000008.11', '10': 'NC_000010.11', '11': 'NC_000011.10', '12': 'NC_000012.12',
            '13': 'NC_000013.11', '15': 'NC_000015.10', '16': 'NC_000016.10', '19': 'NC_000019.10'
        }

        # Convert to HGVS genomic notation
        hgvs_list = []
        for var in lov:
            chrom, pos, ref, alt = var.split('_')
            if chrom in chr_to_refseq:
                hgvs = f"{chr_to_refseq[chrom]}:g.{pos}{ref}>{alt}"
                hgvs_list.append(hgvs)
            else:
                hgvs_list.append(f"Unknown chromosome {chrom}")
        res = dict(zip(lov, hgvs_list))
        return res

    def hgvs2variantid(hgvs_list):
        # Dictionary mapping chromosomes to RefSeq genomic accessions (GRCh38)
        chr_to_refseq = {
            '1': 'NC_000001.11', '2': 'NC_000002.12', '3': 'NC_000003.12', '6': 'NC_000006.12',
            '8': 'NC_000008.11', '10': 'NC_000010.11', '11': 'NC_000011.10', '12': 'NC_000012.12',
            '13': 'NC_000013.11', '15': 'NC_000015.10', '16': 'NC_000016.10', '19': 'NC_000019.10'
        }
        refseq_to_chr = {v: k for k, v in chr_to_refseq.items()}

        res = {}


        for hgvs in hgvs_list:
            try:
                refseq, posrefalt = hgvs.split(":g.")

                if refseq not in refseq_to_chr:
                    res[hgvs] = f"Unknown RefSeq {refseq}"
                    continue

                chrom = refseq_to_chr[refseq]

                # Match typical "123456A>T"
                m = re.match(r"(\d+)([ACGT]+)>([ACGT]+)", posrefalt)
                if not m:
                    res[hgvs] = f"Cannot parse {hgvs}"
                    continue

                pos, ref, alt = m.groups()
                variantid = f"{chrom}_{pos}_{ref}_{alt}"

                res[hgvs] = variantid

            except Exception as e:
                res[hgvs] = f"Error parsing {hgvs}: {e}"

        return res

    def hgvs2rsids(hgvs_list):
        """
        Convert a list of HGVS genomic variants to rsIDs using Ensembl REST API.
        """
        server = "https://rest.ensembl.org"
        ext = "/variant_recoder/homo_sapiens"
        headers = {"Accept": "application/json"}  # Content-Type handled by json=

        # Make POST request
        r = requests.post(server + ext, headers=headers, json={"ids": hgvs_list})

        # Check for errors
        if not r.ok:
            r.raise_for_status()

        # Decode JSON response
        decoded = r.json()
        mapping_dict = {list(element.values())[0]['input']:list(element.values())[0]['id'][0] for element in decoded}
        converted_ids = list(map(mapping_dict.get,hgvs_list))
        res = dict(zip(hgvs_list,converted_ids))
        
        return res


    def variantid2rsid(lovids):
        hgvs_intermediate = variantid2hgvs(lovids)
        rsids = hgvs2rsids(list(hgvs_intermediate.values()))
        final_res = {}
        for var in lovids:
            try:
                final_res[var] = rsids[hgvs_intermediate[var]]
            except:
                final_res[var] = None
        
        return final_res
    
    if (source == 'variantid') & (target == 'rsid'):
        results = variantid2rsid(list_of_variants)
        
    elif (source == 'rsid') & (target == 'variantid'):
        MappingDict = {}
        info = get_variants_info(list_of_variants)
        
        for v in list_of_variants:
            try:
                variant_positions = extract_positions(info[v])
                MappingDict[v] = variant_positions
            except Exception as e:
                print(e, f"Couldn't Convert Variant {v}")
        mapped_variants = list(map(MappingDict.get, list_of_variants))
        results = list(zip(list_of_variants, mapped_variants))
    
    elif (source == 'hgvs') & (target == 'rsids'):
        results = hgvs2rsids(list_of_variants)
    
    elif (source == 'rsids') & (target == 'hgvs'):
        results = {}
        info = get_variants_info(list_of_variants)
        variantids = {v:extract_positions(info[v]) for v in list_of_variants}
        for var in list_of_variants:
            try:
                lofhgvs = list(variantid2hgvs(variantids[var]).values())
                results[var] = lofhgvs
            
            
            except Exception as e:
                print(e, f"Couldn't Convert Variant {rsid}")
                results[var] = None

    elif (source == 'variantid') & (target == 'hgvs'):
        results = variantid2hgvs(list_of_variants)
    
    elif (source == 'hgvs') & (target == 'variantid'):
        results = hgvs2variantid(list_of_variants)
        
    
    
    return results

        
    
    
    

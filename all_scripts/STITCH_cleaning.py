import pandas as pd
import urllib.parse
import urllib.request
import io


def human_check(uniprot_id_string, from_db, to_db):
    '''
    :param uniprot_id_string: string with the ID to map to Uniprot
    :param from_db: type of identifier (from what database comes from -> ID (Uniprot ID), STRING_ID
    :param to_db: Map to which database (ID ->To Uniprot ID with gene name; ACC -> To Uniprot Accession number)
    :return: Dataframe containing in one column the input IDs and the other the respective mapping
    '''
    url = 'https://www.uniprot.org/uploadlists/'

    params = {
        'from': from_db,
        'to': to_db,
        'format': 'tab',
        'query': uniprot_id_string
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()

    answer = response.decode('utf-8')

    data = io.StringIO(answer)
    dataf = pd.read_csv(data, sep="\t")

    return dataf


def stitch_cleaning(link_file, chemical_file, inchi_file):
    links = pd.read_csv(link_file, sep='\t')  # file containig the relationship between proteins and compounds

    links_cleaned = links[links['combined_score'] >= 800]  # Stitch cutoff - higher more reliable

    association = pd.DataFrame()
    chunksize = 10 ** 6
    for chunk in pd.read_csv(chemical_file, chunksize=chunksize, sep='\t', usecols=['chemical',
                                                                                    'name',
                                                                                    'SMILES_string']):
        # chemical file is too big to be loaded all together, so we load it in batch and map it to the link file
        association = pd.concat([association, links_cleaned.merge(chunk, how='inner', on='chemical')])

    association = association[['chemical', 'protein', 'name', 'SMILES_string']].drop_duplicates()
    association['name'] = association.name.str.upper()
    uniprot_ids_string = " ".join(list(set(association['protein'].to_list())))
    mapped_uniprots = human_check(uniprot_ids_string, 'STRING_ID', 'ACC')
    association = association.merge(mapped_uniprots, left_on='protein', right_on='From', how='inner')[['chemical', 'name', 'To', 'SMILES_string']]

    # map to add the inchi key, as before the file is too big to added all together, and, even worse, the keys are
    # splitted in to columns instead of one, so is needed a double mapping

    res_inchi = pd.DataFrame()
    res_inchi_2 = pd.DataFrame()
    for chunk in pd.read_csv(inchi_file, chunksize=chunksize, sep='\t'):
        res_inchi = pd.concat(
            [res_inchi, association.merge(chunk, how='inner', left_on='chemical', right_on='flat_chemical_id')])
        res_inchi_2 = pd.concat(
            [res_inchi_2, association.merge(chunk, how='inner', left_on='chemical', right_on='stereo_chemical_id')])

    final = res_inchi.append(res_inchi_2, ignore_index=True)[['chemical', 'name', 'To', 'inchikey', 'SMILES_string']].drop_duplicates()

    final['Database_STITCH'] = 'STITCH'

    final = final.rename(columns={'name': 'compound_name',
                                  'To': 'target_id',
                                  'inchikey': 'standard_inchi_key'})

    final.to_csv('STITCH_cleaned.input', sep='\t', index=False)


stitch_cleaning('STITCH/9606.protein_chemical.links.v5.0.tsv',
                'STITCH/chemicals.v5.0.tsv',
                'STITCH/chemicals.inchikeys.v5.0.tsv')


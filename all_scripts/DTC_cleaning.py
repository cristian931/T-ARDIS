#!/usr/bin/env python
# coding: utf-8

"""This script is used to clean the drug-target database, the idea is to have as final output a pairwise
relationship between the drug and the targets as Uniprot ids filtered out for significant activation assays like the
IC50 value"""


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

#### DTC #####

def DTC_cleaning(DTC_file):

    df = pd.read_csv(DTC_file,
                 usecols=['standard_inchi_key',
                          'compound_name',
                          'target_id',
                          'standard_type',
                          'standard_relation',
                          'standard_value',
                          'standard_units'])

    type_list = ['IC50', 'EC50', 'POTENCY']

    dfIC50 = df[df['standard_type'].isin(type_list)]
    dfIC50.dropna(subset=['standard_relation'], inplace=True)

    # remove entries that indicates values above a threshold
    dfIC50 = dfIC50[~dfIC50['standard_relation'].str.contains('>')]

    dfIC50_active = dfIC50[(dfIC50['standard_value'] <= 100)
                       &
                       (dfIC50['standard_units'] == 'NM')
                       &
                       (dfIC50['standard_relation'].str.contains('=|<'))]  # remove other possible operands (like ~)

    dfIC50_active = dfIC50_active[['standard_inchi_key',
                                   'compound_name',
                                   'target_id']]

    dfIC50_active.dropna(subset=['target_id'], inplace=True)

    dfIC50_active['target_id'] = dfIC50_active['target_id'].apply(lambda x: x.split(', '))  # some uniprot id are listed
    # together as strings, transform them into list

    dfIC50_active = dfIC50_active.explode('target_id')  # get single pairwise drug uniprot id relationship
    # for those listed together

    inchi_key_grouped = dfIC50_active.groupby(['standard_inchi_key', 'compound_name'], dropna=False).agg(set).reset_index()  # Collapse drugs with the same INCHI key (same compound)

    inchi_key_grouped.dropna(subset=['compound_name'], inplace=True)

    inchi_key_grouped['len'] = inchi_key_grouped['compound_name'].apply(lambda x: len(x))  # Check how many different drugs name map on the same INCHI key

    final = inchi_key_grouped[~(inchi_key_grouped['compound_name'] == 'None')]
    # Be sure of removing possible missing entries and drugnames
    final = final.drop(columns=['len']).explode('target_id')  # obtain drug-target pairwise
    # relationship
    final['Database_DTC'] = 'Drug_Target_Commons'

    print(final)

    # check the human proteins

    uniprot_ids_string = " ".join(list(set(final['target_id'].to_list())))

    mapped_uniprots = human_check(uniprot_ids_string, 'ID', 'ID')

    final = final.merge(mapped_uniprots, left_on='target_id', right_on='From', how='inner')

    final = final[final['To'].str.contains('HUMAN')].drop_duplicates()

    final.to_csv('DTC_cleaned.input')


DTC_cleaning('DRUG_TARGETS_COMMONS/DTC_data.csv')


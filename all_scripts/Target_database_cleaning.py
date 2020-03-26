#!/usr/bin/env python
# coding: utf-8

"""This script is used to clean the different drug-target databases, the idea is to have as final output a pairwise
relationship between the drug and the targets as Uniprot ids """


import pandas as pd
import numpy as np


# DRUGBANK: From Drugbank two different file have been downloaded, the drugbank vocabulary, containing the drugs
# names and their respective codes and the Target file, containg the information of the drugs' targets and the
# relationship with the former

drug_drugbank_file = pd.read_csv('DRUGBANK/drugbank vocabulary.csv', dtype=object)
drug_drugbank = drug_drugbank_file[['DrugBank ID', 'Common name']].rename(columns={'DrugBank ID': 'Drug IDs'})


target_drugbank_file = pd.read_csv('DRUGBANK/all.csv', dtype=object)

target_drugbank = target_drugbank_file[['UniProt ID', 'Drug IDs']]

target_drugbank['Drug IDs'] = target_drugbank['Drug IDs'].apply(lambda x: x.split(';'))

target_exploded = target_drugbank.explode('Drug IDs')

drugbank_merged = pd.merge(drug_drugbank, target_exploded, how='left', on='Drug IDs').dropna()

drugbank_merged['Common name'] = drugbank_merged['Common name'].str.upper()

drugbank_final = drugbank_merged.drop_duplicates()
drugbank_final.to_csv('Drugbank_relationship.input', sep='\t', index=False)

########################################################################################################################

# DGiDB: unique File from which extract the information. However the conversion from Gene name to UNIPROT id is
# needed. the latter has been downloaded from uniprot searching for the human proteome and selecting the uniprot
# entries together with the gene names and also the gene names synonyms

dgidb_file = pd.read_csv('DGidb/interactions.tsv', sep='\t')

dgidb_interactions = dgidb_file[['drug_claim_primary_name','gene_claim_name']].dropna()
dgidb_interactions['drug_claim_primary_name'] = dgidb_interactions['drug_claim_primary_name']\
                                                .apply(lambda x: x.split(', '))

dgidb_interactions = dgidb_interactions.explode('drug_claim_primary_name')

dgidb_interactions['drug_claim_primary_name'] = dgidb_interactions['drug_claim_primary_name'].str.upper()


uniprot = pd.read_csv('human_proteome_gene_names_11_03.tab',
                      sep='\t',
                      names=['Uniprot ID', 'Primary Gene', 'Synonym'],
                      header=0).fillna('')

uniprot['Gene names'] = uniprot[['Primary Gene','Synonym']].agg(' '.join, axis=1)
uniprot['Gene names'] = uniprot['Gene names'].apply(lambda x: x.strip().split())

uniprot_exploded = uniprot[['Uniprot ID', 'Gene names']].explode('Gene names')


dgidb_merged = pd.merge(dgidb_interactions,
                        uniprot_exploded,
                        left_on=['gene_claim_name'],
                        right_on=['Gene names'],
                        how='left'
                        )\
                        .dropna().drop(columns=['gene_claim_name'])

dgidb_final = dgidb_merged.drop_duplicates()

dgidb_final.to_csv('DGidb_relationships.input', sep='\t', index=False)

########################################################################################################################

# DRUGCENTRAL

drugcentral_file = pd.read_csv('DRUG_CENTRAL/drug.target.interaction.tsv', sep='\t', dtype=object)

drugcentral = drugcentral_file.loc[drugcentral_file['ORGANISM'] == 'Homo sapiens'][['DRUG_NAME', 'ACCESSION', 'GENE']]

drugcentral['DRUG_NAME'] = drugcentral['DRUG_NAME'].str.upper()

drugcentral['ACCESSION'] = drugcentral['ACCESSION'].apply(lambda x: x.split('|'))

drugcentral['GENE'] = drugcentral['GENE'].apply(lambda x: x.split('|'))

drugcentral_final = drugcentral.explode('ACCESSION')

drugbank_final.to_csv('drugcentral_merged.input', sep='\t', index=False)


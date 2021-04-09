#!/usr/bin/env python
# coding: utf-8

# This three function are used to clean the drug se database MEDEFFECT, SIDER and off side, The faers will have a
# different procedure. All three need a slightly different approach since the different database construction
import pandas as pd

# MEDEFFECT: The info are separated in two file, the first contain the drug name and id, the second the side
# effect with the identifying id and the side effect name in english and french and the respective system affected

header_drug = ['REPORT_DRUG_ID',
               'REPORT_ID',
               'DRUG_PRODUCT_ID',
               'DRUGNAME',
               'DRUGINVOLV_ENG',
               'DRUGINVOLV_FR',
               'ROUTEADMIN_ENG',
               'ROUTEADMIN_FR',
               'UNIT_DOSE_QTY',
               'DOSE_UNIT_FR',
               'DOSE_UNIT_ENG',
               'FREQUENCY',
               'FREQ_TIME',
               'FREQUENCY_TIME_ENG',
               'FREQUENCY_TIME_FR',
               'FREQ_TIME_UNIT_ENG',
               'FREQ_TIME_UNIT_FR',
               'THERAPY_DURATION',
               'THERAPY_DURATION_UNIT_ENG',
               'THERAPY_DURATION_UNIT_FR',
               'DOSAGEFORM_ENG',
               'DOSAGEFORM_FR'
               ]

header_reaction = ['REACTION_ID',
                   'REPORT_ID',
                   'DURATION',
                   'DURATION_UNIT_ENG',
                   'DURATION_UNIT_FR',
                   'PT_NAME_ENG',
                   'PT_NAME_FR',
                   'SOC_NAME_ENG',
                   'SOC_NAME_FR',
                   'MEDDRA_VERSION'
                   ]

header_cleaned_drugs = ['DRUG_PRODUCT_ID',
                        'STANDARD_CONCEPT_ID',
                        'DRUGNAME_CLEANED']

medeffect_drugs = pd.read_csv('MEDEFFECT/report_drug.txt', sep='$', names=header_drug, dtype=object)

side_effect_medeffect = pd.read_csv('MEDEFFECT/reactions.txt', sep='$', names=header_reaction, dtype=object)

drugs_cleaned = pd.read_csv('MEDEFFECT/MEDEFFECT_DRUG_CLEANED.csv', sep=',', names=header_cleaned_drugs, dtype=object)

# remove entries with multiple drugs at the same time
drugs_cleaned = drugs_cleaned[~drugs_cleaned['DRUGNAME_CLEANED'].str.contains('/')]
drugs_cleaned['DRUGNAME_CLEANED'] = drugs_cleaned['DRUGNAME_CLEANED'].str.replace(' HYDROCLORIDE', '')

medeffect_report_related = pd.merge(drugs_cleaned[['DRUG_PRODUCT_ID', 'DRUGNAME_CLEANED']],
                                    medeffect_drugs[['DRUG_PRODUCT_ID', 'REPORT_ID']],
                                    how='inner',
                                    on='DRUG_PRODUCT_ID')

medeffect_ADR_related = pd.merge(medeffect_report_related[['REPORT_ID', 'DRUGNAME_CLEANED']],
                                 side_effect_medeffect[['REPORT_ID', 'PT_NAME_ENG']],
                                 how='inner',
                                 on='REPORT_ID')

medeffect_ADR_related['DRUGNAME_CLEANED'] = medeffect_ADR_related['DRUGNAME'].str.upper()
medeffect_ADR_related['PT_NAME_ENG'] = medeffect_ADR_related['PT_NAME_ENG'].str.capitalize()

medeffect_related = medeffect_ADR_related.rename(columns={'REPORT_ID': 'MEDEFFECT_REPORT_ID'}).drop_duplicates()

medeffect_related['Database'] = 'MEDEFFECT'

medeffect_related.to_csv('MEDEFFECT_DRUG_SE.input', sep='\t', index=False)

# SIDER: The file are separated in two file, drug_namse.tsv and meddra_all_se.tsv. Using the unique ID we're able to
# map the drug to their se


sider_drug = pd.read_csv('SIDER_4.1/drug_names.tsv', sep='\t', names=['STICH_ID_1', 'DRUGNAME'], dtype=object)

sider_side_effect = pd.read_csv('SIDER_4.1/meddra_all_se.tsv',
                                sep='\t',
                                names=['STICH_ID_1',
                                       'STICH_ID_2',
                                       'UMLS_ID',
                                       'MEDDRA_CONCEPT',
                                       'UMLS_CONCEPT_ID',
                                       'SIDEEFFECT'
                                       ],
                                dtype=object)


sider_related = pd.merge(sider_drug,
                         sider_side_effect[['STICH_ID_1', 'SIDEEFFECT']],
                         how='left',
                         on='STICH_ID_1'
                         )

sider_related['DRUGNAME'] = sider_related['DRUGNAME'].str.upper()
sider_related['SIDEEFFECT'] = sider_related['SIDEEFFECT'].str.capitalize()

sider_related = sider_related.rename(columns={'STICH_ID_1': 'SIDER_ID'}).drop_duplicates()

sider_related['Database'] = 'SIDER'

sider_related.to_csv('SIDER_DRUG_SE.input', sep='\t', index=False)


# OFFSIDE

offside_db = pd.read_csv('OFFSIDE/OFFSIDES.csv', dtype=object)

offside_db['drug_concept_name'] = offside_db['drug_concept_name'].str.upper()
offside_db['condition_concept_name'] = offside_db['condition_concept_name'].str.capitalize()

offside_clean = offside_db[['drug_rxnorn_id', 'drug_concept_name', 'condition_concept_name']].drop_duplicates()

offside_clean['Database'] = 'OFFSIDE'

offside_clean.to_csv('OFFSIDE_DRUG_SE.input', sep='\t', index=False)

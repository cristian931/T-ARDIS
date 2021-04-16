#!/usr/bin/env python
# coding: utf-8

""" This script is used to combine, merge and extrapolate information from the different drug-se, drug-targets
databases. The aim is to retrieve the pairwise relationship between side effects and drug target and statistically
validate them through fisher exact test and q-value correction """

import pandas as pd
import scipy.stats as stats
from multipy.fdr import qvalue
from pandarallel import pandarallel

pandarallel.initialize()

# Load the drug-se databases previously cleaned

faers = pd.read_csv('relationship_analysis_input_files/Significant_interaction_FAERS.input',
                    sep='\t',
                    dtype=object,
                    usecols=['drugname',
                             'adverse_event',
                             'Database'
                             ]
                    )
faers = faers.rename(columns={'drugname': 'drug',
                              'adverse_event': 'se',
                              'Database': 'Database_FAERS'
                              }
                     ).sort_values('drug')

medeffect = pd.read_csv('relationship_analysis_input_files/Significant_interaction_MEDEFFECT.input',
                        sep='\t',
                        dtype=object,
                        usecols=['drugname',
                                 'adverse_event',
                                 'Database'
                                 ]
                        )
medeffect = medeffect.rename(columns={'drugname': 'drug',
                                      'adverse_event': 'se',
                                      'Database': 'Database_MEDEFFECT'
                                      }
                             ).sort_values('drug')

offside = pd.read_csv('relationship_analysis_input_files/OFFSIDE_DRUG_SE.input',
                      sep='\t',
                      dtype=object,
                      usecols=['drug_concept_name',
                               'condition_concept_name',
                               'Database'
                               ]
                      )
offside = offside.rename(columns={'drug_concept_name': 'drug',
                                  'condition_concept_name': 'se',
                                  'Database': 'Database_OFFSIDE'
                                  }
                         ).sort_values('drug')

sider = pd.read_csv('relationship_analysis_input_files/SIDER_DRUG_SE.input',
                    sep='\t',
                    dtype=object,
                    usecols=['DRUGNAME',
                             'SIDEEFFECT',
                             'Database'
                             ]
                    )
sider = sider.rename(columns={'DRUGNAME': 'drug',
                              'SIDEEFFECT': 'se',
                              'Database': 'Database_SIDER'
                              }
                     ).sort_values('drug')


########################################################################################################################

# Defining the different functions

def meddra_cleaning(dataset):
    # LOAD MEDDRA  DB TO COVERT POSSIBLE LLT TO PT
    meddra_db = pd.read_csv('MeDDRA_complete_LLT', sep=';').drop(columns="Primary SOC").drop_duplicates()

    # Create a dictionary with LLT as keys and PT as values
    meddra_dic = dict(zip(meddra_db['English'], meddra_db['PT']))

    # Map possible LLT in DRUG/ADRs dataframe to PT
    dataset['se'] = dataset['se'].map(meddra_dic).fillna(
        dataset['se'])  # note that there possible be some ADR that are neither PT ot LLT

    # Select only the row whose side effects corresppond to PT in MEDDRA
    dataset = dataset[dataset['se'].isin(set(meddra_db['PT'].to_list()))]

    # Exclude DRUG/ADRs pair if ADRs fall in this particular SOCs
    Excluding_SOC_list = ['General disorders and administration site conditions',
                          'Injury, poisoning and procedural complications',
                          'Investigations',
                          'Neoplasms benign, malignant and unspecified (incl cysts and polyps)',
                          'Product issues',
                          'Social circumstances',
                          'Surgical and medical procedures']

    list_adr_to_remove = meddra_db[meddra_db['SOC'].isin(Excluding_SOC_list)]['PT'].to_list()

    # Exclude ADR being part of particular SOC
    dataset = dataset[~dataset.se.isin(set(list_adr_to_remove))]

    dataset = dataset.groupby('drug').agg(lambda x: list(set(x.tolist()))).reset_index()

    return dataset


def overlap(x):
    """Function to define the overlap between two lists contained in a pandas dataframe
    1) return the length of the drug shared list

    Parameters
    ----------
    x : Dataframe
        Dataframe containing the drug list related to the side effect and the drug list related to the target
    """

    return len(list(set(x['se_drug']) & set(x['tg_drug'])))


def target_se_merging(drug_adr_database, drug_target_database):
    interaction = pd.merge(drug_adr_database,
                           drug_target_database,
                           how='inner',
                           on='drug')

    # Extrapolate the number of drugs that present a particular side effect reconstructing the pairwise relationship
    # between the mapped drugs and sideffect and grouping by the latter. In this way we can compute the number of drugs
    # simply applying the len function.

    dr_se_pairwise = interaction[['drug', 'se']].explode('se').groupby('se') \
        .agg(lambda x: list(set(x.tolist()))).reset_index()

    dr_se_pairwise['se_drug_len'] = dr_se_pairwise['drug'].apply(lambda x: len(x))

    # Same thing can be done for defining the number of drugs that present a particular target.

    dr_tg_pairwise = interaction[['drug', 'target']].explode('target').groupby('target') \
        .agg(lambda x: list(set(x.tolist()))).reset_index()

    dr_tg_pairwise['tg_drug_len'] = dr_tg_pairwise['drug'].apply(lambda x: len(x))

    # extract from dataset the pairwise interaction between side effect and targets, since both are lists, we need first
    # to explode one list (se) extract the column of interest and then explode the target list.

    se_tg_pairwise_part_1 = interaction.explode('se')[['se', 'target']].explode('target')

    se_tg_pairwise_part_2 = pd.merge(se_tg_pairwise_part_1,
                                     dr_se_pairwise.rename(columns={'drug': 'se_drug'
                                                                    }
                                                           ),
                                     how='left',
                                     on='se'
                                     )

    se_tg_pairwise_part_3 = pd.merge(se_tg_pairwise_part_2,
                                     dr_tg_pairwise.rename(columns={'drug': 'tg_drug'
                                                                    }
                                                           ),
                                     how='left',
                                     on='target'
                                     )

    lists = se_tg_pairwise_part_3[['se_drug', 'tg_drug']]  # split the dataframe to reduce memory consumption

    values = se_tg_pairwise_part_3[['se_drug_len', 'tg_drug_len']]

    values['overlap_len'] = lists.parallel_apply(overlap, axis=1)

    # Now define a function for the p-value computation using the fisher exact test.

    interaction_len = len(interaction)  # extract the total number of drugs (len of the interaction dataframe)

    lighter = se_tg_pairwise_part_3[['se', 'target']]

    return interaction, interaction_len, values, lighter


def fisher(x, interaction_len):
    """Function that compute the p-value applying the fisher exact test.
    1) Return the computed p-value

    Parameters
    ----------
    x : Dataframe
        Dataframe containing the following information:

        se_binds : int
            represent the number of drugs that present at the same time the SE AND the target

        se_no_binds : int
            is the number of drugs that present only the side effect

        no_se_binds : int
            is the number of drugs that only bind the target

        no_se_no_binds : int
            is the number of drug that do not present the se and not bind the target

    interaction_len : int
        The total number of drugs found
    """
    se_binds = x['overlap_len']
    se_no_binds = x['se_drug_len'] - x['overlap_len']
    no_se_binds = x['tg_drug_len'] - x['overlap_len']
    no_se_no_binds = interaction_len - se_binds - se_no_binds - no_se_binds
    _, pvalue = stats.fisher_exact([[se_binds, se_no_binds], [no_se_binds, no_se_no_binds]])
    return pvalue


def final_adjustements(drug_adr_database, drug_target_database, database_type):
    # Select the pair column only in order to simplify the computation

    interaction, interaction_len, values, lighter = target_se_merging(drug_adr_database, drug_target_database)

    lighter['pvalue'] = values.parallel_apply(fisher, interaction_len=interaction_len, axis=1)

    lighter.to_csv('p_value_computed' + database_type + '.csv',
                   sep='\t',
                   index=False
                   )

    # obtain the q value correction on the p-value computed

    _, qvals = qvalue(lighter['pvalue'].tolist())

    lighter['qvals'] = qvals

    Final = lighter.drop_duplicates()

    Final.to_csv('qvalues_interactions' + database_type,
                 sep='\t',
                 index=False
                 )

    accepted = Final.loc[Final['qvals'] <= 0.05]
    accepted.to_csv('accepted_interactions' + database_type,
                    sep='\t',
                    index=False
                    )

    return interaction, accepted


def TARDIS_tables(interaction, accepted, database_type):
    exploded_interaction = interaction.explode('se').explode('target')

    drug_tg_se_stats = pd.merge(accepted,
                                exploded_interaction,
                                on=[
                                    'se',
                                    'target'
                                ],
                                how='inner'
                                )

    for tg_dataframe in [drugtargetcommons]:
        drug_tg_se_stats = drug_tg_se_stats.merge(
            tg_dataframe.drop_duplicates(),
            on=[
                'drug',
                'target'
            ],
            how='left'
        )

    if database_type == 'community':
        for se_dataframe in [faers, medeffect]:
            drug_tg_se_stats = drug_tg_se_stats.merge(
                se_dataframe.drop_duplicates(),
                on=[
                    'drug',
                    'se'
                ],
                how='left'
            )
        drug_tg_se_stats.sort_values('drug').to_csv('TARDIS_TG_SE_DRUG_STATS_TABLE_COMMUNITY', sep='\t', index=False)

    else:
        for se_dataframe in [offside, sider]:
            drug_tg_se_stats = drug_tg_se_stats.merge(
                se_dataframe.drop_duplicates(),
                on=[
                    'drug',
                    'se'
                ],
                how='left'
            )
        drug_tg_se_stats.sort_values('drug').to_csv('TARDIS_TG_SE_DRUG_STATS_TABLE_CONTROLLED', sep='\t', index=False)


########################################################################################################################
# Create two distinct datasets, one containing the information derived from community uploaded data (MEDEFFECT, FAERS)
# less controlled, the other from more reliable databases (SIDER, OFFSIDE)

faers_grouped = faers.groupby('drug').agg(set).reset_index()
medeffect_grouped = medeffect.groupby('drug').agg(set).reset_index()

df_se_community = pd.merge(faers_grouped, medeffect_grouped, on='drug', how='inner')
df_se_community['union_se'] = df_se_community.apply(lambda row: row['se_x'].union(row['se_y']), axis=1)
df_se_community = df_se_community[['drug', 'union_se']].explode('union_se')
df_se_community = df_se_community.rename(columns={'union_se': 'se'})

df_se_controlled = offside.append(sider, ignore_index=True)[['drug', 'se']]

# Now we load the datasets regarding drug target relationships
drugtargetcommons = pd.read_csv('relationship_analysis_input_files/Drug_Target_Commons.input',
                                sep='\t',
                                usecols=['compound_name',
                                         'target_id',
                                         'Database'
                                         ]
                                )

drugtargetcommons = drugtargetcommons.rename(columns={'Common name': 'drug',
                                                      'UniProt ID': 'target',
                                                      'Database': 'Database_DRUGBANK'
                                                      }
                                             )

df_target = drugtargetcommons[['drug', 'target']] \
    .groupby('drug') \
    .agg(lambda x: list(set(x.tolist()))) \
    .reset_index()

df_se_community = meddra_cleaning(df_se_community)
df_se_controlled = meddra_cleaning(df_se_controlled)

interaction_community, accepted_community = final_adjustements(df_se_community, df_target, 'community')
interaction_controlled, accepted_controlled = final_adjustements(df_se_controlled, df_target, 'controlled')

TARDIS_tables(interaction_community, accepted_community, 'community')
TARDIS_tables(interaction_controlled, accepted_controlled, 'controlled')

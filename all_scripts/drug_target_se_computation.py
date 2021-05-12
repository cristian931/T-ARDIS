#!/usr/bin/env python
# coding: utf-8

""" This script is used to combine, merge and extrapolate information from the different drug-se, drug-targets
databases. The aim is to retrieve the pairwise relationship between side effects and drug target and statistically
validate them through fisher exact test and q-value correction """

import pandas as pd
import glob
import sys
import scipy.stats as stats
from multipy.fdr import qvalue
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from pandarallel import pandarallel

pandarallel.initialize(progress_bar=True)

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
    meddra_file = glob.glob("**/*/mdhier.asc")
    llt_file = glob.glob("**/*/llt.asc")

    meddra_header = ['pt_code',
                     'hlt_code',
                     'hlgt_code',
                     'soc_code',
                     'pt_name',
                     'hlt_name',
                     'hlgt_name',
                     'soc_name',
                     'soc_abbrev',
                     'null_field',
                     'pt_soc_code',
                     'primary_soc_fg',
                     'tmp']

    llt = pd.read_csv(llt_file[0], sep='$',
                      names=['llt_code', 'llt_name', 'pt_code', 'llt_whoart_code ', 'llt_harts_code ',
                             'llt_costart_sym ', 'llt_icd9_code ', 'llt_icd9cm_code ', 'llt_icd10_code ',
                             'llt_currency', 'llt_jart_code', 'tmp'],
                      usecols=['llt_code', 'llt_name', 'pt_code'])

    meddra_db = pd.read_csv(meddra_file[0], sep='$', names=meddra_header)

    meddra_db = meddra_db[meddra_db['primary_soc_fg'] == 'Y']

    meddra_db_with_LLT = pd.merge(meddra_db, llt, on='pt_code', how='inner')

    # Create a dictionary with LLT as keys and PT as values
    meddra_dic = dict(zip(meddra_db_with_LLT['llt_name'], meddra_db_with_LLT['pt_name']))

    # Map possible LLT in DRUG/ADRs dataframe to PT
    dataset['se'] = dataset['se'].map(meddra_dic).fillna(
        dataset['se'])  # note that there possible be some ADR that are neither PT ot LLT

    # Select only the row whose side effects correspond to PT in MEDDRA
    dataset = dataset[dataset['se'].isin(set(meddra_db_with_LLT['pt_name'].to_list()))]

    # Exclude DRUG/ADRs pair if ADRs fall in this particular SOCs, HLGT and HLT
    Excluding_SOC_list = ['General disorders and administration site conditions',
                          'Injury, poisoning and procedural complications',
                          'Investigations',
                          'Neoplasms benign, malignant and unspecified (incl cysts and polyps)',
                          'Product issues',
                          'Social circumstances',
                          'Surgical and medical procedures',
                          'Infections and infestations',
                          'Psychiatric disorders']

    list_PT_SOC_to_remove = meddra_db[meddra_db['soc_name'].isin(Excluding_SOC_list)]['pt_name'].to_list()
    # Exclude ADR being part of particular SOC
    dataset = dataset[~dataset.se.isin(set(list_PT_SOC_to_remove))]

    dataset = dataset.groupby('drug').agg(set).reset_index()

    return dataset


def overlap(x):
    """Function to define the overlap between two lists contained in a pandas dataframe
    1) return the length of the drug shared list

    Parameters
    ----------
    x : Dataframe
        Dataframe containing the drug list related to the side effect and the drug list related to the target
    """

    return len(x['se_drug'].intersection(x['tg_drug']))


def target_se_merging(drug_adr_database, drug_target_database):
    interaction = pd.merge(drug_adr_database,
                           drug_target_database,
                           how='inner',
                           on='drug')

    # Compute the TANIMOTO Score on mapped drugs to remove the one to similar
    daf_smiles = interaction[['SMILES_string']].drop_duplicates().dropna()
    df_smiles = daf_smiles['SMILES_string']

    c_smiles = []
    for ds in df_smiles:
        try:
            Chem.CanonSmiles(ds)
            c_smiles.append(ds)
        except:
            continue

    ms = [Chem.MolFromSmiles(x) for x in c_smiles]

    fps = [FingerprintMols.FingerprintMol(x) for x in ms]

    qu, ta, sim = [], [], []

    for n in range(len(fps) - 1):  # -1 so the last fp will not be used
        s = DataStructs.BulkTanimotoSimilarity(fps[n], fps[n + 1:])  # +1 compare with the next to the last fp
        # collect the SMILES and values
        for m in range(len(s)):
            qu.append(c_smiles[n])
            ta.append(c_smiles[n + 1:][m])
            sim.append(s[m])

    # build the dataframe
    d = {'query': qu, 'target': ta, 'Similarity': sim}
    tanimoto_scores = pd.DataFrame(data=d)

    tanimoto_to_remove = tanimoto_scores[tanimoto_scores['Similarity'] >= 0.7]

    tanimoto_smiles = list(set(tanimoto_to_remove['query'].to_list()))

    interaction = interaction[~interaction['SMILES_string'].isin(tanimoto_smiles)]

    return interaction
    

def pairwiser(interaction):
    # Extrapolate the number of drugs that present a particular side effect reconstructing the pairwise relationship
    # between the mapped drugs and side effects and grouping by the latter.
    # In this way we can compute the number of drugs simply applying the len function.

    dr_se_pairwise = interaction[['drug', 'se']].explode('se').groupby('se').agg(set).reset_index()

    dr_se_pairwise['se_drug_len'] = dr_se_pairwise['drug'].apply(lambda x: len(x))

    # Same thing can be done for defining the number of drugs that present a particular target.

    dr_tg_pairwise = interaction[['drug', 'target']].explode('target').groupby('target').agg(set).reset_index()

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

    se_tg_pairwise_part_3['overlap_len'] = se_tg_pairwise_part_3.parallel_apply(overlap, axis=1)

    return se_tg_pairwise_part_3


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


def final_adjustements(interaction, database_type):
    
    values_df = pairwiser(interaction)

    values_df['pvalue'] = values_df.parallel_apply(fisher, interaction_len=len(interaction), axis=1)

    values_df[['se', 'target', 'pvalue']].to_csv('p_value_computed_' + database_type + '.csv', sep='\t', index=False)

    # obtain the q value correction on the p-value computed

    _, qvals = qvalue(values_df['pvalue'].tolist())

    values_df['qvals'] = qvals

    Final = values_df[['se', 'target', 'pvalue', 'qvals']].drop_duplicates()

    Final.to_csv('qvalues_interactions_' + database_type,
                 sep='\t',
                 index=False
                 )

    accepted = Final.loc[Final['qvals'] <= 0.05]
    accepted.to_csv('accepted_interactions_' + database_type,
                    sep='\t',
                    index=False
                    )

    return accepted


def TARDIS_tables(interaction, accepted, database_type):
    exploded_interaction = interaction.explode('se').explode('target').drop(columns=['SMILES_string']).drop_duplicates()
    #print(exploded_interaction)

    drug_tg_se_stats = pd.merge(accepted,
                                exploded_interaction,
                                on=[
                                    'se',
                                    'target'
                                ],
                                how='inner'
                                )
    
    print(drug_tg_se_stats)

    for tg_dataframe in [
        dtc[['drug', 'target', 'Database_DTC']].drop_duplicates(),
        stitch[['drug', 'target', 'Database_STITCH']].drop_duplicates()
    ]:
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
                how='inner'
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
        drug_tg_se_stats.sort_values('drug')\
            .dropna(subset=['Database_OFFSIDE', 'Database_SIDER'], how='all')\
            .to_csv('TARDIS_TG_SE_DRUG_STATS_TABLE_CONTROLLED', sep='\t', index=False)


########################################################################################################################
# Create two distinct datasets, one containing the information derived from community uploaded data (MEDEFFECT, FAERS)
# less controlled, the other from more reliable databases (SIDER, OFFSIDE)


faers_grouped = faers.groupby('drug').agg(set).reset_index()
medeffect_grouped = medeffect.groupby('drug').agg(set).reset_index()

df_se_community = pd.merge(faers_grouped, medeffect_grouped, on='drug', how='inner')
df_se_community['union_se'] = df_se_community.apply(lambda row: row['se_x'].union(row['se_y']), axis=1)
df_se_community = df_se_community[['drug', 'union_se']].explode('union_se', ignore_index=True)
df_se_community = df_se_community.rename(columns={'union_se': 'se'})

df_se_controlled = offside.append(sider, ignore_index=True)[['drug', 'se']]

# Now we load the datasets regarding drug target relationships
stitch = pd.read_csv('relationship_analysis_input_files/STITCH_cleaned.input',
                     sep='\t')

dtc = pd.read_csv('relationship_analysis_input_files/DTC_cleaned.input', sep='\t')


stitch = stitch.rename(columns={'compound_name': 'drug',
                                'target_id': 'target',
                                      }
                             )
dtc = dtc.rename(columns={'compound_name': 'drug',
                          'target_id': 'target',
                                      }
                             )

df_target = stitch.append(dtc, ignore_index=True)\
    .groupby(['standard_inchi_key', 'drug'], dropna=False)\
    .agg(set).reset_index()


df_target = df_target[['drug', 'target', 'Database_STITCH', 'Database_DTC', 'SMILES_string']].explode('target').explode('SMILES_string')


df_target = df_target[['drug', 'target', 'SMILES_string']] \
    .groupby('drug') \
    .agg(set) \
    .reset_index().explode('SMILES_string')


df_se_community = meddra_cleaning(df_se_community)
df_se_controlled = meddra_cleaning(df_se_controlled)

interaction_community = target_se_merging(df_se_community, df_target)
interaction_controlled = target_se_merging(df_se_controlled, df_target)
accepted_community = final_adjustements(interaction_community, 'community')
accepted_controlled = final_adjustements(interaction_controlled, 'controlled')

TARDIS_tables(interaction_community, accepted_community, 'community')
TARDIS_tables(interaction_controlled, accepted_controlled, 'controlled')
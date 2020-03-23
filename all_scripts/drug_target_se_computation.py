#!/usr/bin/env python
# coding: utf-8

''' This script is used to combine, merge and extrapolate information from the differetn drug-se, drug-tragets databases. The aim is to retrieve the pairwise relationship between side effects and drug target and statistically validate them through fisher exact test and q-value correction'''

import pandas as pd
import scipy.stats as stats
from multipy.fdr import qvalue


# Load the drug-se databases previusly cleaned

faers = pd.read_csv('relationship_analysis_input_files/faers_pairwise.input', sep='\t', names=['drug', 'se'], dtype=object)
faers = faers.rename(columns={'lookup_value' : 'drug', 'reac_pt_list':'se'}).sort_values('drug')

medeffect = pd.read_csv('relationship_analysis_input_files/MEDEFFECT_DRUG_SE.input', sep='\t', dtype=object, usecols=['DRUGNAME','PT_NAME_ENG'])
medeffect= medeffect.rename(columns={'DRUGNAME' : 'drug', 'PT_NAME_ENG':'se'}).sort_values('drug')

offside = pd.read_csv('relationship_analysis_input_files/OFFSIDE_DRUG_SE.input', sep='\t', dtype=object, usecols=['drug_concept_name','condition_concept_name'])
offside = offside.rename(columns={'drug_concept_name' : 'drug', 'condition_concept_name':'se'}).sort_values('drug')

sider = pd.read_csv('relationship_analysis_input_files/SIDER_DRUG_SE.input', sep='\t', dtype=object, usecols=['DRUGNAME','SIDEEFFECT'])
sider = sider.rename(columns={'DRUGNAME' : 'drug', 'SIDEEFFECT':'se'}).sort_values('drug')


total = faers.append([offside, sider, medeffect], ignore_index=True)
total.sort_values('drug')

total = total.groupby('drug').agg(lambda x: list(set(x.tolist()))).reset_index()

# Now we load the datasets regarding drug target relationships


drugbank = pd.read_csv('relationship_analysis_input_files/Drugbank_relationship.input', sep='\t', usecols = ['Common name', 'UniProt ID'])
drugbank = drugbank.rename(columns={'Common name':'drug','UniProt ID':'target' })

dgidb = pd.read_csv('relationship_analysis_input_files/DGidb_relationships.input', sep='\t', usecols=['drug_claim_primary_name','Uniprot ID'])
dgidb= dgidb.rename(columns={'drug_claim_primary_name':'drug', 'Uniprot ID':'target'})

drugcentral = pd.read_csv('relationship_analysis_input_files/drugcentral_merged.input', sep='\t', usecols=['DRUG_NAME', 'ACCESSION'])
drugcentral = drugcentral.rename(columns={'DRUG_NAME':'drug', 'ACCESSION':'target'})

df_target = drugbank.append([dgidb, drugcentral], ignore_index=True)

df_target_collapsed = df_target.groupby('drug').agg(lambda x: list(set(x.tolist()))).reset_index()


# merge the drug - se database and the drug - target database mapping using the drug name, of course use an inner join to exclude the unmapped matches

interaction = pd.merge(total, df_target_collapsed, how='inner', on='drug')


# Extrapolate the number of drugs that present a particular side effect reconstructing the pairwise relationship between the mapped drugs and sideffect and grouping by the latter.  
# In this way we can compute the number of drugs simply applying the len function.

dr_se_pairwise = interaction[['drug', 'se']].explode('se').groupby('se').agg(lambda x: list(set(x.tolist()))).reset_index()
dr_se_pairwise['se_drug_len'] = dr_se_pairwise['drug'].apply(lambda x: len(x))


# Same thing can be done for defying the number of drugs that present a particular target.

dr_tg_pairwise = interaction[['drug', 'target']].explode('target').groupby('target').agg(lambda x: list(set(x.tolist()))).reset_index()
dr_tg_pairwise['tg_drug_len'] = dr_tg_pairwise['drug'].apply(lambda x: len(x))

# extract from dataset the pairwise interaction between side effect and targets, since both are lists, we need first to explode one list (se) extract the column of interst and then explode the target list.
se_tg_pairwise = interaction.explode('se')[['se', 'target']].explode('target')
new1 = pd.merge(se_tg_pairwise, dr_se_pairwise.rename(columns={'drug':'se_drug'}), how='left', on='se')
new2 = pd.merge(new1, dr_tg_pairwise.rename(columns={'drug':'tg_drug'}), how='left', on='target')


# Define a function to find the number of drug that present at the both time the side effect and the target.  
# Simply for each couple of se and targets take the computed list of drugs for both and check for the shared ones with set.

from pandarallel import pandarallel
pandarallel.initialize()

def overlap(x):
    return len(list(set(x['se_drug']) & set(x['tg_drug'])))

lists = new2[['se_drug', 'tg_drug']] # split the dataframe to reduce memory consumption

values = new2[['se_drug_len', 'tg_drug_len']]

values['overlap_len'] = lists.parallel_apply(overlap, axis=1)

# Now define a function for the p-value computation using the fisher exact test.  
# Given a pair of SE target:  
# *interaction_len* is the total number of drugs found  
# *se_binds* represent the number of drugs that present at the same time the SE AND the target  
# *se_no_binds* is the number of drugs that present only the side effect  
# *no_se_binds* is the number of drugs that only bind the target  
# *no_se_no_binds* is the number of drug that do not present the se and not bind the target  


interaction_len = len(interaction)

def fisher(x, interaction_len=interaction_len):
    se_binds = x['overlap_len']
    se_no_binds = x['se_drug_len'] - x['overlap_len']
    no_se_binds = x['tg_drug_len'] - x['overlap_len']
    no_se_no_binds = interaction_len - se_binds - se_no_binds - no_se_binds
    _, pvalue = stats.fisher_exact([[se_binds, se_no_binds], [no_se_binds, no_se_no_binds]])
    return pvalue


# Select the pair column only in order to simplyfy the computation

lighter = new2[['se', 'target']] # extract only the se and the target

lighter['pvalue'] = values.parallel_apply(fisher, axis=1)


lighter.to_csv('p_value_computed.csv', sep='\t', index=False)


# obtain the q value correction on the p-value computed

_, qvals = qvalue(lighter['pvalue'].tolist())


lighter['qvals'] = qvals


lighter.to_csv('qvalues_all_interactions', sep='\t', index=False)


accepted = lighter.loc[lighter['qvals'] <= 0.05]
accepted.to_csv('accepted_interactions', sep='\t', index=False)

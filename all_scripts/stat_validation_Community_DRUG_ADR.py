import sys
import numpy.random
import pandas as pd
import numpy as np
from numpy.random import normal
from pandarallel import pandarallel
pandarallel.initialize(nb_workers=8, progress_bar=True)


def create_cross_table(pandas_df):
    cross_table = pd.crosstab(pandas_df.iloc[:, 2], pandas_df.iloc[:, 1], margins=True, margins_name='Total_Reports')
    return cross_table


def parallel_log_likelihood_ratio(entries_df, cross_table):
    drugname = entries_df['drugname']
    adverse_event = entries_df['adverse_event']
    report_value = cross_table.loc[adverse_event][drugname]
    if report_value == 0:
        logLR = np.nan
    else:
        total_drug_reports = cross_table[drugname]['Total_Reports']
        total_adverse_event_reports = cross_table.loc[adverse_event]['Total_Reports']
        all_events_reports = cross_table['Total_Reports']['Total_Reports']

        logLR = report_value \
                * (np.log10(report_value) - np.log10(total_adverse_event_reports)) \
                + total_drug_reports \
                * (np.log10(total_drug_reports) - np.log10(all_events_reports - total_adverse_event_reports)) \
                - (report_value + total_drug_reports) \
                * (np.log10(report_value + total_drug_reports) - np.log10(all_events_reports))

    return logLR

# def log_likelihood_ratio(cross_table):
#     result_list = []
#     for adverse_event in list(cross_table.index):
#         for drugname in list(cross_table.columns):
#             report_value = cross_table.loc[adverse_event][drugname]
#             if report_value == 0:
#                 logLR = np.nan
#             else:
#                 total_drug_reports = cross_table[drugname]['Total_Reports']
#                 total_adverse_event_reports = cross_table.loc[adverse_event]['Total_Reports']
#                 all_events_reports = cross_table['Total_Reports']['Total_Reports']
#
#                 logLR = report_value\
#                         * (np.log10(report_value) - np.log10(total_adverse_event_reports))\
#                         + total_drug_reports\
#                         * (np.log10(total_drug_reports) - np.log10(all_events_reports - total_adverse_event_reports))\
#                         - (report_value + total_drug_reports)\
#                         * (np.log10(report_value + total_drug_reports) - np.log10(all_events_reports))
#
#             result_list.append((drugname, adverse_event, logLR))
#     LLR_dataframe = pd.DataFrame(result_list, columns=['DRUGNAME', 'ADVERSE EVENT', 'LLR'])
#
#     return LLR_dataframe


def multinomial_distribution_and_MonteCarlo_sampling(cross_table):

    dist_list = []

    Total_Adverse_Event_reports = cross_table['Total_Reports'].to_list()[:-1]
    Total_Adverse_Event_reports_probabilities = list(map(lambda x:
                                                         x/cross_table['Total_Reports']['Total_Reports'],
                                                         Total_Adverse_Event_reports))
    rng = np.random.default_rng(seed=42)
    for drugname in list(cross_table.columns):
        mult_dis = rng.multinomial(cross_table[drugname]['Total_Reports'],
                                   Total_Adverse_Event_reports_probabilities)
        dist_list.append((drugname, mult_dis))

    dist_df = pd.DataFrame(dist_list, columns=['drugname', 'Distribution'])

    dist_df['mu'] = dist_df['Distribution'].apply(lambda x: np.mean(x, axis=0))
    dist_df['sigma'] = dist_df['Distribution'].apply(lambda x: np.std(x, axis=0))

    dist_df['MC'] = dist_df.apply(lambda x: normal(x['mu'], x['sigma'], 1000), axis=1)

    dist_df['5th_perc'] = dist_df['MC'].apply(lambda x: np.percentile(x, 5))

    return dist_df


if __name__ == '__main__':
    df = pd.read_csv(sys.argv[1], sep='\t')
    database = sys.argv[2]
    crosstable = create_cross_table(df)
    data = [(x, y) for x in list(crosstable.columns)[:-1] for y in list(crosstable.index)[:-1]]  # create every possible relationship between drug and SE
    LLR_dataframe = pd.DataFrame(
        data=data,
        columns=['drugname', 'adverse_event']
    )
    # print(LLR_dataframe)

    LLR_dataframe['logLR'] = LLR_dataframe.parallel_apply(parallel_log_likelihood_ratio,
                                                          cross_table=crosstable,
                                                          axis=1
                                                          )
    # LLR_dataframe.to_csv('LLR_dataframe_temp', sep='\t', index=False)
    LLR_dataframe = LLR_dataframe.dropna()
    distributions_df = multinomial_distribution_and_MonteCarlo_sampling(crosstable)
    # distributions_df.to_csv('MonteCarlo_sampling_distribution', sep='\t', index=False)

    compare_dist_to_LLR = pd.merge(LLR_dataframe,
                                   distributions_df,
                                   on='drugname',
                                   how='left')

    compare_dist_to_LLR['significant'] = np.where(
        compare_dist_to_LLR['logLR'] >= compare_dist_to_LLR['5th_perc'],
        'yes',
        'no'
    )
    positives = compare_dist_to_LLR[compare_dist_to_LLR['significant'] == 'yes']
    positives['Database'] = database
    filtered_positives = positives[['drugname', 'adverse_event', 'logLR', '5th_perc', 'Database']]
    filtered_positives.to_csv('Significant_interaction_' + database + '.input', sep='\t')


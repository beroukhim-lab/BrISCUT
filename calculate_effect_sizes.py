import pandas as pd
import numpy as np
import os
import scipy.stats as stats
import itertools



def calculate_effect_sizes(thedate,ci, infoloc, amptela, amptelb, deltela, deltelb, ampcenta, ampcentb, delcenta, delcentb):
    resultsfolder = '_'.join(['results',thedate,str(ci)+'/'])
    breakpointfolder = 'breakpoint_files_'+thedate+'/'
    tts = [i for i in os.listdir(resultsfolder) if os.path.isdir(i)]
    results = pd.read_csv(resultsfolder+'all_BrISCUT_results.txt',sep='\t')
    results['ID'] = zip(results.type,results.arm, results.direction, results.telcent, results.code)
    all_dedup = results.drop_duplicates('ID')
    #dedup = pd.read_csv(resultsfolder+'PANCAN_BrISCUT_dedup_210319.txt',sep='\t')
    info = pd.read_csv(infoloc, sep='\t',
                       index_col='chromosome_info').transpose().to_dict()

    beta_params = {('amp','tel'): (amptela, amptelb), ('del','tel'): (deltela,deltelb),
              ('amp','cent'): (ampcenta, ampcentb), ('del','cent'): (delcenta, delcentb)}

    #arms = ['13','14','15','21','22'] + [str(a)+'q' for a in range(1,23) if a not in [13,14,15,21,22]] + [str(a)+'p' for a in range(1,23) if a not in [13,14,15,21,22]]

    def coords(arm):
        if arm in ['13', '14', '15', '21', '22']:
            coord = (info[int(arm)]['q_start'], info[int(arm)]['q_end'])
        elif arm.endswith('q'):
            coord = (info[int(arm[:-1])]['q_start'], info[int(arm[:-1])]['q_end'])
        elif arm.endswith('p'):
            coord = (info[int(arm[:-1])]['p_start'], info[int(arm[:-1])]['p_end'])
        else:
            coord = (
            info[int(arm)]['p_start'], info[int(arm)]['p_end'], info[int(arm)]['q_start'], info[int(arm)]['q_end'])
        return coord

    def add_fractions(dedup):
        newstarts = []
        newends = []
        for i in dedup.index:
            arm = dedup.loc[i, 'arm']
            start = dedup.loc[i, 'Peak.Start']
            end = dedup.loc[i, 'Peak.End']
            telcent = dedup.loc[i, 'telcent']
            if (arm.endswith('p') and telcent == 'tel') or (not arm.endswith('p') and telcent == 'cent'):
                newstart = (start - coords(arm)[0] + 1) / (coords(arm)[1] - coords(arm)[0] + 1)
                newend = (end - coords(arm)[0] + 1) / (coords(arm)[1] - coords(arm)[0] + 1)
            else:
                newstart = 1 - ((end - coords(arm)[0] + 1) / (coords(arm)[1] - coords(arm)[0] + 1))
                newend = 1 - ((start - coords(arm)[0] + 1) / (coords(arm)[1] - coords(arm)[0] + 1))
            newstarts.append(newstart)
            newends.append(newend)
        dedup['Peak.Start1'] = newstarts
        dedup['Peak.End1'] = newends
        return dedup

    def get_selection_effect_sizes(all_dedup, thedate):
        all_dedup = add_fractions(all_dedup)
        all_groups = all_dedup.groupby(['type', 'arm', 'direction', 'telcent'])

        massive = []
        for type, arm, direction, telcent in all_groups.groups:
            df = all_groups.get_group((type, arm, direction, telcent))
            df = df.sort('Peak.Start.1')
            df = df.reset_index()
            armlength = coords(arm)[1] - coords(arm)[0] + 1
            starts = df['Peak.Start.1'].tolist()
            ends = df['Peak.End.1'].tolist()
            segments = {i: {'start': 0, 'end': 0, 'empirical': 0, 'expected': 0, 'rel_selection': 1} for i in
                        range(0, len(starts) + 1)}
            for i in range(0, len(starts) + 1):
                if i == 0:
                    segments[i]['start'] = 0
                    segments[i]['end'] = ends[i]
                elif i == len(starts):
                    segments[i]['start'] = starts[i - 1]
                    segments[i]['end'] = 1
                else:
                    segments[i]['start'] = starts[i - 1]
                    segments[i]['end'] = ends[i]

            segments[0]['empirical'] = df.loc[0, 'n_left']
            segments[0]['expected'] = df.loc[0, 'n_left']
            segments[len(starts)]['empirical'] = df.loc[len(starts) - 1, 'n_right']

            for i in range(1, len(segments) - 1):  # 1, 2, 3
                # print i
                if df.loc[i - 1, 'iter'] > df.loc[i, 'iter']:  # if this one is higher iter, like in 2q del, then give
                    # this is the case in rows 0 to 1. correspond to segment 1
                    segments[i]['empirical'] = df.loc[i - 1, 'n_right']
                else:
                    segments[i]['empirical'] = df.loc[i, 'n_left']

            alpha, beta = beta_params[(direction, telcent)]
            pos_total = 0
            neg_total = 0

            for i in range(1, len(segments)):
                segments[i]['expected'] = (stats.beta.cdf(segments[i]['end'], alpha, beta) - stats.beta.cdf(
                    segments[i]['start'], alpha, beta)) * \
                                          (segments[i - 1]['empirical'] / (
                                          stats.beta.cdf(segments[i - 1]['end'], alpha, beta) - stats.beta.cdf(
                                              segments[i - 1]['start'], alpha, beta)))
                segments[i]['rel_selection'] = segments[i]['empirical'] / segments[i]['expected']
            print type, arm, direction, telcent
            print segments

            total_expected = sum([segments[i]['expected'] for i in segments])


            df['first_segment_score'] = segments[0]['empirical'] / (
            (stats.beta.cdf(segments[0]['end'], alpha, beta)) * armlength) * 1000000
            df['rel_selection'] = df.index.map(lambda x: segments[x + 1]['rel_selection'])
            # df['total_expected'] = total_expected
            df['arm_length'] = armlength
            massive.append(df)

        newdf = pd.concat(massive).drop('index', axis=1)
        newdf.to_csv('all_BrISCUT_selection_effect_sizes_'+thedate+'.txt', sep='\t', index=False)
        return newdf

    return get_selection_effect_sizes(all_dedup,thedate)
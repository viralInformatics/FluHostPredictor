# !/usr/bin/python
# -*- coding: utf-8 -*-
# Author:lihuiru

import os
from collections import defaultdict, Counter

import pandas as pd
def parse_hits(hits,low= False):
    query_evalues = defaultdict(float)
    if not hits:
        return 'Unknown'

    processed_hits = []
    hits = hits[0].split("\n")
    for hit in hits:

        prefix = hit.split('_')[0] if not low else hit.split('_')[1].replace("-"," ")
        e_value = float(hit.strip().split(' ')[-1])
        query_evalues[prefix] += e_value
        processed_hits.append(prefix)

    hit_count = Counter(processed_hits)
    most_common_hits = hit_count.most_common()
    highest_count = most_common_hits[0][1]

    tied_hits = [hit for hit, count in most_common_hits if count == highest_count]

    if len(tied_hits) > 1:
        return min(tied_hits, key = lambda hit: query_evalues[hit])
    return most_common_hits[0][0]


def getMostCommonHitProtein(result, dir):

    print("getMostCommonHitProtein ...")
    mostCommonClass = []
    similarVirusIDs = []

    file_path = os.path.join(dir, result)

    df = pd.read_csv(file_path, sep = '\t', header = None)
    df.columns = [
        'query_id', 'subject_id', 'percent_identity', 'alignment_length', 'mismatches',
        'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score'
    ]
    df.loc[:, 'query_id'] = df.loc[:, 'query_id'].apply(lambda x: x.split('_')[0])

    grouped = df.groupby('query_id')

    for queryName, group in grouped:
        top_hits = group.nsmallest(5, 'evalue')
        hits = top_hits['subject_id'].tolist()

        hit_types = [hit.split('_')[0] for hit in hits]
        most_common_hit = max(set(hit_types), key = hit_types.count)

        mostCommonClass.append((queryName, most_common_hit))

        matching_hits = top_hits[top_hits['subject_id'].apply(lambda x: x.split('_')[0] == most_common_hit)]
        if not matching_hits.empty:
            best_hit = matching_hits.nsmallest(1, 'evalue').iloc[0]
            virus_id = best_hit['subject_id'].split('_')[-1]
            similarVirusIDs.append((queryName, virus_id))
        else:
            similarVirusIDs.append((queryName, "Unknown"))

    return mostCommonClass, similarVirusIDs

def getMostCommonHitProteinLowLevelHost(result, dir, Host):
    print("getMostCommonHitProteinLowLevelHost ...")

    Host = {i[0]: i[1] for i in Host}
    mostCommonClass = []

    file_path = os.path.join(dir, result)

    
    df = pd.read_csv(file_path, sep = '\t', header = None)
    df.columns = [
        'query_id', 'subject_id', 'percent_identity', 'alignment_length', 'mismatches',
        'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score'
    ]

    
    grouped = df.groupby('query_id')

    for queryName, group in grouped:
        high_level_host = Host.get(queryName, '')  
        hits = group[group['subject_id'].str.startswith(high_level_host)]

        if not hits.empty:
            
            hit_types = hits['subject_id'].apply(lambda x: x.split('_')[1]).tolist()
            most_common_hit = max(set(hit_types), key = hit_types.count)
            mostCommonClass.append((queryName, most_common_hit))
        else:
            
            mostCommonClass.append((queryName, 'Unknown'))

    return mostCommonClass


if __name__ == '__main__':
    Host,similarVirusIDs = getMostCommonHitProtein("temp/querySeqToHostDB", "./")
    print(Host,similarVirusIDs)
    HostLowLevel = getMostCommonHitProteinLowLevelHost("temp/querySeqToHostDB","./",Host)
    s = getMostCommonHitProtein("temp/querySeqToProteinTypeDB", "./")
    

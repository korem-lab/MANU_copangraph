import sys
import os
import pandas as pd
import Utils.parse_seq as ps

CONTIG_NAME_FIELD = -1

RELEVANT_BUGS = [
    's__Enterococcus_B faecium',
    's__Escherichia coli',
    's__Klebsiella pneumoniae'
]

def contig_to_node_map(gfa):
    with open(gfa) as f:
        segments = [e for e in ps.parse_gfa(f) if e.type == ps.GFATypes.S]
    
    contigs = sorted(list({e.rest.split(':')[CONTIG_NAME_FIELD] for e in segments}))
    nodes = sorted(list({e.nid for e in segments}))
    df = pd.DataFrame(index=nodes, columns=contigs)
    for s in segments:
        contig = s.rest.split('')[CONTIG_NAME_FIELD]
        node = s.nid
        df.loc[node, contig] = 1
    return df

def extract_relevant_classifications(classified_bins):
    filt =  classified_bins.classification.apply(lambda x: any(rbug in x for rbug in RELEVANT_BUGS))
    classified_bins = classified_bins.loc[filt, :]
    relevant_bins = pd.DataFrame(['bin', 'E_faecium', 'E_coli', 'K_pneumoniae'])
    for i, idx in enumerate(classified_bins.index):
        classification = classified_bins.loc[idx, 'classification']
        relevant_bins.loc[i, :] = [classified_bins.loc[idx, 'user_genome']] + [rbug in classification for rbug in RELEVANT_BUGS]
    return relevant_bins
    
def count_state(ncolor, node, persistance_table, persistence):
    if persistence:
        return (ncolor.loc[node, :] & persistance_table).sum()
    else:
        return (ncolor.loc[node, :] & ~persistance_table).sum()
    
def nodes_of_interest(ncolor, contig_node_map, relevant_bins, sample_map):
    
    df = pd.DataFrame(columns=['node', 'num_persistance', 'num_clearance', 'E_faecium', 'E_coli', 'K_pneumoniae'])
    for idx in relevant_bins.index:
        contig_keys = get_contig_keys(relevant_bins.loc[idx, 'bin'], sample_map)
        for c in contig_keys:
            # for every contig from the relevant bin, look up the nodes it's in
            for node in contig_node_map.index[contig_node_map.loc[:, c]]:
                df.loc[len(df), :] = [
                    node, 
                    count_state(ncolor, node, persistence=True), 
                    count_state(ncolor, node, persistence=False),
                    relevant_bins.loc[idx, 'E_feacium'],
                    relevant_bins.loc[idx, 'E_coli'],
                    relevant_bins.loc[idx, 'K_pneumoniae']
                ]
    return df
                
        
    
if __name__ == '__main__':
    
    
    
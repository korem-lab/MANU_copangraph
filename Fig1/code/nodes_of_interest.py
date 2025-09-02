import sys
import os
import glob
import re
from scipy.sparse import lil_matrix, find
from collections import defaultdict
import pandas as pd
import numpy as np
import tqdm
import Utils.parse_seq as ps

CONTIG_NAME_FIELD = 2
SAMPLE_NUMERAL_FIELD = 3

DATA = './data/PanMetagenomeViz/'

RELEVANT_BUGS = [
    's__Enterococcus_B faecium',
    's__Escherichia coli',
    's__Klebsiella pneumoniae'
]
RELEVANT_BUGS_NAME = [e.replace(' ', '_') for e in RELEVANT_BUGS]

def get_contig_sample_key(seg):
    dat = seg.rest.split(':')
    return f'{dat[CONTIG_NAME_FIELD]}-{dat[SAMPLE_NUMERAL_FIELD]}'
    
def contig_to_node_map(gfa):
    print('constructing map...') 
    with open(gfa) as f:
        segments = [e for e in ps.parse_gfa(f) if e.type == ps.GFATypes.S]
    # get the set of unique contig-sample keys 
    contigs = sorted(list({get_contig_sample_key(e) for e in segments}))
    nodes = sorted(list({e.nid for e in segments}))
    contigs_idx = {e:i for i, e in enumerate(contigs)}
    nodes_idx = {e:i for i, e in enumerate(nodes)}
    sp_df = lil_matrix((len(nodes), len(contigs)), dtype=np.uint8)
    
    for s in segments:
        node = s.nid
        contig = get_contig_sample_key(s)
        sp_df[nodes_idx[node], contigs_idx[contig]] = 1
    df = pd.DataFrame.sparse.from_spmatrix(sp_df, index=nodes, columns=contigs)
    return df

def extract_relevant_classifications(gtdbtk, all_species=False):
    global RELEVANT_BUGS_NAME
    if not all_species:
        filt =  gtdbtk.classification.apply(lambda x: any(rbug in x for rbug in RELEVANT_BUGS))
        gtdbtk = gtdbtk.loc[filt, :]
        relevant_bins = pd.DataFrame(columns=['bin'] + RELEVANT_BUGS_NAME)
        for i, idx in enumerate(gtdbtk.index):
            classification = gtdbtk.loc[idx, 'classification']
            relevant_bins.loc[i, :] = [gtdbtk.loc[idx, 'user_genome']] + [rbug in classification for rbug in RELEVANT_BUGS]
        return relevant_bins
    else:
        # filter classifications without species
        print('all species')
        filt = gtdbtk.classification.apply(lambda x: len(re.findall(r'(s__[\w ]+)', x)) != 0)
        gtdbtk = gtdbtk.loc[filt, :]
        species = gtdbtk.classification.apply(lambda x: re.findall(r'(s__[\w ]+)', x)[0].strip().replace(' ', '_'))
        cols = sorted(list(set(species)))
        RELEVANT_BUGS_NAME = cols
        relevant_bins = pd.DataFrame(columns = ['bin'] + cols)
        for i, idx in enumerate(gtdbtk.index):
            s = species.loc[idx]
            relevant_bins.loc[i, :] = [gtdbtk.loc[idx, 'user_genome']] + [s == c for c in cols]
        return relevant_bins
        
def count_state(ncolor, node, persistance_table, persistence):
    if persistence:
        return (ncolor.loc[node, :].values & persistance_table.persistence.values).sum()
    else:
        return (ncolor.loc[node, :].values & ~persistance_table.persistence.values).sum()
    
def nodes_of_interest(ncolor, contig_node_map, relevant_bins, c2b_map, persistence_table):
    
    row, col, _ = find(contig_node_map.sparse.to_coo()) 
    columns=['node', 'contig', 'bin', 'num_persistence', 'num_clearence', 'total_samples', 'is_interesting'] + RELEVANT_BUGS_NAME
    df = pd.DataFrame(np.full((len(row), len(columns)), pd.NA), columns=columns)
    
    for i, (r, c) in tqdm.tqdm(enumerate(zip(row, col)), total=len(row)):
        node = contig_node_map.index[r]
        contig = contig_node_map.columns[c]
        bin = c2b_map.get(contig, None)
        is_interesting = any(relevant_bins.bin == bin)
        nump = count_state(ncolor, int(node), persistence_table, persistence=True)
        numc = count_state(ncolor, int(node), persistence_table, persistence=False)
        if is_interesting:
            bin_index  = relevant_bins.index[relevant_bins.bin == bin][0]
            df.loc[i, :] = [node, contig, pd.NA if bin is None else bin, nump, numc, numc+nump, is_interesting] + [relevant_bins.loc[bin_index, rbug] for rbug in RELEVANT_BUGS_NAME]
        else:
            df.loc[i, :] = [node, contig, pd.NA if bin is None else bin, nump, numc, numc+nump, is_interesting] + [False for rbug in RELEVANT_BUGS_NAME]
            
    return df

def filter_gfa(gfa, noi, min_samples=0, min_degree=0, always_include_interesting=True, attached_to_interesting=False, k_away=10):
   
    segments, links = list(), list() 
    with open(gfa) as f:
        for e in ps.parse_gfa(f):
            if e.type == ps.GFATypes.S:
                segments.append(e)
            else:
                links.append(e)
    node_index = {str(e):i for i, e in enumerate(sorted(list(set(noi.node))))}
    node_name = {i:str(e) for i, e in enumerate(sorted(list(set(noi.node))))}
    adjM = lil_matrix((len(node_index), len(node_index)), dtype=bool)
    
        
        
    print('building adjM')
    for l in tqdm.tqdm(links):
        if l.l_nid == l.r_nid:
            continue
        adjM[node_index[l.l_nid], node_index[l.r_nid]] = True
        adjM[node_index[l.r_nid], node_index[l.l_nid]] = True
    print('nnz: ', adjM.nnz) 
    
    # calculate connectivity k away
    print('get k away')
    adjM_csr = adjM.tocsr()
    adjM_k  = adjM_csr.copy()
    print(k_away)
    for p in tqdm.tqdm(range(2, k_away+1)):
        adjM_k = (adjM_k>0) + (adjM_csr**p > 0)
        
    filtered_nodes = set()
    degree = adjM.sum(axis=1)
    is_interesting = noi.is_interesting
    is_interesting.index = noi.node
    is_interesting = is_interesting[~is_interesting.index.duplicated()]
    print('filtering nodes')
    for idx in tqdm.tqdm(set(noi.index)):
        
        # prioritize interesting nodes
        n = str(noi.loc[idx, 'node'])
        if is_interesting[n] and always_include_interesting:
            filtered_nodes.add(n)
            continue
        
        # otherwise, see if it passes other filters
        if degree[node_index[n]] < min_degree:
            continue
        if noi.loc[idx, 'total_samples'] < min_samples:
            continue
        
        # if it's attached to something interesting
        if attached_to_interesting:
            _, c, _ = find(adjM_k[node_index[n], :])
            for i in c:
                if is_interesting[node_name[i]]:
                    filtered_nodes.add(n)
        else:
            filtered_nodes.add(n)
    
    # use filtered nodes to filter gfa
    filtered_links = list()
    filtered_segs = list()
    print('writing data')
    for l in links:
        if l.l_nid in filtered_nodes and l.r_nid in filtered_nodes:
            filtered_links.append(l)
    for s in segments:
        if s.nid in filtered_nodes:
            filtered_segs.append(s)
    print('num nodes before filter: ',  len(node_index))
    print('num nodes after filter', len(set(e.nid for e in filtered_segs)))
            
    with open(os.path.join(DATA, 'filtered_gfa.gfa'), 'w') as f:
        for e in filtered_segs + filtered_links:
            e.write(f)
            

def contig_to_bin_map(bin_fastas, sample_lookup):
   
    c2b_map = dict()
    for fl in tqdm.tqdm(bin_fastas):
        sample_name = re.findall('(3[0-9]+).bin', fl)[0]
        sample_num = sample_lookup.loc[int(sample_name), 0]
        with open(fl) as fasta:
            for e in ps.parse(fasta, ps.Fasta):
                c2b_map[f'{e.hdr}-{sample_num}'] = os.path.basename(fl.replace('.fa', ''))
    print('c2b_map size: ', len(c2b_map))
    return c2b_map
    
if __name__ == '__main__':
    
    sample_to_numeral = os.path.join(DATA, 'sample_map.csv')
    bin_fastas = glob.glob(os.path.join(DATA, '*.bin.*.fa'))
    gtdbtk = os.path.join(DATA, 'classify/gtdbtk.bac120.summary.tsv')
    ncolor = os.path.join(DATA, 'mdro_pos_02_mo1000_ms100.ncolor.csv')
    gfa = os.path.join(DATA, 'mdro_pos_02_mo1000_ms100.gfa')
    persistence_table = os.path.join(DATA, 'persistence_table.csv')
    noi_fl = os.path.join(DATA, 'nodes_of_interest.pkl')
    
    
    if os.path.exists(noi_fl):
        noi = pd.read_pickle(noi_fl) 
        filter_gfa(gfa, noi, min_samples=8, min_degree=2, always_include_interesting=True, attached_to_interesting=True, k_away=5)
        sys.exit()
    # map sample name to numeral
    sample_to_numeral = pd.read_csv(sample_to_numeral, header=None, index_col=0)
    sample_to_numeral[0] = range(len(sample_to_numeral))
   
    gtdbtk = pd.read_csv(gtdbtk, delimiter='\t')
    relevant_bins = extract_relevant_classifications(gtdbtk, all_species=True)
    
    # make sure the persistence table is ordered according to copanbgraph
    persistence_table = pd.read_csv(persistence_table, index_col=0)
    persistence_table = persistence_table.loc[sample_to_numeral.index, :]
    ncolor = pd.read_csv(ncolor, index_col=0)
    # map every contig to its bin    
    c2b_map  = contig_to_bin_map(bin_fastas, sample_to_numeral)
    # map every contig to node
    #node2contig = contig_to_node_map(gfa)
    node2contig =pd.read_pickle(os.path.join(DATA, 'node2contig.pkl'))
    
    
    # get bins
    noi = nodes_of_interest(ncolor, node2contig, relevant_bins, c2b_map, persistence_table)
    
    # write data
    noi.to_pickle(os.path.join(DATA, 'nodes_of_interest.pkl'))
    filter_gfa(gfa, noi, min_samples=8, min_degree=2, always_include_interesting=True, k_away=5)

    
    
    

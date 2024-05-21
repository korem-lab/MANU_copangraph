import networkx as nx
import Utils.parse_seq as ps
import pandas as pd
import os
import tqdm
DATA = './data/PanMetagenomeViz/'
RELEVANT_BUGS = [
    's__Enterococcus_B faecium',
    's__Escherichia coli',
    's__Klebsiella pneumoniae'
]
if __name__ == '__main__':
    #G = nx.MultiGraph()
    #G.add_nodes_from([(4, {"Color": "red", "Label": "Klieb"}), (5, {"Color": "green", "Label":"Kleib"})])
    #G.add_edge(4, 5)
    #nx.write_gexf(G, 'testing.gexf')
    
    
    noi = pd.read_pickle(os.path.join(DATA, 'nodes_of_interest.pkl'))
    with open(os.path.join(DATA, 'filtered_gfa.gfa')) as f:
        gfa = list(ps.parse_gfa(f))
    
    G = nx.MultiGraph()
    # add all the nodes first
    print('construct nodes')
    for e in tqdm.tqdm(gfa):
        if e.type == ps.GFATypes.S:
            node_tab = noi.loc[noi.node == e.nid, :]
            label = str()
            if node_tab.iloc[0, :]['is_interesting']:
                bugs = node_tab.loc[:, RELEVANT_BUGS].any(axis=0)
                for i, bug in enumerate(bugs):
                    if bug:
                        label += RELEVANT_BUGS[i]
            else:
                label = '-'
            num_p = node_tab.iloc[0,:]['num_persistence']
            num_c = node_tab.iloc[0,:]['num_clearence']
            G.add_node(e.nid, label=label, num_c=num_c, num_p = num_p, total=num_p + num_c)
    print('construct edges')
    for e in tqdm.tqdm(gfa):
        if e.type == ps.GFATypes.L:
            G.add_edge(e.l_nid, e.r_nid)
    nx.write_gexf(G, os.path.join(DATA, 'filtered_gfa.gexf'))
    
    
    
import sys
import os
import pandas as pd
from collections import defaultdict
import Utils.parse_seq as ps

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage: <exe> <cpg.gfa> <cpg_edge_table> <cpg_node_table> <sample_list>')
        sys.exit()
        
    # load node and edge table
    ntable = pd.read_csv(sys.argv[3], index_col=0)
    etable = pd.read_csv(sys.argv[2], index_col=0)

    # load sample names
    sample_names = pd.read_csv(sys.argv[4], header=None)[0]
    
    # load gfa
    sample_seg_map, sample_link_map = defaultdict(list), defaultdict(list)
    with open(sys.argv[1]) as f:
        for e in ps.parse_gfa(f):
            if e.type == ps.GFATypes.S:
                sample_id = e.rest.split(':')[3]
                sample_seg_map[int(sample_id)].append(e)
            else:
                assert int(e.l_nid) <= int(e.r_nid)
                edge_idx = f'{e.l_nid}({"+" if e.l_ori == 1 else "-"}) -> {e.r_nid}({"+" if e.r_ori == 1 else "-"})'
                for sample_id, occ in enumerate(etable.loc[edge_idx, :]):
                    if occ:
                        sample_link_map[sample_id].append(e)
    
    # print each samples gfa
    coasm_sz = len(sample_names)
    print(sample_names)
    for sample_id in sample_seg_map.keys():
        sample_name = sample_names[sample_id].replace('GHERIG_', '')
        print(sample_name)
        with open(os.path.join('./data/GherigGraphQuality/extractions', f'{coasm_sz}_sample_ssasm_sscpg_sslr_{sample_name}.gfa'), 'w') as f:
            for seg in sample_seg_map[sample_id]:
                seg.write(f)
            for link in sample_link_map[sample_id]:
                link.write(f)
                
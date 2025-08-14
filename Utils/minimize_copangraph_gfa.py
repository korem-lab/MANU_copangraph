import os 
import sys
import parse_seq as ps
from collections import defaultdict

if __name__ == '__main__':
    if len(sys.argv) != 2: 
        print('Usage: <exe> <gfa>')
        sys.exit()

    _, gfa = sys.argv

    links = list() 
    by_node = defaultdict(list)
    with open(gfa) as f:
        for e in ps.parse_gfa(f):
            if e.type == ps.GFATypes.L:
                links.append(e)
            else:
                by_node[e.nid].append(e)
        
    with open(gfa.replace('.gfa', '.minimized.gfa'), 'w') as f:
        for k, v in by_node.items():
            v[0].write(f) # write only one sequence per node
        for e in links:
            e.write(f)
        

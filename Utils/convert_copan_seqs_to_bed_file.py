import os
import sys
import parse_seq as ps


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Usage: <exe> <gfa>')
        sys.exit()

    # stream gfa and convert to bed
    gfa = sys.argv[1]
    assert os.path.exists(gfa)
    base_name = os.path.splitext(gfa)[0]
    with open(f'{base_name}.bed', 'w') as b, open(gfa) as g:
        for e in ps.parse_gfa(g):
            if e.type != ps.GFATypes.S:
                continue
            node = e.nid
            _, _, cn, sn, _, lb, rb, ori = e.rest.split(':')
            b.write(f'{sn}_{cn}\t{lb}\t{rb}\t{node}\n')
    print('DONE')
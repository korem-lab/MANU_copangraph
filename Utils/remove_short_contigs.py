import os
import sys
import parse_seq as ps

if __name__ == '__main__':

    if len(sys.argv) != 3:
        print('Usage: <exe> <contigs> <remove_less_than_or_equal_to>')
        sys.exit(-1)
    min_len = int(sys.argv[2])

    with open(sys.argv[1]) as fin, open(sys.argv[1].replace('.fa', '.short.fa'), 'w') as fout:
        for e in ps.parse(fin, ps.Fasta):
            if len(e.seq) <= min_len:
                continue
            e.write(fout)
            


# takes in a metacarvel gfa and a contigs file
# creates two gfa files:
# 1. reformatted gfa file for only contigs included, but adding the sequence
# 2. gfa file with all contigs included.
import sys
import parse_seq as ps

if __name__ == "__main__":
   
    if len(sys.argv) != 3:
        print("usage: <exe> <gfa> <contigs>")
        sys.exit()
    _, gfa_file, contig_file = sys.argv
    contigs = {e.hdr.split()[0]:e.seq for e in ps.parse(open(contig_file), ps.Fasta)}
    with open(gfa_file, 'r') as f:
        data = [l.strip() for l in f]
        data_s = [e.split() for e in data if e[0] != 'L']
        data_l = [e for e in data if e[0] == 'L']

    data_s_seq = ['\t'.join(data_s.pop(0))]
    for e in data_s:
        seg, name, star, rest = e
        assert(name in contigs.keys())
        data_s_seq.append('\t'.join([seg, name, contigs[name], rest]))
    with open(gfa_file.replace('.gfa', '.seq.gfa'), 'w') as f:
        f.write('\n'.join(data_s_seq))
        f.write('\n')
        f.write('\n'.join(data_l))

import os
import sys
import parse_seq as ps

if __name__ == '__main__':

    if len(sys.argv) != 4:
        print('<exe> <long-read-asm> <min contig len> <outdir>')
        sys.exit()
    fa = sys.argv[1]
    lower_bound = int(sys.argv[2])
    out_dir = sys.argv[3]

    with open(fa) as f:
        data = [e for e in ps.parse(f, ps.Fasta) if len(e) >= lower_bound * 10**6] 
    base = fa.split('/')[-2]
    print(base)
    for e in data:
        with open(os.path.join(out_dir, f'{base}_{e.hdr}.fa'), 'w') as f:
            e.write(f)




import parse_seq as ps
import sys
import os


if __name__ == '__main__':

    if len(sys.argv) != 4:
        print('exe : <fasta> <outdir> <sample_name>')
        sys.exit()
    _, fa, outdir, sn = sys.argv

    with open(fa) as fi, open(os.path.join(outdir, f'{sn}.pe_ext.fasta'), 'w') as fo:
        for e in ps.parse(fi, ps.Fasta):
            contig_name, *rest = e.hdr.split(' ')
            e.hdr = f'{sn}:{contig_name}:0:0'
            e.write(fo)

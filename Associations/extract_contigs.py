import parse_seq as ps
import sys
import os

if __name__ == '__main__':
    
    if len(sys.argv) != 3:
        print('usage: <exe> <pe_ext> <out_dir>')
        sys.exit()
    
    _, pe_ext, out_dir = sys.argv
    with open(pe_ext) as f:
        data = [e for e in ps.parse(f, ps.Fasta)]
    processed_contigs = set()
    contigs = list()
    for e in data:
        _, contig, lb, rb = e.hdr.split(':')
        if contig in processed_contigs:
            continue
        else:
            e.hdr = contig
            e.seq = e.seq[int(lb): -int(rb) if int(rb) != 0 else None]
            contigs.append(e)
            processed_contigs.add(contig)
    
    basename = os.path.basename(pe_ext)
    basename = os.path.splitext(basename)[0]
    basename = os.path.splitext(basename)[0]
    print(basename)
    with open(os.path.join(out_dir, f'{basename}.final_contigs.fa'), 'w') as f:
        for e in contigs:
            e.write(f)


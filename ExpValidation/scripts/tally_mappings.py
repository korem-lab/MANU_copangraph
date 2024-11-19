import sys
import os
import pickle
from collections import defaultdict
QNAME=0
RNAME=2

if __name__ == '__main__':

    _, m1, m2, r1, num_reads, out = sys.argv
    s = os.path.basename(r1).replace('PLT', '').replace('_1.fastq.gz', '')
    
    ref_counts = defaultdict(int)
    for sam in [m1, m2]:
        with open(sam) as f:
            prev = None
            for record in f:
                if record.startswith('@'):
                    continue
                record = record.strip().split()

                # if its a new read, empty the ref_seqs set
                if prev != record[QNAME]:
                    prev = record[QNAME]
                    ref_seqs = set()

                if record[RNAME] in ref_seqs:
                    # then we the read has multi-mapped to the reference.
                    # Count it only once, so continue
                    continue
                
                # increment the number of reads mapping to this reference
                ref_counts[record[RNAME]] += 1
                ref_seqs.add(record[RNAME])
    with open(os.path.join(out, f'{s}_ref_counts.pkl'), 'wb') as f:
        pickle.dump((num_reads, ref_counts), f)






   




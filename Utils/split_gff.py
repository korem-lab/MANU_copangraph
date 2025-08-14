import sys
import os
import glob
DATA = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/plt/TaskI/prokka'

if __name__ == '__main__':
    gffs = glob.glob(os.path.join(DATA, '*.gff'))
    for gff in gffs:
        print(gff)
        with open(gff) as f:
            data = f.read()
        try:
            _gff = data[:data.index('##FASTA')]
        except:
            print('token not in ', gff)
            continue
        _fa = data[data.index('##FASTA') + len('##FASTA\n'):]
        with open(gff, 'w') as f:
            f.write(_gff)
        with open(os.path.splitext(gff)[0] + '.fasta', 'w') as f:
            f.write(_fa)

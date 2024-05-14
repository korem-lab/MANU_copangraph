import glob
import re
import os

MANITOU_PATH = '/manitou/pmg/users/ic2465/Projects/MANU_copangraph_2022/data/'
BURG_PATH = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/reads'
if __name__ == '__main__':
    files = glob.glob(os.path.join(MANITOU_PATH + '*.fq.gz'))
    samples = {re.findall('_(3[0-9]+|3[0-9]+-MG)_', f)[0] for f in files}
    for f in files:
        s = re.findall('_(3[0-9]+|3[0-9]+-MG)_', f)[0]
        mate = 'R1' if 'R1' in f else 'R2'
        os.symlink(f, os.path.join(BURG_PATH, f'{s}_{mate}.fastq.gz'))
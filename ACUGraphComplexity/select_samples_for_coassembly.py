import random
import pandas as pd
import re
import glob

REPLICATES = 3
SIZES = [1, 5, 10, 20]

if __name__ == '__main__':
    
    extensions = glob.glob('./data/ACUGraphComplexity/*.pe_ext.fasta')
    coassembly_table = pd.DataFrame(columns=['replicate', 'N', 'sample_name', 'path'])
    for rep in range(REPLICATES):
        for n in SIZES:
            samples = set()
            while len(samples) < n:
                samples.add(random.choice(extensions))
            for s in samples:
                name = re.findall('(3[0-9]+).pe_ext.fasta', s)[0]
                coassembly_table.loc[len(coassembly_table), :] = [rep, n, name, s]
    
    print(coassembly_table.to_string()) 
    coassembly_table.to_csv('./data/ACUGraphComplexity/coassembly_table.csv', index=None)
    
                
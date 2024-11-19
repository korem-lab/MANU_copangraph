import pandas as pd
import sys
import os 
import pickle
import glob


if __name__  == '__main__':

    _, feature, outdir, persistence_tab, feature_occ, seq_len = sys.argv
    df = pd.DataFrame(columns=['feature', 'sequence', 'seq_len', 'library_size', 'sample', 'in_feature', 'persistence', 'n_mapped', 'n_mapped_norm_ls', 'n_mapped_norm_ls_len'])
    feature_occ = pd.read_csv(feature_occ, index_col=0, header=None)
    feature_occ.columns = ['occurence']
    seq_len = pd.read_csv(seq_len, index_col=0, header=None)
    seq_len.columns = ['len']
    persistence_tab = pd.read_csv(persistence_tab, index_col=0)

    # get feature counts
    sequences = [e.replace('"', '') for e in seq_len.index]
    seq_len.index = sequences
    pickle_dumps = glob.glob(os.path.join(outdir, '*_ref_counts.pkl'))
    for file in pickle_dumps:
        lib_sz, counts = pickle.load(open(file, 'rb'))
        lib_sz = int(lib_sz)
        sample = os.path.basename(file).replace('_ref_counts.pkl', '')
        persistence = persistence_tab.loc[int(sample), 'persistence']
        in_feature = feature_occ.loc[int(sample), 'occurence']
        for sequence in sequences:
            slen = seq_len.loc[sequence, 'len']
            n_mapped = counts[sequence]
            n_mapped_norm_ls = n_mapped / float(lib_sz)
            n_mapped_norm_ls_len = n_mapped_norm_ls / float(slen)
            df.loc[len(df), :] = [feature, sequence, slen, lib_sz, sample, in_feature,  persistence, n_mapped, n_mapped_norm_ls, n_mapped_norm_ls_len]
    
    # write df
    df.to_csv(os.path.join(outdir, 'feature_counts.csv'))


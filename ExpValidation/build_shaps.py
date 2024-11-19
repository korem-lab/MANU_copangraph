import pickle
import glob
import os
import sys
import pandas as pd
import shap
import re
import matplotlib.pyplot as plt

COPAN_DIR = '/burg/pmg/users/az2732/copan_MOMS_PI_out'
SHAP_DIR = '/burg/pmg/users/az2732/copan_prediction_out/shap/'
OUT_DIR = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/Predictions/moms-pi'
NUM_FEATURES = 20

if __name__ == '__main__':
    
    shap_pickles = glob.glob(os.path.join(SHAP_DIR, '*.pkl'))
    
    for sp in shap_pickles:
        print(sp)
        
        
        # MDRO filter
        #base = re.findall('nodes_and_edges_(.*)\\.pkl', sp)[0]
        #if '05' not in base:
        #    base = 'mdro_pos_02_' + base
        #else:
        #    base = 'mdro_pos_' + base
        
        # CRC 
        #if 'CRC' not in sp:
        #    continue
        # MOMS-PI 
        if 'MOMS' not in sp:
            continue
        
        sv = pickle.load(open(sp, 'rb'))
        nmap = pd.read_csv(os.path.join(COPAN_DIR, 'MOMS_PI250' + '.ncolor.feature_map.csv'), header=None)
        emap = pd.read_csv(os.path.join(COPAN_DIR, 'MOMS_PI250' + '.ecolor.feature_map.csv'), header=None)
        print(nmap)
        print(emap)
        print(sv)
        
        for model_number in range(10):
            m = sv[model_number]
            plt.clf()
            shap.plots.beeswarm(m, max_display=NUM_FEATURES)
            bn = os.path.basename(sp).replace('.pkl', '')
            plt.tight_layout()
            plt.savefig(os.path.join(OUT_DIR, bn + f'_{model_number}.pdf'))
            # Extract the names of the top features
            tf = pd.DataFrame({'feature': m.feature_names, 'importance': m.abs.mean(0).values})
            tf.loc[:, 'name'] = pd.NA
            for i in tf.index:
                fn = tf.loc[i, 'feature']
                if fn.startswith('n'):
                    fn = int(fn[1:])
                    tf.loc[i, 'name'] = nmap.loc[fn,0]
                else:
                    fn = int(fn[1:])
                    tf.loc[i, 'name'] = emap.loc[fn,0]
            tf.sort_values(by='importance', ascending=False)[:NUM_FEATURES].to_csv(
                os.path.join(OUT_DIR, bn + f'_{model_number}.tf.csv')
            )
                    

import pickle
import glob
import os
import sys
import pandas as pd
import shap
import matplotlib.pyplot as plt

SHAP_DIR = '/manitou/pmg/users/ym2877/copan-mdro-predictions-april2024/shap'
OUT_DIR = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/Predictions/mdro'
NUM_FEATURES = 20

if __name__ == '__main__':
    
    shap_pickles = glob.glob(os.path.join(SHAP_DIR, '*.pkl'))
    
    for sp in shap_pickles:
        sv = pickle.load(open(sp, 'rb'))
        for model_number in range(10):
            m = sv[model_number]
            print(sp, model_number)
            #plt.clf()
            #shap.plots.beeswarm(m, max_display=NUM_FEATURES)
            bn = os.path.basename(sp).replace('.pkl', '')
            #plt.tight_layout()
            #plt.savefig(os.path.join(OUT_DIR, bn + f'_{model_number}.pdf'))
            # Extract the names of the top features
            tf = pd.DataFrame({'feature': m.feature_names, 'importance': m.abs.mean(0).values})
            tf.sort_values(by='importance', ascending=False)[:NUM_FEATURES].feature.to_csv(
                os.path.join(OUT_DIR, bn + f'_{model_number}.tf.csv')
            )

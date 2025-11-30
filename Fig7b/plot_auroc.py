# REQUIRES: 
from scipy.stats import zscore
from sklearn.metrics import roc_curve, roc_auc_score
import pandas as pd 
import dill
import numpy as np
import os
import matplotlib.pyplot as plt

NESTED = True
CLEAN_FIG = True
STUDY = 'ACU_VRE'
FEATURES = ['copan', 'spacegraphcats', 'mOTUs', 'humann', 'clinical']
COLORS = ['tab:blue', 'tab:purple', 'tab:orange', 'tab:green', 'tab:red']  # colors for the plot
LABELS = ['Copangraph', 'spacegraphcats', 'mOTUs3', 'HUMAnN3', 'Clinical']  

def load(pred, ignore=True):
    if isinstance(pred,str):
        print('Loading prediction from %s...'%pred)
        pred = dill.load(open(pred,'rb'),ignore=ignore)
    try:
        nests = len(pred.res[0]['results'])
        robust = len(pred.res[0]['results'][0])
        folds = len(pred.res[0]['results'][0][0])

        print('Detected %s iterations, %s nests (outer folds), %s inner folds, and %s train/test splits!'%(len(pred.res),nests,folds,robust))
    except:
        pass
    return pred

DPI = 1400
FONTSIZE = 10
if CLEAN_FIG:
    TARGET_FILENAME = f"../data/Fig6b/combined_auroc_update_UNLBL.pdf"
else:
    TARGET_FILENAME = f"../data/Fig6b/combined_auroc_update_LBL.pdf"

def find_pickle_filename(feature):
    if feature == 'copan':
        filename = "../data/Fig6b/copan.prediction.pkl"
    elif feature == 'spacegraphcats':
        filename='../data/Fig6b/spacegraphcats.prediction.pkl'
    elif feature == 'mOTUs':
        filename = "../data/Fig6b/motus.prediction.pkl"
    elif feature == 'humann':
        filename = "../data/Fig6b/humann.prediction.pkl"
    elif feature == 'clinical':
        filename = "../data/Fig6b/clinical.prediction.pkl"
    else:
        raise ValueError(f"Couldn't find .pkl file for {feature}")
    return filename

def plot_roc_curve(y_test,y_pred,standardize=False,label=None,color='blue'):
    if standardize and np.std(y_pred).round(7) != 0:
        y_pred = zscore(y_pred)

    auc = roc_auc_score(y_test,y_pred)
    fpr, tpr, _ = roc_curve(y_test,y_pred)

    label = f"auROC={np.round(auc, 2)}" if label is None else f"{label}, auROC={np.round(auc, 2)}"
    plt.plot(fpr, tpr, label=label, color=color, lw=3)

    if CLEAN_FIG:
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])
    else:
        plt.ylabel('True Positive Rate')
        plt.xlabel('False Positive Rate')
        plt.tick_params(axis='both', which='major', labelsize=FONTSIZE)
        plt.legend(loc=4, fontsize=FONTSIZE)
    plt.plot([0, 1], [0, 1], linestyle='--', lw=1, color='black')

def main():
    # Find .pkl prediction files 
    src_files = []
    for feature in FEATURES:
        if NESTED:
            filename = find_pickle_filename(feature) 
            if not os.path.exists(filename):
                raise FileNotFoundError(f"File {filename} does not exist.")
            src_files.append(filename)
            print(f"Found {filename}")
        else:
            raise NotImplementedError("Not implemented for non-nested runs.")
    
    plt.figure(figsize=(4,4))
    for filename,color,label in zip(src_files,COLORS,LABELS):
        print("Loading", filename)
        pred = load(filename)
        if hasattr(pred, 'nested_res'):
            plot_roc_curve(pred.nested_res['y_test'], pred.nested_res['y_pred'], label=label, color=color)
        else:
            raise ValueError(f"Something went wrong... the prediction object does not have nested_res attribute.")

    if not CLEAN_FIG:
        plt.title(f"VRE Prediction task")
    plt.savefig(TARGET_FILENAME, dpi=DPI, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()

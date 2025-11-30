from PredictionPipelineV3 import *
import sys
import pandas as pd 
import os
import matplotlib.pyplot as plt

NESTED = True
CLEAN_FIG = False
STUDY = 'MOMPIST'
FEATURES = ['copan', 'mOTUs', 'humann']
COLORS = ['tab:blue', 'tab:orange', 'tab:green']  # colors for the plot
LABELS = ['Copangraph', 'mOTUs3', 'HUMAnN3']  
DPI = 1400
FONTSIZE = 10
def find_pickle_filename(feature):
    if feature == 'copan':
        filename = "../data/Fig6gh/momspi/graph_sd02_mh1000_ms100_abund_deep.prediction.pkl"
    elif feature == 'mOTUs':
        filename = "../data/Fig6gh/momspi/motus.prediction.pkl"
    elif feature == 'humann':
        filename = "../data/Fig6gh/momspi/humann.prediction.pkl"
    else:
        raise ValueError(f"Couldn't find .pkl file for {feature}")
    return filename

def plot_roc_curve(y_test,y_pred,standardize=False,label=None,color='blue',is_clean=False):
    if standardize and np.std(y_pred).round(7) != 0:
        y_pred = zscore(y_pred)

    auc = roc_auc_score(y_test,y_pred)
    fpr, tpr, _ = roc_curve(y_test,y_pred)

    label = f"auROC={np.round(auc, 2)}" if label is None else f"{label}, auROC={np.round(auc, 2)}"
    plt.plot(fpr, tpr, label=label, color=color, lw=3)

    if is_clean:
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
    CLEAN_FIG = True if sys.argv[1] == 'clean' else False
    if CLEAN_FIG:
        TARGET_FILENAME = f"../data/Fig6gh/Fig6h_UNLBL.pdf"
    else:
        TARGET_FILENAME = f"../data/Fig6gh/FIg6h_LBL.pdf"
    
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
        pred = Prediction.load(filename)
        if hasattr(pred, 'nested_res'):
            plot_roc_curve(pred.nested_res['y_test'], pred.nested_res['y_pred'], label=label, color=color,is_clean=CLEAN_FIG)
        else:
            raise ValueError(f"Something went wrong... the prediction object does not have nested_res attribute.")

    if not CLEAN_FIG:
        plt.title(f"MOMSPI Prediction task")
    plt.savefig(TARGET_FILENAME, dpi=DPI, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()

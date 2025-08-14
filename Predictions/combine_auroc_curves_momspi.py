from PredictionPipelineV3 import *
import pandas as pd 
import os
import matplotlib.pyplot as plt

NESTED = True
CLEAN_FIG = False
STUDY = 'MOMS-PI250'  # CRC-CN, MOMS-PI250, ACU
FEATURES = ['copan', 'mOTUs', 'humann-orgspc', 'sgc']  # copan or copan_0.05 (ACU only); humann or humann-orgspc
COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']  # colors for the plot
LABELS = ['Copangraph', 'mOTUs 3', 'HUMAnN 3', 'Spacegraphcats']  # what is shown on the plot. $\mathbf{*}$ for significance (humann, except in MOMS-PI250)
DPI = 1400
FONTSIZE = 16
if CLEAN_FIG:
    TARGET_FILENAME = f"/burg/pmg/users/az2732/Copan_Runs/plots/{STUDY}_new_combinedROCs_dpi{str(DPI)}_CLEAN.pdf"
else:
    TARGET_FILENAME = f"/burg/pmg/users/az2732/Copan_Runs/plots/{STUDY}_new_combinedROCs_dpi{str(DPI)}.pdf"
ACU_FEATURE_MAP = {
    'copan': 'NODES_AND_EDGES_mo1000_ms75',
    'copan_0.05': 'NODES_AND_EDGES_05_mo1000_ms100',
    'mOTUs': 'MOTUS_CLEANV2_02',
    'humann-orgspc': 'HUMANN_GxO_CLEANV2_02',
    'humann-g': 'HUMANN_G_CLEANV2_02',
    'sgc': 'SGC_CLEANV2_02',
    'clinical': 'CLINICAL'
}

def find_pickle_filename(feature):
    if STUDY == 'ACU':
        SRC_DIR = "/manitou/pmg/users/ym2877/copan-mdro-predictions-april2024/preds"
        feature_coded = ACU_FEATURE_MAP[feature]
        filename = f"{SRC_DIR}/mdro_pos_{feature_coded}.pkl"
    elif STUDY == 'MOMS-PI250' and feature == 'copan':
        filename = "/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/moms-pi/results/occ/momspi-random_it10_sd005.prediction.pkl"
    else:
        SRC_DIR = "/burg/pmg/users/az2732/copan_prediction/copan_prediction_out"
        filename = f"{SRC_DIR}/{STUDY}_{feature}_nest.pkl"
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
            filename = find_pickle_filename(feature) # f"{SRC_DIR}/{STUDY}_{feature}_nest.pkl"
            if not os.path.exists(filename):
                raise FileNotFoundError(f"File {filename} does not exist.")
            src_files.append(filename)
            print(f"Found {filename}")
        else:
            raise NotImplementedError("Not implemented for non-nested runs.")
    
    plt.figure(figsize=(8, 8))
    for filename,color,label in zip(src_files,COLORS,LABELS):
        print("Loading", filename)
        pred = Prediction.load(filename)
        if hasattr(pred, 'nested_res'):
            plot_roc_curve(pred.nested_res['y_test'], pred.nested_res['y_pred'], label=label, color=color)
        # spacegraphcat's output on MOMS-PI is not 
        else:
            if STUDY == 'MOMS-PI250' and 'Spacegraphcats' in label:
                nested_res = pd.read_csv(f"/burg/pmg/users/az2732/Copan_Runs/results/SGC_MOMS-PI/SGC_MOMSPI_nested_res.csv")
                plot_roc_curve(nested_res.y_test, nested_res.y_pred, label=label, color=color)
            else:
                raise ValueError(f"Something went wrong... the prediction object does not have nested_res attribute.")

    if not CLEAN_FIG:
        plt.title(f"MOMS-PI Prediction task")
    plt.savefig(TARGET_FILENAME, dpi=DPI, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()

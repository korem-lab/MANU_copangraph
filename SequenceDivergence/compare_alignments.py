import random

import evaluation.utils.parse_seq as ps
from calign import aligner
import sys
import json
import tqdm
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def align(s1, s2):
    aln_obj = aligner(s1, s2)[0]
    minlen = min(len(s1), len(s2))
    if minlen == 0:
        return 1
    score = 1 - aln_obj.n_mismatches / min(len(s1), len(s2))
    return score


def build_copan_comparison(copan_dat, gt):
    gt = {i: score for i, score in enumerate(gt)}
    copan = [1 - e['sd'] for e in copan_dat]
    gt = [gt[e['s_cid'] - 1] for e in copan_dat]
    return gt, copan


def filter_pass(e, sd, ori):
    return (e['sd'] <= sd) and (e['ori'] == ori)


if __name__ == '__main__':

    SD_THRESH = 0.10

    fa_fl = sys.argv[1]
    with open(fa_fl) as f:
        fastas = [fa for fa in ps.parse(f, ps.Fasta)]

    copan_dat = list()
    print('extracting copan alignments')
    total = 0
    with open(sys.argv[2]) as f:
        for i, e in enumerate(f):
            total += 1
            # output has a header
            e = json.loads(e.strip())
            if filter_pass(e, sd=SD_THRESH, ori='forward'):
                copan_dat.append(e)

    aligner_scores = []
    print(f'remaining after filtering: {len(copan_dat)/total}')
    print('brentp aligning regions')
    for aln in tqdm.tqdm(copan_dat):
        try:
            q = fastas[aln['q_cid']].seq[aln['q_beg']: aln['q_end']]
            s = fastas[aln['s_cid']].seq[aln['s_beg']: aln['s_end']]
            aligner_scores.append(align(q, s))
        except IndexError:
            print(aln, 'couldnt find index')

    # rgx = '_([01]\.[0-9]+$)'
    # gt_clipped = [float(re.findall(rgx, fa.hdr)[0]) for fa in fastas[:5]]
    # gt = [float(re.findall(rgx, fa.hdr)[0]) for fa in fastas]

    # x, y = np.array(gt_clipped), np.array(aligner_scores)
    # sns.scatterplot(x=x, y=y)
    # plt.title('brentp seqid vs gt')
    # plt.savefig('brentp-gt.pdf')
    # plt.show()
    # plt.clf()

    # x, y = build_copan_comparison(copan_dat, gt)
    # sns.scatterplot(x=x, y=y)
    # plt.title('Max jump = l 1000')
    # plt.savefig('copan-gt-mj1000.pdf')
    # plt.show()
    # plt.clf()

    ax = sns.scatterplot(x=np.array(aligner_scores), y=1 - np.array([e['sd'] for e in copan_dat]), edgecolor=None, s=10, alpha=0.5)
    ax.set_ylim(1 - SD_THRESH * 5, 1)
    plt.axhline(y=1 - SD_THRESH, color='red')
    plt.axvline(x=1 - SD_THRESH, color='red')
    plt.axvline(x=np.mean(np.array(aligner_scores)), color='blue')
    plt.axvline(x=np.median(np.array(aligner_scores)), color='green'
    plt.title(f'copan sd={SD_THRESH} mj=1000 vs brentp')
    plt.savefig(f'copan-brentp-mj1000-sd{SD_THRESH}.pdf')
    plt.show()
    plt.clf()

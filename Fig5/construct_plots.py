
import pandas as pd
import os 
from matplotlib import pyplot as plt
import seaborn as sns

def plot_unlabelled_version(ax, name, tight_layout=True):
    name = os.path.splitext(name)[0]
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    frame1 = plt.gca()
    frame1.legend().set_visible(False)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(f'{name}_UNLBLD.pdf', dpi=1400, bbox_inches='tight')
    plt.savefig(f'{name}_UNLBLD.png', dpi=900, bbox_inches='tight')
    plt.clf()


if __name__ == '__main__':

     # megahit mean genomes per node
    plot_dat = pd.read_csv('../data/Fig5/resource_dat.csv', index_col=0)
    plt.figure(figsize=(4,4))
    hue_order = ['copangraph','megahit', 'metaspades']
    pal = ['#4c72b0', '#c44e52', '#55a868']
    ax=sns.lineplot(x=plot_dat.coasm_sz, y=plot_dat['time(s)']/3600, hue=plot_dat.tool, hue_order=hue_order, palette=pal, errorbar='sd')
    ax=sns.scatterplot(x=plot_dat.coasm_sz, y=plot_dat['time(s)']/3600, hue=plot_dat.tool, hue_order=hue_order, palette=pal, ax=ax)
    y_max = (plot_dat['time(s)']/3600).max()
    ax.set_ybound(0, y_max * 1.2)
    ax.set_ylabel('time (hrs)')
    name = '../data/Fig5/time_resource_plot'
    plt.savefig(name + '.pdf', dpi=1400, bbox_inches='tight')
    plot_unlabelled_version(ax, name)

    ax=sns.lineplot(x=plot_dat.coasm_sz, y=plot_dat['mem(kb)']/2**20, hue=plot_dat.tool, hue_order=hue_order, palette=pal, errorbar='sd')
    ax=sns.scatterplot(x=plot_dat.coasm_sz, y=plot_dat['mem(kb)']/2**20, hue=plot_dat.tool, hue_order=hue_order, palette=pal, ax=ax)
    y_max = (plot_dat['mem(kb)']/2**20).max()
    ax.set_ybound(0, y_max * 1.2)
    ax.set_ylabel('memory (GB)')
    name = '../data/Fig5/mem_resource_plot'
    plt.savefig(name + '.pdf', dpi=1400, bbox_inches='tight')
    plot_unlabelled_version(ax, name)


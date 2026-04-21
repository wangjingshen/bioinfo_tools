import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import wilcoxon, mannwhitneyu, ranksums
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hcluster
from scipy.spatial.distance import squareform

def cohens_d(x, y):
    pooled_std = np.sqrt(((len(x)-1) * np.var(x, ddof=1) 
                          + (len(y)-1) * np.var(y, ddof=1)) / 
                             (len(x) + len(y) - 2))
    return (np.mean(x) - np.mean(y)) / pooled_std
    

def wilcoxon_test(consistencies_matrix, group_A_cells, group_B_cells):
    """
        Performs an unpaired wilcoxon test for each reaction between group_A and group_B
    """
    group_A = consistencies_matrix.loc[:,group_A_cells]
    group_B = consistencies_matrix.loc[:,group_B_cells]
    results = pd.DataFrame(index = consistencies_matrix.index, columns = ['wilcox_stat', 'wilcox_pval', 'cohens_d'], dtype='float64')
    for rxn in consistencies_matrix.index:
        A, B = group_A.loc[rxn].to_numpy().ravel(), group_B.loc[rxn].to_numpy().ravel()
        stat, pval = mannwhitneyu(A, B, alternative='two-sided')
        c_d = cohens_d(A, B)
        results.loc[rxn, ['wilcox_stat', 'wilcox_pval', 'cohens_d']] = stat, pval, c_d
    results['adjusted_pval'] = np.array(multipletests(results['wilcox_pval'], method='fdr_bh')[1], dtype='float64')
    return results

def get_reaction_consistencies(compass_reaction_penalties, min_range=1e-3):
    """
    Converts the raw penalties outputs of compass into scores per reactions where higher numbers indicate more activity
    """
    df = -np.log(compass_reaction_penalties + 1)
    df = df[df.max(axis=1) - df.min(axis=1) >= min_range]
    df = df - df.min().min()
    return df

def get_group_barcode(group_info):
    '''
    Cell    Group
    Barcode1        G1
    Barcode2        G1
    Barcode3        G2
    '''
    group_info = pd.read_csv(group_info,sep="\t", header=0)
    g1_name = group_info.iloc[:,1].unique()[0]
    g2_name = group_info.iloc[:,1].unique()[1]
    g1_barcode_list = group_info.iloc[:,0][ group_info.iloc[:,1] == g1_name]
    g2_barcode_list = group_info.iloc[:,0][ group_info.iloc[:,1] == g2_name]
    return g1_name, g2_name, g1_barcode_list, g2_barcode_list

def get_wilcox_result(reaction_penalties_data, metadata, g1_cells, g2_cells):
    reaction_consistencies = get_reaction_consistencies(reaction_penalties_data)
    reaction_metadata = pd.read_csv(metadata, index_col = 0)

    wilcox_results = wilcoxon_test(reaction_consistencies, g1_cells, g2_cells)
    wilcox_results['metadata_r_id'] = ""
    for r in wilcox_results.index:
        if r in reaction_metadata.index:
            wilcox_results.loc[r, 'metadata_r_id'] = r
        elif r[:-4] in reaction_metadata.index:
            wilcox_results.loc[r, 'metadata_r_id'] = r[:-4]
        else:
            print("Should not occur")

    W = wilcox_results.merge(reaction_metadata, how='left', left_on='metadata_r_id', right_index=True, validate='m:1')
    W = W[W['confidence'].isin([0,4])]
    W = W[~W['EC_number'].isna()]
    W.loc[(W['formula'].map(lambda x: '[m]' not in x)) & (W['subsystem'] == "Citric acid cycle"), 'subsystem'] = 'Other'
    return(W)

def plot_differential_scores(plot_data, group1, group2, title, plot_type, outdir, name):
    plt.figure(figsize=(7,6))
    axs = plt.gca()
    axs.scatter(plot_data['cohens_d'], -np.log10(plot_data['adjusted_pval']), c='#695D73')
    axs.set_xlabel("Cohen's d", fontsize=16)
    axs.set_ylabel("-log10 (Wilcoxon-adjusted p)", fontsize=16)
    axs.set_xlim(-2.2, 2.2)
    axs.axvline(0, dashes=(3,3), c='black')
    axs.axhline(1, dashes=(3,3), c='black')
    axs.set_title(title, fontdict={'fontsize':20})
    axs.annotate('', xy=(0.5, -0.06), xycoords='axes fraction', xytext=(0, -0.06), arrowprops=dict(arrowstyle="<-", color='#348C73', linewidth=4))
    axs.annotate(group1, xy=(0.7, -0.12), xycoords='axes fraction', fontsize=16)
    axs.annotate('', xy=(0.5, -0.06), xycoords='axes fraction', xytext=(1, -0.06), arrowprops=dict(arrowstyle="<-", color='#E92E87', linewidth=4))
    axs.annotate(group2, xy=(0.2, -0.12), xycoords='axes fraction', fontsize=16)
    plt.savefig(f'{outdir}/{name}.pdf', dpi=200)
    
    if(plot_type == "anno"):
        for r in plot_data.index:
            x = plot_data.loc[r, 'cohens_d']
            y = -np.log10(plot_data.loc[r, 'adjusted_pval'])
            if(plot_data.loc[r, 'cohens_d'] <= 0):
                offset = (-150, 0)
                axs.annotate(r + plot_data.loc[r, 'reaction_name'], (x,y), xytext = offset, textcoords='offset pixels', arrowprops={'arrowstyle': "-"})
            else:
                offset = (10,-20)
                axs.annotate(r + plot_data.loc[r, 'reaction_name'], (x,y), xytext = offset, textcoords='offset pixels', arrowprops={'arrowstyle': "-"})
            plt.savefig(f'{outdir}/{name}_anno.pdf', dpi=200)

def plot(wilcox_res, plot_type, subsystem, group1, group2, outdir, name):
    plot_data = wilcox_res[wilcox_res['subsystem'] == subsystem]
    plot_differential_scores(plot_data = plot_data, plot_type = plot_type, group1 = group1, group2= group2, title=subsystem, outdir = outdir, name = name)

def run(reaction_penalties, metadata, group_info, subsystem, plot_type, outdir, name):
    reaction_penalties = pd.read_csv(reaction_penalties, sep="\t", index_col = 0)
    g1_name, g2_name, g1_barcode_list, g2_barcode_list = get_group_barcode(group_info)
    print(f'wilcox test: {g1_name} vs {g2_name}')
    wilcox_res = get_wilcox_result(reaction_penalties, metadata, g1_barcode_list, g2_barcode_list)
    wilcox_res.to_csv(f'{outdir}/wilcox_res.tsv',sep="\t")
    plot(wilcox_res, plot_type = plot_type, subsystem = subsystem, group1 = g1_name, group2 = g2_name, outdir = outdir, name = name)
    print("Compass plot done.")

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--reaction_penalties', help = 'reactions.tsv from the result of compass')
    parsers.add_argument('--metadata', help = 'reaction metadata provided by the author')
    parsers.add_argument('--group_info', help = 'a tsv contains group info. The first column corresponding to barcodes, The second column corresponding to group')
    parsers.add_argument('--subsystem', help = 'subsystem used to plot, such as Glycolysis/gluconeogenesis. More subsystem can be found in the reaction metadata provided by the author')
    parsers.add_argument('--plot_type', help = 'plot type, anno or no(t plot anno). Default: no', default = 'no')
    parsers.add_argument('--outdir', help = 'outdir', default = './')
    parsers.add_argument('--name', help = 'name, just like subsystem, special symbols need to be removed')
    args = parsers.parse_args()

    if not os.path.exists(f"{args.outdir}/"):
        os.system(f"mkdir -p {args.outdir}/")
    
    run(args.reaction_penalties, args.metadata, args.group_info, args.subsystem, args.plot_type, args.outdir, args.name) 

if __name__ == '__main__':
    main()
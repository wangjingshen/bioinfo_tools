# need fix
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from compass_analysis import cohens_d, wilcoxon_test, get_reaction_consistencies, get_metareactions

def get_reaction_consistencies(compass_reaction_penalties, min_range=1e-3):
    """
    Converts the raw penalties outputs of compass into scores per reactions where higher numbers indicate more activity
    """
    df = -np.log(compass_reaction_penalties + 1)
    df = df[df.max(axis=1) - df.min(axis=1) >= min_range]
    df = df - df.min().min()
    return df

def get_group_barcode(group_metadata):
    '''
    Cell    Group
    Barcode1        G1
    Barcode2        G1
    Barcode3        G2
    Barcode4        G2
    Barcode5        G1
    '''
    group_metadata = pd.read_csv(group_metadata,sep="\t", header=0)
    g1_name = group_metadata.iloc[:,1].unique()[0]
    g2_name = group_metadata.iloc[:,1].unique()[1]
    g1_barcode_list = group_metadata.iloc[:,0][ group_metadata.iloc[:,1] == g1_name]
    g2_barcode_list = group_metadata.iloc[:,0][ group_metadata.iloc[:,1] == g2_name]
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

    W = wilcox_results.merge(reaction_metadata, how='left',
                             left_on='metadata_r_id', right_index=True, validate='m:1')
    W = W[W['confidence'].isin([0,4])]
    W = W[~W['EC_number'].isna()]
    W.loc[(W['formula'].map(lambda x: '[m]' not in x)) & (W['subsystem'] == "Citric acid cycle"), 'subsystem'] = 'Other'
    return(W)

def plot_differential_scores(plot_data, group1, group2, title, plot_type, outdir):
    plt.figure(figsize=(7,6))
    axs = plt.gca()
    axs.scatter(plot_data['cohens_d'], -np.log10(plot_data['adjusted_pval']), c='#695D73')
    axs.set_xlabel("Cohen's d", fontsize=16)
    axs.set_ylabel("-log10 (Wilcoxon-adjusted p)", fontsize=16)
    #Everything after this should be tweaked depending on your application
    axs.set_xlim(-2.2, 2.2)
    axs.axvline(0, dashes=(3,3), c='black')
    axs.axhline(1, dashes=(3,3), c='black')
    axs.set_title(title, fontdict={'fontsize':20})
    axs.annotate('', xy=(0.5, -0.06),
                xycoords='axes fraction',
                xytext=(0, -0.06), 
                arrowprops=dict(arrowstyle="<-", color='#348C73', linewidth=4))
    axs.annotate(group1,
                xy=(0.7, -0.12),
                xycoords='axes fraction',
                fontsize=16)
    axs.annotate('', xy=(0.5, -0.06),
                xycoords='axes fraction',
                xytext=(1, -0.06),
                arrowprops=dict(arrowstyle="<-", color='#E92E87', linewidth=4))
    axs.annotate(group2,
                xy=(0.2, -0.12),
                xycoords='axes fraction',
                fontsize=16)
    plt.savefig(f'{outdir}/compass_plot.pdf', dpi=200)
    
    if(plot_type == "anno"):
        for r in plot_data.index:
            x = plot_data.loc[r, 'cohens_d']
            y = -np.log10(plot_data.loc[r, 'adjusted_pval'])
            if(plot_data.loc[r, 'cohens_d'] <= 0):
                offset = (-150, 0)
                axs.annotate(r + plot_data.loc[r, 'reaction_name'], 
                        (x,y),
                        xytext = offset,
                        textcoords='offset pixels',
                        arrowprops={'arrowstyle': "-"})
            else:
                offset = (10,-20)
                axs.annotate(r + plot_data.loc[r, 'reaction_name'], 
                        (x,y),
                        xytext = offset,
                        textcoords='offset pixels',
                        arrowprops={'arrowstyle': "-"})
            plt.savefig(f'{outdir}/compass_plot_anno.pdf', dpi=200)


def plot(wilcox_res, plot_type, subsystem, group1, group2, outdir):
    plot_data = wilcox_res[wilcox_res['subsystem'] == subsystem]
    plot_differential_scores(
        plot_data = plot_data,                     
        plot_type = plot_type,
        group1 = group1,
        group2= group2,
        title=subsystem,
        outdir = outdir)

def run(reaction_penalties, metadata, group_metadata, subsystem, outdir):
    reaction_penalties = pd.read_csv(reaction_penalties, sep="\t", index_col = 0)
    g1_name, g2_name, g1_barcode_list, g2_barcode_list = get_group_barcode(group_metadata)
    wilcox_res = get_wilcox_result(reaction_penalties, metadata, g1_barcode_list, g2_barcode_list)
    wilcox_res.to_csv(f'{outdir}/wilcox_res.tsv',sep="\t")
    plot(wilcox_res, plot_type = "no", subsystem = subsystem, group1 = g1_name, group2 = g2_name, outdir = outdir)
    #plot(wilcox_res, plot_type = "anno", subsystem = subsystem, group1 = g1_name, group2 = g2_name, outdir = outdir)  # add anno, not recommended


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--reaction_penalties', help = 'reaction_penalties')
    parsers.add_argument('--metadata', help = 'metadata')
    parsers.add_argument('--group_metadata', help = 'group_metadata')
    parsers.add_argument('--subsystem', help = 'subsystem')
    parsers.add_argument('--outdir', help = 'outdir', default = './')
    args = parsers.parse_args()

    if not os.path.exists(f"{args.outdir}/"):
        os.system(f"mkdir -p {args.outdir}/")
    
    run(args.reaction_penalties, args.metadata, args.group_metadata, args.subsystem, args.outdir) 

if __name__ == '__main__':
    main()
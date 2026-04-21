import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from scipy.stats import wilcoxon, mannwhitneyu, ranksums
from statsmodels.stats.multitest import multipletests

class ScFEA_wilcoxon_test():
    """
    run wilcoxon test like compass.Currently only two cells are supported for comparison.
    Multiple(>=3) groups can be split into two groups and run the process.

    Test path: /SGRNJ06/randd/USER/wangjingshen/project/scFEA1.1.2/test/test_neuroblastoma_data/run_wilcoxon_test.sh
        
        1)flux: flux from scFEA
        2)metadata:
            group
        C1      A
        C2      B

        3)gname:  A,B
        4)reaction: reaction.Name in scFEA

    """
    def __init__(self, args):
        self.flux = pd.read_csv(args.flux, sep=",", index_col = 0)
        self.metadata = pd.read_csv(args.metadata,sep="\t")
        self.gname = str.split(args.gname,',')
        self.supermodule = args.supermodule
        self.outdir = args.outdir
        if not os.path.exists(self.outdir):
            os.system(f'mkdir -p {self.outdir}')
    
    def cohens_d(self, x, y):
        pooled_std = np.sqrt(((len(x)-1) * np.var(x, ddof=1) + (len(y)-1) * np.var(y, ddof=1)) / (len(x) + len(y) - 2))
        return (np.mean(x) - np.mean(y)) / pooled_std


    def wilcoxon_test(self, consistencies_matrix, group_A_cells, group_B_cells):
        """
                Performs an unpaired wilcoxon test (or mann-whitney U test) for each reaction between group_A and group_B
        """
        #per reaction/meta-reaction, conduct wilcoxon test between group_A and group_B
        group_A = consistencies_matrix.loc[:,group_A_cells]
        group_B = consistencies_matrix.loc[:,group_B_cells]
        results = pd.DataFrame(index = consistencies_matrix.index, columns = ['wilcox_stat', 'wilcox_pval', 'cohens_d'], dtype='float64')
        for rxn in consistencies_matrix.index:
                A, B = group_A.loc[rxn].to_numpy().ravel(), group_B.loc[rxn].to_numpy().ravel()
                #sometimes there's a solitary value, and we don't want to test then
                if len(np.unique(A)) == 1 and len(np.unique(B)) == 1:
                        if np.unique(A) == np.unique(B):
                                #we've got no data. set p-value to 1 and skip!
                                #(p-value needs to be 1 so multipletests doesn't cry)
                                results.loc[rxn, ['wilcox_pval']] = 1
                                continue
                stat, pval = mannwhitneyu(A, B, alternative='two-sided')
                c_d = self.cohens_d(A, B)
                results.loc[rxn, ['wilcox_stat', 'wilcox_pval', 'cohens_d']] = stat, pval, c_d
        results['adjusted_pval'] = np.array(multipletests(results['wilcox_pval'], method='fdr_bh')[1], dtype='float64')
        return results
    
    def plot_differential_scores(self, data, anno1, anno2, title, c):
        plt.figure(figsize=(7,6))
        axs = plt.gca()
        axs.scatter(data['cohens_d'], -np.log10(data['adjusted_pval']), c=c)
        axs.set_xlabel("Cohen's d", fontsize=16)
        axs.set_ylabel("-log10 (Wilcoxon-adjusted p)", fontsize=16)
        axs.set_xlim(-2.2, 2.2)
        axs.axvline(0, dashes=(3,3), c='black')
        axs.axhline(1, dashes=(3,3), c='black')
        axs.set_title(title, fontdict={'fontsize':20})
        axs.annotate('', xy=(0.5, -0.06), xycoords='axes fraction', xytext=(0, -0.06),arrowprops=dict(arrowstyle="<-", color='#348C73', linewidth=4))
        axs.annotate(anno1, xy=(0.7, -0.12), xycoords='axes fraction', fontsize=16)
        axs.annotate('', xy=(0.5, -0.06), xycoords='axes fraction', xytext=(1, -0.06),arrowprops=dict(arrowstyle="<-", color='#E92E87', linewidth=4))
        axs.annotate(anno2, xy=(0.2, -0.12), xycoords='axes fraction', fontsize=16)
        plt.savefig(self.outdir + "/" + self.supermodule + "_wilcoxon_test.pdf")
    
    def analysis(self):
        g1_cells = self.metadata.index[self.metadata['group'] == self.gname[0]]
        g2_cells = self.metadata.index[self.metadata['group'] == self.gname[1]]
        # run wilcox
        wilcox_results = self.wilcoxon_test(self.flux.T, g1_cells, g2_cells)
        module_info = pd.read_csv("/SGRNJ03/randd/user/wangjingshen/bio_soft/scFEA_v1.1.2/data/module_info_add_supermodule.tsv",sep="\t", index_col=0)
        wilcox_results.index = module_info.loc[wilcox_results.index,'Name']

        module_plot = module_info[ module_info.Supermodule_name== self.supermodule]
        plot_sub = wilcox_results.loc[module_plot.Name,]
        self.plot_differential_scores(plot_sub, self.gname[0], self.gname[1], self.supermodule, "#695D73")

    def run(self):
        self.analysis()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--flux', help='flux', required = True)
    parsers.add_argument('--metadata', help="metadata", required = True)
    parsers.add_argument('--gname', help="gname", required = True)
    parsers.add_argument('--supermodule', help="supermodule", required = True)
    parsers.add_argument('--outdir', help="outdir", required = True)    

    args = parsers.parse_args()
    runner = ScFEA_wilcoxon_test(args)
    runner.run()

if __name__ == '__main__':
    main()
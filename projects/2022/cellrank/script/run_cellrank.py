import sys
import argparse
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import warnings
import matplotlib.pyplot as plt

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
scv.settings.plot_prefix = ""
cr.settings.verbosity = 2

warnings.simplefilter("ignore", category=UserWarning)
warnings.simplefilter("ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=DeprecationWarning)


class Cellrank():
    def __init__(self,args):   
        self.sc_file = args.sc_file
        self.loom_file = args.loom_file
        self.cluster_key = args.cluster_key
        self.outdir = args.outdir
        self.reduction_key = args.reduction_key if args.reduction_key != None else "umap"
        self.driver_genes = args.driver_genes       
        self.n_initial_states = int(args.n_initial_states) if args.n_initial_states != None else args.n_initial_states
        self.n_terminal_states = int(args.n_terminal_states) if args.n_terminal_states != None else args.n_terminal_states
        self.save_h5ad = args.save_h5ad

        sc.settings.figdir = self.outdir
        scv.settings.figdir = self.outdir
        cr.settings.figdir = self.outdir

    def run(self):
        adata = sc.read_h5ad(self.sc_file)
        velocity_loom = scv.read(self.loom_file)
        velocity_loom.var_names_make_unique()
        adata = scv.utils.merge(adata, velocity_loom)

        # run scvelo
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
        scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
        scv.tl.recover_dynamics(adata, n_jobs=8)
        scv.tl.velocity(adata, mode="dynamical")
        scv.tl.velocity_graph(adata)
        scv.pl.velocity_embedding_stream(adata, basis=self.reduction_key, color=self.cluster_key, legend_fontsize=12, title="", 
            smooth=0.8, legend_loc="right",min_mass=4, show=False, save="scvelo.png")

        # identify terminal and initial states
        cr.tl.terminal_states(adata, cluster_key=self.cluster_key, n_states=self.n_terminal_states, weight_connectivities=0.2)
        cr.pl.terminal_states(adata, basis=self.reduction_key, show=False, save="terminal_states")
        cr.tl.initial_states(adata, cluster_key=self.cluster_key, n_states=self.n_initial_states)
        cr.pl.initial_states(adata, basis=self.reduction_key, discrete=True, show=False, save="initial_states")

        # Compute fate maps
        cr.tl.lineages(adata)
        cr.pl.lineages(adata, basis=self.reduction_key, same_plot=False, show=False, save="fate_maps")
        scv.tl.recover_latent_time(adata, root_key="initial_states_probs", end_key="terminal_states_probs")
        scv.pl.scatter(adata, color=[self.cluster_key, "latent_time"], fontsize=16, cmap="viridis", 
            perc=[2, 98], colorbar=True, rescale_color=[0, 1], title=["clusters", "scvelo latent time"],
            basis=self.reduction_key,show=False, save="scvelo_latent_time")
        scv.tl.paga(adata, groups="cell_type", root_key="initial_states_probs", end_key="terminal_states_probs",
            use_time_prior="velocity_pseudotime",)
        cr.pl.cluster_fates(adata, mode="paga_pie", cluster_key= self.cluster_key, legend_kwargs={"loc": "top right out"}, 
            legend_loc="top left out", node_size_scale=1, edge_width_scale=1, max_edge_width=4, title="directed PAGA",
            basis=self.reduction_key, show=False, save='fate_maps_in_paga')

        # Compute lineage drivers
        cr.tl.lineage_drivers(adata)
        terminal_lineage_drivers_file = f'{self.outdir}/terminal_lineage_drivers.csv'
        adata.varm['terminal_lineage_drivers'].to_csv(terminal_lineage_drivers_file)
        
        # Gene expression trends
        model = cr.ul.models.GAM(adata)
        if (self.driver_genes != None):
            cr.pl.gene_trends(adata, model=model, data_key="X", ncols=2, time_key="latent_time", same_plot=True,
                genes=tuple(pd.read_csv(self.driver_genes, header=None).iloc[0,]),  
                hide_cells=True, n_test_points=200, save="driver_genes_trends")  
        else:
            driver_genes_trends_plot = [cr.pl.gene_trends(adata, model=model, data_key="X", ncols=2, time_key="latent_time",
                genes=adata.varm['terminal_lineage_drivers'].sort_values(i+'_corr', ascending=False).index[0:4],
                same_plot=True, hide_cells=True, n_test_points=200, save=i+"_top4_driver_genes_trends") 
                for i in adata.obs['terminal_states'].cat.categories.values] 
        
        # heatmap
        driver_genes_trends_heatmap = [cr.pl.heatmap(adata,model,
            genes=adata.varm['terminal_lineage_drivers'][i+"_corr"].sort_values(ascending=False).index[:100],
            show_absorption_probabilities=True, lineages=i, n_jobs=1, backend="loky", save=i+"_driver_genes_heatmap")
            for i in adata.obs['terminal_states'].cat.categories.values]
        
        # save h5ad file
        if self.save_h5ad != None:
            del adata.uns['coarse_fwd']
            results_file = f'{self.outdir}/cellrank.h5ad'
            adata.write(results_file)

        # done
        print("##-------------------------")
        print("CellRank done.")


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--sc_file', help='h5ad from scanpy or seurat', required=True)
    parsers.add_argument('--loom_file', help='loom file from velocyto', required=True)
    parsers.add_argument('--cluster_key', help='name of variable of cell type', required=True)
    parsers.add_argument('--outdir', help='output dir', required=True)
    parsers.add_argument('--reduction_key', help='Optional, name of reduction. Default: umap')
    parsers.add_argument('--n_initial_states', help="Optional, number of initial states. Default: number from cellrank's inference")
    parsers.add_argument('--n_terminal_states', help="Optional, number of terminal states. Default: number from cellrank's inference")
    parsers.add_argument('--driver_genes', help='Optional, driver genes(.csv) for plot. Default: top4 driver genes')
    parsers.add_argument('--save_h5ad', help='Optional, save file to disk. Default: False')

    args = parsers.parse_args()
    cellrank_obj = Cellrank(args) 
    cellrank_obj.run()

if __name__ == '__main__':
    main()
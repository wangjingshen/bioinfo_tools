import sys
import argparse
import scanpy as sc
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')     # method 1
#import os
#os.environ['MPLBACKEND'] = 'Agg'   # method 2
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from cell2location.plt import plot_spatial


def add_root():
    root = Path(__file__).resolve().parent
    while not (root / "utils").exists():
        parent = root.parent
        if parent == root:
            raise FileNotFoundError("utils not found！")
        root = parent
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    return root

bioinfo_root = add_root()
pipeline_root = Path(__file__).resolve().parents[1]

from utils.utils import mkdir, logger, timer, execute_cmd

group_size = 5


@timer
def Cell2location_space_mod(space_input, sc_mod_path, spname, N_cells_per_location = 3, detection_alpha = 20, space_max_epochs = 30000):
    '''
    # the expected average cell abundance: tissue-dependent hyper-prior which can be estimated from paired histology
    # hyperparameter controlling normalisation of within-experiment variation in RNA detection
    '''    
    space_mod_path = f'{spname}/02.cell2location_space/'
    outs = f'{spname}/outs/'
    mkdir(f'{space_mod_path}')
    mkdir(f'{outs}')

    # read sc ref mod
    adata_ref = sc.read_h5ad(f"{sc_mod_path}/sc_ref.h5ad")
    mod = cell2location.models.RegressionModel.load(f"{sc_mod_path}", adata_ref)

    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        df = adata_ref.varm['means_per_cluster_mu_fg']
        df_cols = [f"means_per_cluster_mu_fg_{ct}" for ct in adata_ref.uns['mod']['factor_names']]
        df_ref = df[df_cols].copy()
        df_ref.columns = adata_ref.uns['mod']['factor_names']
    else:
        df_ref = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
            for i in adata_ref.uns['mod']['factor_names']]].copy()

    df_ref = df_ref.loc[adata_ref.var_names]

    ## read space
    adata_vis = sc.read_visium(f"{space_input}/")
    adata_vis.var_names_make_unique()
    adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]

    # remove MT genes for spatial mapping
    adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var_names]
    adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_vis.var_names, df_ref.index)
    adata_vis = adata_vis[:, intersect].copy()
    df_ref = df_ref.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata = adata_vis, batch_key = "sample")

    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df = df_ref, N_cells_per_location = N_cells_per_location, detection_alpha = detection_alpha)

    mod.train(max_epochs = space_max_epochs, batch_size = None, train_size = 1)

    # plot ELBO loss history during training, removing first 100 epochs from the plot
    mod.plot_history(1000)
    plt.legend(labels=['full data training'])

    adata_vis = mod.export_posterior(adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs})

    mod.plot_QC()
    plt.savefig(f'{space_mod_path}/QC.pdf', bbox_inches="tight", dpi=300)
    plt.close()

    # Save model and anndata object with results
    mod.save(f"{space_mod_path}", overwrite=True)
    adata_vis.write(f"{space_mod_path}/space.h5ad")

    # plot cell_abundance
    cell_abundance = adata_vis.obsm["q05_cell_abundance_w_sf"]
    cell_abundance.to_csv(f"{outs}/cell_abundance_q05.tsv", sep="\t", index=True)

    cell_abundance_columns = cell_abundance.columns.tolist()
    cluster_list = [i.split("sf_")[-1] for i in cell_abundance_columns]

    for orig_col, cluster in zip(cell_abundance_columns, cluster_list):
        adata_vis.obs[cluster] = cell_abundance[orig_col].copy()

    with mpl.rc_context({'axes.facecolor':  'black', 'figure.figsize': [4.5, 5]}):
        sc.pl.spatial(adata_vis, cmap='magma', color=cluster_list, ncols=4, size=1.3, show=False,
                      img_key='hires', vmin=0, vmax='p99.2')  # limit color scale at 99.2% quantile of cell abundance

        plt.savefig(f"{outs}/cell2location_multi_panel.pdf", bbox_inches="tight", dpi=300)

    clust_col = [str(i) for i in cluster_list]  # avoid the digital cluster, 0,1,2,3...

    for idx, start in enumerate(range(0, len(clust_col), group_size)):
        end = start + group_size
        sub_col = clust_col[start:end]
        sub_label = cluster_list[start:end]

        with mpl.rc_context({'figure.figsize': (12, 10)}):
            fig = plot_spatial(
                adata=adata_vis, color=sub_col, labels=sub_label, show_img=True, style='fast', 
                img_alpha=0.3, max_color_quantile=0.992, circle_diameter=6, colorbar_position='right'
            )
        fig.savefig(f"{outs}/cell2loc_group_{idx+1}.pdf", bbox_inches="tight", pad_inches=0.2, dpi=300)
        plt.close(fig)
        print(f"plot {idx+1} group, cluster: {sub_label}")



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--space_input', required=True, help='space input')
    parser.add_argument('--sc_mod_path', required=True, help='sc mod path')
    parser.add_argument('--spname', help='spname')
    parser.add_argument('--N_cells_per_location', default = 3, type = int, help='N_cells_per_location')
    parser.add_argument('--detection_alpha', default = 20, type = int, help='detection alpha')
    parser.add_argument('--space_max_epochs', default = 30000, type = int, help='space_max_epochs')
    args = parser.parse_args()
    
    Cell2location_space_mod(
        space_input = args.space_input, sc_mod_path = args.sc_mod_path, spname = args.spname,
        N_cells_per_location = args.N_cells_per_location, detection_alpha = args.detection_alpha, 
        space_max_epochs = args.space_max_epochs
    )


if __name__ == '__main__':
    main()
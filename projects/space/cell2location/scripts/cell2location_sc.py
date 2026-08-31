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


def add_root():
    root = Path(__file__).resolve().parent
    while not (root / "utils").exists():
        parent = root.parent
        if parent == root:
            raise FileNotFoundError("utils not found!")
        root = parent
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    return root

bioinfo_root = add_root()
pipeline_root = Path(__file__).resolve().parents[1]


from utils.utils import mkdir, logger, timer

@timer
def Cell2location_sc_mod(sc_ref, spname, labels_key = 'cluster', sc_max_epochs = 250, 
        batch_key = None, categorical_covariate_keys = None,
        cell_count_cutoff = 5, cell_percentage_cutoff2 = 0.03, nonz_mean_cutoff = 1.12):
    '''
    train sc reference signatures
    '''
    space_mod_path = f'{spname}/01.cell2location_sc/'
    mkdir(f'{space_mod_path}')
    # read data
    adata_ref = sc.read(sc_ref)

    if "raw" in adata_ref.layers:
        print(">>> Use layer 'raw' as count matrix")
        adata_ref.X = adata_ref.layers["raw"].copy()
    else:
        print(">>> layer 'raw' not found, keep original adata_ref.X")

    # filter gene
    selected = filter_genes(adata_ref, 
        cell_count_cutoff = cell_count_cutoff, 
        cell_percentage_cutoff2 = cell_percentage_cutoff2, 
        nonz_mean_cutoff = nonz_mean_cutoff)
    plt.savefig(f'{space_mod_path}/gene_filter_hist.pdf', bbox_inches="tight", dpi=300)
    plt.close()

    adata_ref = adata_ref[:, selected].copy()

    # prepare anndata for the regression model
    setup_kwargs = {"adata": adata_ref, "labels_key": labels_key}
    if batch_key is not None:
        setup_kwargs["batch_key"] = batch_key
    if categorical_covariate_keys is not None:
        setup_kwargs["categorical_covariate_keys"] = categorical_covariate_keys

    cell2location.models.RegressionModel.setup_anndata(**setup_kwargs)

    #cell2location.models.RegressionModel.setup_anndata(adata = adata_ref, labels_key = labels_key,
    #    batch_key = batch_key,categorical_covariate_keys = categorical_covariate_keys)  

    # create the regression model
    mod = RegressionModel(adata_ref)
    mod.train(max_epochs = sc_max_epochs)

    # plot
    mod.plot_history(20)
    plt.legend(labels=['full data training'])
    plt.savefig(f'{space_mod_path}/elbo_loss_curve.pdf', bbox_inches="tight", dpi=300)
    plt.close()

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})
    mod.plot_QC()
    plt.savefig(f'{space_mod_path}/QC.pdf', bbox_inches="tight", dpi=300)
    plt.close()

    # Save model
    mod.save(f'{space_mod_path}', overwrite=True)

    # Save anndata object with results
    adata_file = f'{space_mod_path}/sc_ref.h5ad'
    adata_ref.write(adata_file)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sc_ref', required=True, help='sc reference h5ad')
    parser.add_argument('--spname', help='spname')
    parser.add_argument('--labels_key', default = 'cluster', help='labels_key')
    parser.add_argument('--sc_max_epochs', default = 250, type = int, help='sc_max_epochs')
    parser.add_argument('--batch_key', default = None, help='batch_key')
    parser.add_argument('--categorical_covariate_keys', default = None, help='categorical_covariate_keys')
    parser.add_argument('--cell_count_cutoff', default = 5, type = int, help='cell_count_cutoff')
    parser.add_argument('--cell_percentage_cutoff2', default = 0.03, type = float, help='cell_percentage_cutoff2')
    parser.add_argument('--nonz_mean_cutoff', default = 1.12, type = float, help='nonz_mean_cutoff')
    args = parser.parse_args()

    if args.batch_key == "None":
        args.batch_key = None
    if args.categorical_covariate_keys == "None":
        args.categorical_covariate_keys = None


    Cell2location_sc_mod(
        sc_ref = args.sc_ref, spname = args.spname, labels_key = args.labels_key, sc_max_epochs = args.sc_max_epochs,
        batch_key = args.batch_key, categorical_covariate_keys = args.categorical_covariate_keys, 
        cell_count_cutoff = args.cell_count_cutoff, cell_percentage_cutoff2 = args.cell_percentage_cutoff2, 
        nonz_mean_cutoff = args.nonz_mean_cutoff
    )


if __name__ == '__main__':
    main()
import os
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


@timer
def Cell2location_plot(spname, group_size = 5):

    space_mod_path = f'{spname}/02.cell2location_space/'
    outs = f'{spname}/outs/'
    mkdir(f'{space_mod_path}')
    mkdir(f'{outs}')
    
    adata_vis = sc.read_h5ad(f"{space_mod_path}/space.h5ad")
    mod = cell2location.models.Cell2location.load(f"{space_mod_path}", adata_vis)

    # add 5% quantile, representing confident cell abundance, 'at least this amount is present',
    # to adata.obs with nice names for plotting
    ##adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    cell_abundance = adata_vis.obsm["q05_cell_abundance_w_sf"]
    cell_abundance.to_csv(f"{outs}/cell_abundance_q05.tsv", sep="\t", index=True)



    cell_abundance_columns = cell_abundance.columns.tolist()
    cluster_list = [i.split("sf_")[-1] for i in cell_abundance_columns]

    for orig_col, cluster in zip(cell_abundance_columns, cluster_list):
        adata_vis.obs[cluster] = cell_abundance[orig_col].copy()

    # plot in spatial coordinates
    with mpl.rc_context({'axes.facecolor':  'black', 'figure.figsize': [4.5, 5]}):
        sc.pl.spatial(adata_vis, cmap='magma', color=cluster_list, ncols=4, size=1.3, show=False,
                      img_key='hires', vmin=0, vmax='p99.2')  # limit color scale at 99.2% quantile of cell abundance

        plt.savefig(f"{outs}/cell2location_multi_panel.pdf", bbox_inches="tight", dpi=300)
        plt.savefig(f"{outs}/cell2location_multi_panel.png", bbox_inches="tight", dpi=300)

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
        fig.savefig(f"{outs}/cell2loc_group_{idx+1}.png", bbox_inches="tight", pad_inches=0.2, dpi=300)
        plt.close(fig)
        print(f"plot {idx+1} group, cluster: {sub_label}")
    
    max_idx = cell_abundance.idxmax(axis=1)
    adata_vis.obs["max_cell_type"] = [x.split("sf_")[-1] for x in max_idx]

    with mpl.rc_context({'figure.figsize': (12, 10), "axes.facecolor":"black"}):
        sc.pl.spatial(
            adata_vis,
            color="max_cell_type",
            img_key="hires",
            size=1.3,
            show=False,
        )
    plt.savefig(f"{outs}/cell2location_max_cell_type.pdf", bbox_inches="tight", dpi=300)
    plt.savefig(f"{outs}/cell2location_max_cell_type.png", bbox_inches="tight", dpi=300)
    plt.close()




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--spname', required=True, help='spname')
    parser.add_argument('--group_size', default = 5, type = int, help='group_size')
    args = parser.parse_args()
    
    Cell2location_plot(
        space_mod_path = args.spname, group_size = args.group_size)


if __name__ == '__main__':
    main()
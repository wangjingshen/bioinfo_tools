# import dependencies
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from utils import add_log

class Scanpy_run:
    def __init__(self, args):
        self.mtx = args.mtx
        self.species = args.species
        self.threads = int(args.threads)
        self.n_neighbors = int(args.n_neighbors)
        self.n_pcs = int(args.n_pcs)
        self.resolution = float(args.resolution)
        self.outdir = args.outdir
        self.name = args.name

        if not os.path.exists(self.outdir):
            os.system(f'mkdir -p {self.outdir}')
        #sc.settings.figdir = self.outdir

    @add_log
    def analysis(self):
        adata = sc.read_10x_mtx(
            self.mtx, 
            var_names = 'gene_symbols', 
            cache = False
        ) # use gene symbols for the variable names (variables-axis index)
        sc.pp.filter_cells(adata, min_genes = 0)  # simply compute the number of genes per cell (computers 'n_genes' column)
        
        if self.species == 'human':
            mito_genes = adata.var_names.str.startswith('MT-')   # mito and genes/counts cuts
        if self.species == 'mouse':
            mito_genes = adata.var_names.str.startswith('mt-') 
        adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis = 1).A1 / np.sum(adata.X, axis = 1).A1   # for each cell compute fraction of counts in mito genes vs. all genes
        adata.obs['n_counts'] = adata.X.sum(axis = 1).A1   # add the total counts per cell as observations-annotation to adata

        # initial cuts
        sc.pp.filter_cells(adata, min_genes=200 )
        sc.pp.filter_genes(adata, min_cells=3 )
        #adata = adata[adata.obs['n_genes'] < self.n_genes_t, :]
        #adata = adata[adata.obs['percent_mito'] < self.percent_mt_t, :]

        # create basic row and column attributes for the loom file:
        row_attrs = { 
            "Gene": np.array(adata.var_names) 
        }
        col_attrs = {
            "CellID": np.array(adata.obs_names) ,
            #"nGene": np.array( np.sum(adata.X.transpose()>0 , axis = 0)).flatten() ,
            #"nUMI": np.array( np.sum(adata.X.transpose() , axis = 0)).flatten() ,
        }
        lp.create(f'{self.outdir}/pyscenic_input.loom', adata.X.transpose(), row_attrs, col_attrs)

        ## sc pipeline
        # Total-count normalize (library-size correct) to 10,000 reads/cell
        sc.pp.normalize_per_cell(adata, counts_per_cell_after = 1e4)
        
        # log transform the data
        sc.pp.log1p(adata)  
        
        # identify highly variable genes
        sc.pp.highly_variable_genes(
            adata,
            layer=None,
            n_top_genes=None,
            min_disp=0.5,
            max_disp=np.inf,
            min_mean=0.0125,
            max_mean=3,
            span=0.3,
            n_bins=20,
            flavor='seurat',
            subset=False,
            inplace=True,
            batch_key=None,
            check_values=True
        )
        #adata = adata[:, adata.var['highly_variable']]
        #sc.pp.regress_out(adata, ['n_counts', 'percent_mito'], n_jobs = self.threads) #, n_jobs=args.threads)  # regress out total counts per cell and the percentage of mitochondrial genes expressed
        
        # scale each gene to unit variance, clip values exceeding SD 10
        sc.pp.scale(
            adata,
            zero_center=True,
            max_value=10,
            copy=False,
            layer=None,
            obsm=None
        )

        # principal component analysis
        sc.pp.pca(
            adata,
            n_comps=50,
            zero_center=True,
            svd_solver='auto',
            random_state=0,
            return_info=False,
            use_highly_variable=True,
            dtype='float32',
            copy=False,
            chunked=False,
            chunk_size=None
        )

        # neighborhood graph of cells
        sc.pp.neighbors(
            adata,
            n_neighbors=self.n_neighbors,
            n_pcs=self.n_pcs,
            use_rep=None,
            knn=True,
            random_state=0,
            method='umap',
            metric='euclidean',
            key_added=None,
            copy=False
        )
        
        # compute t-sne
        sc.tl.tsne(
            adata,
            n_pcs=self.n_pcs,
            n_jobs=self.threads,
            copy=False,
        )
        
        # compute UMAP
        sc.tl.umap(
            adata,
            min_dist=0.5,
            spread=1.0,
            n_components=2,
            maxiter=None,
            alpha=1.0,
            gamma=1.0,
            negative_sample_rate=5,
            init_pos='spectral',
            random_state=0,
            a=None,
            b=None,
            copy=False,
            method='umap',
            neighbors_key=None
        )

        sc.tl.leiden(
            adata,
            resolution=self.resolution,
            restrict_to=None,
            random_state=0,
            key_added='cluster',
            adjacency=None,
            directed=True,
            use_weights=True,
            n_iterations=-1,
            partition_type=None,
            neighbors_key=None,
            obsp=None,
            copy=False
        ) # cluster the neighbourhood graph

        # add cluster to h5ad, for h5ad anno
        #adata.obs['cluster_id'] = adata.obs['leiden']

        # write h5ad
        adata.write(f'{self.outdir}/{self.name}.h5ad')
    
    def run(self):
        self.analysis()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--mtx', help = 'mtx', required = True)
    parsers.add_argument('--species', help = 'species, huamn or mouse, default: human')
    parsers.add_argument('--threads', help = 'threads, default: 10')
    parsers.add_argument('--n_neighbors', help = 'n_neighbors, default: 15')
    parsers.add_argument('--n_pcs', help = 'n_pcs, default: 20')
    parsers.add_argument('--resolution', help = 'resolution, default: 0.3')
    parsers.add_argument('--outdir', help = 'outdir', required = True)
    parsers.add_argument('--name', help = 'name', required = True)
    args = parsers.parse_args()

    #if(args.n_genes_t == None):
    #    args.n_genes_t = 4000
    #if(args.percent_mt_t == None):
    #    args.percent_mt_t = 0.15
    if(args.threads == None):
        args.threads = 10
    if(args.species == None):
        args.species = 'human'
    if(args.n_neighbors == None):
        args.n_neighbors = 15
    if(args.n_pcs == None):
        args.n_pcs = 20
    if(args.resolution == None):
        args.resolution = 0.3
    runner = Scanpy_run(args) 
    runner.run()

if __name__ == '__main__':
    main()



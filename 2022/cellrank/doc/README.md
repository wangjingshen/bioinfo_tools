Here, we use cellrank to uncover cellular dynamics.

## Conda env
cellrank

## Path
/SGRNJ03/randd/user/wangjingshen/project/cellrank/cellrank_pipeline/run_cellrank.py

## Test path

/SGRNJ03/randd/user/wangjingshen/project/cellrank/cellrank_pipeline/test/test_scanpy_h5ad/
    
/SGRNJ03/randd/user/wangjingshen/project/cellrank/cellrank_pipeline/test/test_seurat_h5ad/

## Input

1) h5ad file from scanpy or seurat

2) loom file from velocyto (For example, /SGRNJ03/randd/user/wangjingshen/project/cellrank/10X_1k_mouse_brain/data/run_velocyto.sh)

3) Optional, csv file of driver genes (For example, /SGRNJ03/randd/user/wangjingshen/project/cellrank/cellrank_pipeline/test/test_scanpy_h5ad/driver_genes.csv) 


## Usage
```
python /SGRNJ03/randd/user/wangjingshen/project/cellrank/cellrank_pipeline/run_cellrank.py  \
  --sc_file /SGRNJ03/randd/user/wangjingshen/project/cellrank/10X_1k_mouse_brain/data/neuron_scanpy.h5ad \
  --loom_file /SGRNJ03/randd/user/wangjingshen/project/cellrank/10X_1k_mouse_brain/data/neuron_1k_v3_possorted_genome_sorted_UY8NC.loom \
  --cluster_key cell_type \
  --outdir results/
```
```
help info:
----
  -h, --help            show this help message and exit
  --sc_file             h5ad from scanpy or seurat
  --loom_file           loom file from velocyto
  --cluster_key         name of variable of cell type
  --outdir              dir to save fig
  --n_initial_states    Optional, number of initial states. Default: number from cellrank's inference
  --n_terminal_states   Optional, number of terminal states. Default: number from cellrank's inference  
  --driver_genes        Optional, driver genes(.csv) for plot. Default: top4 driver genes
  --save_h5ad           Optional, save file to disk. Default: False
```

## Results
```
scvelo.png                        figure of scvelo
scvelo_latent_time.pdf            figure of latent time of scvelo
initial_states.pdf                initial states infered by cellrank 
terminal_states.pdf               termainal states infered by cellrank
fate_maps.pdf                     fate maps
fate_maps_in_paga.png             fate maps in PAGA
terminal_lineage_drivers.csv      table of driver genes for all lineages
xxx_driver_genes_trends.png       dynamics of driver genes along pseudotime
xxx_driver_genes_heatmap.png      heatmap of driver genes for xxx lineage
```

# env: pyscenic_env

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import diopy
import h5py

# (seurat to) h5 to scanpy -- 
def h5_to_h5ad(h5, outdir, name):
    h5 = h5py.File(h5, 'r')
    adata = diopy.input.h5_to_adata(h5 = h5, assay_name='RNA')
    h5.close()
    adata.write_h5ad(f'{outdir}/{name}.h5ad') 

# scanpy to h5 (to seurat) --
def h5ad_to_h5(h5ad, outdir, name):
    adata = sc.read_h5ad(h5ad)
    diopy.output.write_h5(adata, file = f'{outdir}/{name}.h5', save_X=False)
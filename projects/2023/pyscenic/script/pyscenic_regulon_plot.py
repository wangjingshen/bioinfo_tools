import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import umap
import json
import base64
import zlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import natsort
from pyscenic.binarization import binarize
import math
import subprocess

class Pyscenic_regulon_plot:
    def __init__(self, args):
        self.pyscenic_metadata_loom = args.pyscenic_metadata_loom
        self.h5ad = args.h5ad
        self.network = args.network
        self.node = args.node
        self.regulons = args.regulons
        self.ncol = int(args.ncol)
        self.outdir = args.outdir
        if not os.path.exists(self.outdir):
            os.system(f'mkdir -p {self.outdir}')
        sc.settings.figdir = self.outdir

    def plot(self):
        lf = lp.connect(self.pyscenic_metadata_loom, mode='r', validate=False )
        meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
        exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
        auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
        binary_mtx, auc_thresholds = binarize( auc_mtx, num_workers = 8)

        # create a dictionary of regulons:
        regulons = {}
        for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
            regulons[i] =  list(r[r==1].index.values)
        
        # show the AUC distribution
        # select regulons:
        r = [i+'_(+)' for i in self.regulons.split(',')]
        nrow = math.ceil(round(len(r)/int(self.ncol)))    
        fig, axs = plt.subplots(nrow, self.ncol, figsize=(4*self.ncol, 4*nrow), dpi=150, sharey=False)
        for i,ax in enumerate(axs):
            if(i >= len(r)):
                break
            sns.distplot(auc_mtx[ r[i] ], ax=ax, norm_hist=True, bins=100)
            ax.plot( [ auc_thresholds[ r[i] ] ]*2, ax.get_ylim(), 'r:')
            ax.title.set_text( r[i] )
            ax.set_xlabel('')
        fig.text(-0.01, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', size='large')
        fig.text(0.5, -0.01, 'AUC', ha='center', va='center', rotation='horizontal', size='large')
        #fig.tight_layout()
        fig.savefig(f'{self.outdir}/pyscenic_auc_distribution.pdf', bbox_inches='tight')


        # tf regulon activity binary --
        # Generate a binary regulon activity matrix
        adata = sc.read_h5ad(self.h5ad)
        for i in r:
            adata.obs[i] = binary_mtx[i]
            adata.obs[i][adata.obs[i]==0] = 'off'
            adata.obs[i][adata.obs[i]==1] = 'on'
            adata.obs[i] = adata.obs[i].astype('category').cat.set_categories(['on','off'], ordered = True)
        sc.pl.umap(adata, color=r, frameon=False, palette='Set1', size=20, show = False, save='_pyscenic_regulon.pdf')
    
    def cytoscape_subset(self):
        cmd1 = (f'/SGRNJ/Public/Software/conda_env/jsr4.1/bin/Rscript /SGRNJ06/randd/USER/wangjingshen/script/pyscenic/script/regulons_cytoscape_subset.R '
                f'--network {self.network} ' 
                f'--node {self.node} ' 
                f'--regulons {self.regulons} ' 
                f'--outdir {self.outdir} ')
        subprocess.check_call(cmd1, shell = True)

    def run(self):
        self.plot()
        self.cytoscape_subset()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--pyscenic_metadata_loom', help='pyscenic add metadata loom', required=True)
    parsers.add_argument('--h5ad', help='h5ad', required=True)
    parsers.add_argument('--network', help='network', required=True)
    parsers.add_argument('--node', help='node', required=True)
    parsers.add_argument('--regulons', help='regulons', required=True)
    parsers.add_argument('--ncol', help='ncol', required=True)
    parsers.add_argument('--outdir', help='outdir', required=True)

    args = parsers.parse_args()
    runner = Pyscenic_regulon_plot(args) 
    runner.run()

if __name__ == '__main__':
    main()

# env: pyscenic_env

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import diopy
import h5py
import subprocess
from h5_h5ad_trans import h5_to_h5ad
from h5_h5ad_trans import h5ad_to_h5


class Seurat_scanpy_trans:
    def __init__(self, args):
        self.h5ad = args.h5ad
        self.rds = args.rds
        self.mode = args.mode
        self.outdir = args.outdir
        self.name = args.name

        if not os.path.exists(self.outdir):
            os.system(f'mkdir -p {self.outdir}')
    
    def run(self):
        if(self.mode == 'rds'):
            # rds to h5
            cmd1 = f'/SGRNJ/Public/Software/conda_env/jsr4.1/bin/Rscript /SGRNJ06/randd/USER/wangjingshen/script/seurat_scanpy_trans/script/h5_seurat_trans.R --rds {self.rds} --input_type rds --outdir {self.outdir} --name {self.name}'
            subprocess.check_call(cmd1, shell=True)
            h5_to_h5ad(h5=f'{self.outdir}/{self.name}.h5', outdir=self.outdir, name=self.name)

        if(self.mode == 'h5ad'):
            h5ad_to_h5(h5ad=self.h5ad, outdir=self.outdir, name=self.name)
            # h5 to rds
            if(self.outdir == '.'):
                cmd1 = f'/SGRNJ/Public/Software/conda_env/jsr4.1/bin/Rscript /SGRNJ06/randd/USER/wangjingshen/script/seurat_scanpy_trans/script/h5_seurat_trans.R --h5 {self.name}.h5 --input_type h5 --outdir {self.outdir} --name {self.name}'
            else:
                cmd1 = f'/SGRNJ/Public/Software/conda_env/jsr4.1/bin/Rscript /SGRNJ06/randd/USER/wangjingshen/script/seurat_scanpy_trans/script/h5_seurat_trans.R --h5 {self.outdir}/{self.name}.h5 --input_type h5 --outdir {self.outdir} --name {self.name}'
            subprocess.check_call(cmd1, shell=True)


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--h5ad', help = 'h5ad')
    parsers.add_argument('--rds', help = 'rds')
    parsers.add_argument('--mode', help = 'mode, rds or h5ad')
    parsers.add_argument('--outdir', help = 'outdir', required = True)
    parsers.add_argument('--name', help = 'name', required = True)
    args = parsers.parse_args()

    runner = Seurat_scanpy_trans(args) 
    runner.run()

if __name__ == '__main__':
    main()




import subprocess
import pandas as pd
import argparse
import os
import sys
from pathlib import Path
import time

from get_biotype import get_biotype_from_gtf
root = Path(__file__).resolve().parents[1]
dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from utils.utils import find_file, mkdir, logger, execute_cmd, timer

class SeuratGtf():
    def __init__(self, args):
        self.matrix_10X = args.matrix_10X
        self.spname = args.spname
        self.gname = args.gname
        self.rm_batch = args.rm_batch
        self.rm_batch_var = args.rm_batch_var
        self.species = args.species
        self.resolution = args.resolution
        self.gtf = args.gtf
        self.gene_type = args.gene_type
        self.outdir = Path(args.outdir)
        mkdir(self.outdir)

    @timer
    def get_biotype(self) -> None:
        cmd = get_biotype_from_gtf(self.gtf, self.outdir)
        execute_cmd(cmd)

    @timer
    def run_seurat(self) -> None:
        cmd = (f'Rscript {root}/script/seurat.R '
              f'--matrix_10X {self.matrix_10X} '
              f'--spname {self.spname} '
              f'--gname {self.gname} '
              f'--rm_batch {self.rm_batch} '
              f'--rm_batch_var {self.rm_batch_var} '
              f'--resolution {self.resolution} '
              f'--species {self.species} '
              f'--gtf_biotype {self.outdir}/gtf_biotype.tsv '
              f'--gene_type {self.gene_type} '
              f'--outdir {self.outdir} '
              )
        execute_cmd(cmd)

    @timer
    def run(self) -> None:
        self.get_biotype()
        self.run_seurat()

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--matrix_10X', help='matrix_10X', required=True)
    parsers.add_argument('--spname', help='spname', required=True)
    parsers.add_argument('--gname', help='gname', required=True)
    parsers.add_argument('--rm_batch', default = "F" ,help='rm_batch')
    parsers.add_argument('--rm_batch_var', default = "sample", help='rm_batch_var')
    parsers.add_argument('--resolution', default = 0.8, type = float, help='resolution, default: 0.8')
    parsers.add_argument('--species', default = "human", help='species, default: human')
    parsers.add_argument('--gtf', help='gtf', required=True)
    parsers.add_argument('--gene_type', default = "totalRNA", help='gene_type, default: totalRNA')
    parsers.add_argument('--outdir', default = "outdir", help='outdir, default: outdir')

    args = parsers.parse_args()
    runner = SeuratGtf(args) 
    runner.run()

if __name__ == '__main__':
    main()

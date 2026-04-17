import subprocess
import pandas as pd
import argparse
import os
import sys
from pathlib import Path
import time

root = Path(__file__).resolve().parents[1]
print(root)
dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
#print(sys.path)
from utils.utils import find_file, mkdir, logger, execute_cmd

class SpacePathseq():
    def __init__(self, space_dir:str, pathseq_df:str, topn_genus:str, outdir:str):
        self.spatial = Path(f'{space_dir}/outs/spatial')
        self.filter_h5 = Path(f'{space_dir}/outs/filtered_feature_bc_matrix.h5')
        self.pathseq_df = Path(f'{pathseq_df}')
        self.topn_genus = Path(f'{topn_genus}')
        self.outdir = Path(outdir)
        mkdir(self.outdir)

    def mkdir_seurat_input(self) -> None:
        mkdir(f"seurat_input/")
        cmds=[f'cp -r {self.spatial} seurat_input/spatial',
              f'cp {self.filter_h5} seurat_input',
              f'mv seurat_input/spatial/positions_list.csv seurat_input/spatial/tissue_positions_list.csv']
        for cmd in cmds:
            execute_cmd(cmd)

    def run_seurat_space(self) -> None:
        cmd = (f'Rscript {root}/script/seurat_space.R '
              f'--pathseq_df {self.pathseq_df} '
              f'--topn {self.topn_genus} '
              f'--outdir {self.outdir} ')
        execute_cmd(cmd)

    def run(self) -> None:
        self.mkdir_seurat_input()
        self.run_seurat_space()

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--space_dir', help='space_dir', required=True)
    parsers.add_argument('--pathseq_df', help='pathseq_df', required=True)
    parsers.add_argument('--topn_genus', help='filter_read_count_json file', required=True)
    parsers.add_argument('--outdir', help='outdir', required=True)

    args = parsers.parse_args()
    runner = SpacePathseq(args.space_dir, args.pathseq_df, args.topn_genus, args.outdir) 
    runner.run()

if __name__ == '__main__':
    main()

import subprocess
import pandas as pd
import argparse
import os
import sys
from pathlib import Path
import time


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

from utils.utils import find_file, mkdir, logger, execute_cmd


class SpaceTag():
    def __init__(self, space_dir:str, tag_dir:str, sample:str):
        self.spatial = Path(f'{space_dir}/outs/spatial')
        self.filter_h5 = Path(f'{space_dir}/outs/filtered_feature_bc_matrix.h5')
        self.sample = sample
        self.tag_umi = Path(f'{tag_dir}/03.count_tag/{sample}_umi_tag.tsv')


    def mkdir_seurat_input(self) -> None:
        mkdir(f"space_input/")
        cmds=[f'cp -r {self.spatial} space_input/spatial',
              f'cp {self.filter_h5} space_input',
              f'mv space_input/spatial/positions_list.csv space_input/spatial/tissue_positions_list.csv']
        for cmd in cmds:
            execute_cmd(cmd)

    def run_plot(self) -> None:
        cmd1 = (f'Rscript {pipeline_root}/scripts/plot.R '
              f'--space_input space_input/  '
              f'--tag_umi {self.tag_umi} '
              f'--outdir {self.sample} ')
        cmd2 = f'rm -rf space_input'
        execute_cmd(cmd1)
        execute_cmd(cmd2)

    def run(self) -> None:
        self.mkdir_seurat_input()
        self.run_plot()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--space_dir', help='space_dir', required=True)
    parsers.add_argument('--tag_dir', help='tag_dir', required=True)
    parsers.add_argument('--sample', help='sample', required=True)

    args = parsers.parse_args()
    runner = SpaceTag(args.space_dir, args.tag_dir, args.sample) 
    runner.run()

if __name__ == '__main__':
    main()

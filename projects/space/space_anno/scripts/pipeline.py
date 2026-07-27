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


class SpaceAnno():
    def __init__(self, space_dir:str, sc:str, score_filter:float, name:str, outdir:str):
        self.spatial = Path(f'{space_dir}/outs/spatial')
        self.filter_h5 = Path(f'{space_dir}/outs/filtered_feature_bc_matrix.h5')
        self.sc = sc
        self.score_filter = score_filter
        self.name = name
        self.outdir = Path(outdir)

    def mkdir_seurat_input(self) -> None:
        mkdir(f"space_input/")
        cmds=[f'cp -r {self.spatial} space_input/spatial',
              f'cp {self.filter_h5} space_input',
              f'mv space_input/spatial/positions_list.csv space_input/spatial/tissue_positions_list.csv']
        for cmd in cmds:
            execute_cmd(cmd)

    def run_transfer_anno(self) -> None:
        cmd = (f'Rscript {pipeline_root}/scripts/space_anno.R '
              f'--space_input space_input/  '
              f'--sc {self.sc} '
              f'--score_filter {self.score_filter} '
              f'--outdir {self.outdir} ')
        execute_cmd(cmd)

    def run_rctd_anno(self) -> None:
        mkdir(f"{self.outdir}/rctd")
        cmd1 = (f'/SGRNJ/Public/Software/conda_env/Seuratv5/bin/Rscript {pipeline_root}/scripts/Seurat_cluster_standard.R '
              f'--h5 {self.filter_h5} '
              f'--spatial {self.spatial} '
              f'--name {self.name} '
              f'--group {self.name} '
              f'--outdir {self.outdir}/rctd_rds ' )

        cmd2 = (f'/SGRNJ/Public/Software/conda_env/Seuratv5/bin/Rscript {pipeline_root}/scripts/visium_RCTDscript2.r '
            f'--sc {self.sc} '
            f'--sp {self.outdir}/rctd_rds/0.Rds/{self.name}.rds '
            f'--slot Spatial '
            f'--prefix {self.name} '
            f'--outdir {self.outdir}/rctd ' )
        execute_cmd(cmd1)
        execute_cmd(cmd2)

    def run(self) -> None:
        self.mkdir_seurat_input()
        self.run_transfer_anno()
        self.run_rctd_anno()

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--space_dir', help='space_dir', required=True)
    parsers.add_argument('--sc', help='sc', required=True)
    parsers.add_argument('--score_filter', help='score_filter', default=0)
    parsers.add_argument('--name', help='name', required=True)
    parsers.add_argument('--outdir', help='outdir', required=True)

    args = parsers.parse_args()
    runner = SpaceAnno(args.space_dir, args.sc, args.score_filter, args.name, args.outdir) 
    runner.run()

if __name__ == '__main__':
    main()

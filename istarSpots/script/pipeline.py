import os
os.environ['OMP_NUM_THREADS'] = '32'
os.environ['MKL_NUM_THREADS'] = '32'
#os.environ['OPENBLAS_NUM_THREADS'] = '1'
#os.environ['NUMEXPR_NUM_THREADS'] = '1'
import sys
import pandas as pd
import argparse
import glob
import subprocess
import logging
import json
import psutil

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

ROOT = Path(__file__).resolve().parent
#ROOT = os.path.dirname(os.path.abspath(__file__))
dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from istar2spots import istar2spots
from utils.utils import mkdir, logger, execute_cmd, timer, run_with_single_thread

CONFIG = {
    "DEVICE": "cpu",  # or "cuda"
    "PixelSize": 0.5,
    "N_HVG": 1000,
    "EPOCHS": 400,
    "FilterSize": 8,
    "MinClusterSize": 20,
}


class IStarSpots:
    def __init__(self, istar_labels: Path, dir: Path, spname: str, k: int, distance_thresh: int, clip: bool):
        self.istar_labels = istar_labels
        self.dir = dir
        self.spname = spname
        self.outdir = f'{self.spname}/outs/'
        self.space_input = f'{self.spname}/space_input/'
        self.k = k
        self.distance_thresh = distance_thresh
        self.clip = clip

        mkdir(self.outdir)
        mkdir(self.space_input)


    @timer
    def make_space_input(self) -> None:
        '''
        generate_input
        '''
        execute_cmd((f'cp {self.dir}/outs/filtered_feature_bc_matrix.h5 {self.space_input}'))
        execute_cmd((f'cp -r {self.dir}/outs/spatial {self.space_input}'))
        execute_cmd((f'mv {self.space_input}/spatial/positions_list.csv {self.space_input}/spatial/tissue_positions_list.csv'))
        self.space_pos = f'{self.space_input}/spatial/tissue_positions_list.csv'

    @timer
    def istar2spots_run(self) -> None:
        '''
        {self.outdir}/istar2spots.png
        {self.outdir}/istar2spots_df.csv
        '''
        logger.info(f"[{self.spname}] Running ln_outs step.")
        istar2spots(self.istar_labels, self.space_pos, self.outdir, self.k, self.distance_thresh, self.clip)

    @timer
    def spot_plot(self) -> None:
        execute_cmd(f'/SGRNJ/Public/Software/conda_env/r4.1_env/bin/Rscript {ROOT}/seurat_plot.R '
                    f'--space_input {self.space_input} '
                    f'--istar2spots {self.outdir}/istar2spots_df.csv '
                    f'--outdir {self.outdir} ')
    
    @timer
    def run(self) -> None:
        logger.info(f'{self.spname} start...')
        self.make_space_input()
        self.istar2spots_run()
        self.spot_plot()
        logger.info(f'{self.spname} completed.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--istar_labels', help='istar labels pickle', required=True)
    parser.add_argument('--dir', help='celescope dir', required=True)
    parser.add_argument('--spname', help='spname', required=True)
    parser.add_argument('--k', default = 3, type = int, help='k istar pixels')
    parser.add_argument('--distance_thresh', default=200, type = int, help='thresh of distance between spots and istar pixels')
    parser.add_argument('--clip', action='store_true', help='clip spots from istar') 
    args = parser.parse_args()

    runner = IStarSpots(args.istar_labels, args.dir, args.spname, args.k, args.distance_thresh, args.clip)
    runner.run()


if __name__ == '__main__':
    main()
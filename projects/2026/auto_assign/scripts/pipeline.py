import os
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

def add_root(levels_up=5):
    root = Path(__file__).resolve()
    for _ in range(levels_up):
        root = root.parent
    if not (root / "utils").exists():
        raise FileNotFoundError(f"utils not found in {root}.")
    sys.path.insert(0, str(root))
    return(root)

root_path = add_root(5)  # Top 5 parent directories of current script (bioinfo_tools)
script_path = Path(__file__).resolve().parent

from utils.utils import mkdir, logger, execute_cmd, timer, run_with_single_thread



class AutoAssign:
    def __init__(self, rds, spname, anno_var, marker_ref, step):
        self.rds = rds
        self.spname = spname
        self.anno_var = anno_var
        self.marker_ref = marker_ref
        self.step = step.strip().split(',')


    @timer
    def auto_assign(self) -> None:
        '''
        auto assign
        '''
        mkdir(self.spname)
        execute_cmd((f'Rscript {script_path}/auto_assign.R '
                f'--rds {self.rds} '
                f'--spname {self.spname} '
                f'--anno_var {self.anno_var} '
                f'--marker_ref {self.marker_ref} '
                ))
        self.anno_file = f'{self.spname}/auto_assign/test_auto_cluster_type.tsv'


    @timer
    def plot(self) -> None:
        '''
        plot
        '''
        execute_cmd((f'Rscript {script_path}/auto_assign_plot.R '
                f'--rds {self.rds} '
                f'--anno_file {self.anno_file} '
                f'--outdir {self.spname}/plot '
                ))

    
    @timer
    def run(self) -> None:
        logger.info(f'sample: {self.spname} start analysis.')
        step_order = ['auto_assign', 'plot']
        for step in step_order:
            if step in self.step:
                getattr(self, step)()
        logger.info(f'sample: {self.spname} completed.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rds', help='rds', required=True)
    parser.add_argument('--spname', help='spname', required=True)
    parser.add_argument('--anno_var', help='anno_var', required=True)
    parser.add_argument('--marker_ref', help='marker_ref', required=True)
    parser.add_argument('--step', default='auto_assign,plot', help='comma-separated step')
    args = parser.parse_args()

    runner = AutoAssign(args.rds, args.spname, args.anno_var, args.marker_ref, args.step)
    runner.run()


if __name__ == '__main__':
    main()
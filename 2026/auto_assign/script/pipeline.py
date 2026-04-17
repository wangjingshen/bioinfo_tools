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

ROOT = Path(__file__).resolve().parent
#ROOT = os.path.dirname(os.path.abspath(__file__))
dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
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
        execute_cmd((f'Rscript {ROOT}/auto_assign.R '
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
        execute_cmd((f'Rscript {ROOT}/auto_assign_plot.R '
                f'--rds {self.rds} '
                f'--anno_file {self.anno_file} '
                f'--outdir {self.spname}/plot '
                ))

    
    @timer
    def run(self) -> None:
        logger.info(f'{self.spname} start...')
        step_order = ['auto_assign', 'plot']
        for step in step_order:
            if step in self.step:
                getattr(self, step)()
        logger.info(f'{self.spname} completed.')


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
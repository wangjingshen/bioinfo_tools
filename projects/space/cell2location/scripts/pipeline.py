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

def add_root():
    root = Path(__file__).resolve().parent
    while not (root / "utils").exists():
        parent = root.parent
        if parent == root:
            raise FileNotFoundError("utils not found!")
        root = parent
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    return root

bioinfo_root = add_root()
pipeline_root = Path(__file__).resolve().parents[1]

from utils.utils import mkdir, rmdir, logger, execute_cmd, timer, run_with_single_thread, make_space_input
from cell2location_sc import Cell2location_sc_mod
from cell2location_space import Cell2location_space_mod
from cell2location_plot import Cell2location_plot


class Cell2location:
    def __init__(self, sc_ref, space_dir, spname, labels_key, N_cells_per_location, detection_alpha, sc_max_epochs, space_max_epochs,
                 batch_key, categorical_covariate_keys, cell_count_cutoff, cell_percentage_cutoff2, nonz_mean_cutoff, step):
        
        self.sc_ref = sc_ref        
        self.space_dir = space_dir
        self.spname = spname
        self.labels_key = labels_key
        self.N_cells_per_location = N_cells_per_location
        self.detection_alpha = detection_alpha
        self.sc_max_epochs = sc_max_epochs
        self.space_max_epochs = space_max_epochs
        self.batch_key = batch_key
        self.categorical_covariate_keys = categorical_covariate_keys
        self.cell_count_cutoff = cell_count_cutoff
        self.cell_percentage_cutoff2 = cell_percentage_cutoff2
        self.nonz_mean_cutoff = nonz_mean_cutoff
        self.step = step.strip().split(',')


    @timer
    def sc_train(self) -> None:
        '''
        train single cell ref mod
        '''
        Cell2location_sc_mod(
            sc_ref = self.sc_ref, labels_key = self.labels_key, sc_max_epochs = self.sc_max_epochs, spname = self.spname, 
            batch_key = self.batch_key, categorical_covariate_keys = self.categorical_covariate_keys, 
            cell_count_cutoff = self.cell_count_cutoff, cell_percentage_cutoff2 = self.cell_percentage_cutoff2, 
            nonz_mean_cutoff = self.nonz_mean_cutoff
        )

    
    @timer
    def space_train(self) -> None:
        '''
        train space mod
        '''
        self.sc_mod_path = f'{self.spname}/01.cell2location_sc'
        self.space_input = make_space_input(self.space_dir, self.spname)

        Cell2location_space_mod(
            space_input = self.space_input, sc_mod_path = self.sc_mod_path, 
            N_cells_per_location = self.N_cells_per_location, detection_alpha = self.detection_alpha, 
            space_max_epochs = self.space_max_epochs, spname = self.spname
        )

        rmdir(f'{self.space_input}')


    @timer
    def plot(self) -> None:
        '''
        plot
        '''
        Cell2location_plot(
            spname = self.spname, group_size = 5
        )


    @timer
    def run(self) -> None:
        logger.info(f'{self.spname} start...')
        step_order = ['sc_train', 'space_train', 'plot']
        for step in step_order:
            if step in self.step:
                getattr(self, step)()
        logger.info(f'{self.spname} completed.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sc_ref', required=True, help='sc reference h5ad')
    parser.add_argument('--space_dir', required=True, help='space dir')
    parser.add_argument('--spname', help='spname')
    parser.add_argument('--labels_key', default = 'cluster', help='labels_key')
    parser.add_argument('--sc_max_epochs', default = 250, type = int, help='sc_max_epochs')
    parser.add_argument('--batch_key', default = None, help='batch_key')
    parser.add_argument('--categorical_covariate_keys', nargs="+", default = None, help='categorical_covariate_keys')
    parser.add_argument('--cell_count_cutoff', default = 5, type = int, help='cell_count_cutoff')
    parser.add_argument('--cell_percentage_cutoff2', default = 0.03, type = float, help='cell_percentage_cutoff2')
    parser.add_argument('--nonz_mean_cutoff', default = 1.12, type = float, help='nonz_mean_cutoff')
    parser.add_argument('--N_cells_per_location', default = 3, type = int, help='N_cells_per_location')
    parser.add_argument('--detection_alpha', default = 20, type = int, help='detection alpha')
    parser.add_argument('--space_max_epochs', default = 30000, type = int, help='space_max_epochs')
    parser.add_argument('--step', default = 'sc_train,space_train,plot', help='step')
    args = parser.parse_args()

    if args.batch_key == "None":
        args.batch_key = None
    #if args.categorical_covariate_keys == "None":
    #    args.categorical_covariate_keys = None

    runner = Cell2location(
        sc_ref = args.sc_ref, space_dir = args.space_dir, spname = args.spname, labels_key = args.labels_key, sc_max_epochs = args.sc_max_epochs, 
        batch_key = args.batch_key,categorical_covariate_keys = args.categorical_covariate_keys, cell_count_cutoff = args.cell_count_cutoff, 
        cell_percentage_cutoff2 = args.cell_percentage_cutoff2, nonz_mean_cutoff = args.nonz_mean_cutoff, 
        N_cells_per_location = args.N_cells_per_location, detection_alpha = args.detection_alpha, space_max_epochs = args.space_max_epochs,
        step = args.step
    )
    runner.run()


if __name__ == '__main__':
    main()
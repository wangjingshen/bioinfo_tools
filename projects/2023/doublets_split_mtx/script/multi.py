import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Doublets_split_mtx_run():
    def __init__(self, mtx, doublets_species, outdir, prefix):
        self.mtx = mtx
        self.doublets_species = doublets_species
        self.outdir = outdir
        self.prefix = prefix
    
    def run(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/doublets_split_mtx/script/doublets_split_mtx.py '
            f'--mtx {self.mtx} '
            f'--doublets_species {self.doublets_species} '
            f'--outdir {self.outdir} '
            f'--prefix {self.prefix} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    mtx_list = df_mapfile['mtx']
    doublets_species_list = df_mapfile['doublets_species']
    outdir_list = df_mapfile['outdir']
    prefix_list = df_mapfile['prefix']
    return mtx_list, doublets_species_list, outdir_list, prefix_list

def run_single(mtx, doublets_species, outdir, prefix):
    runner = Doublets_split_mtx_run(mtx, doublets_species, outdir, prefix)
    runner.run()

def main():
    mapfile = sys.argv[1]
    mtx_list, doublets_species_list, outdir_list, prefix_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, mtx_list, doublets_species_list, outdir_list, prefix_list):
            print('done')


if __name__ == '__main__':
    main()
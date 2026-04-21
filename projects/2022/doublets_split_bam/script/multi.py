import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Doublets_split_bam_run():
    def __init__(self, doublets_species, fj_path, otsu_min_support_read, name, outdir):
        self.doublets_species = doublets_species
        self.fj_path = fj_path
        self.otsu_min_support_read = otsu_min_support_read
        self.name = name
        self.outdir = outdir
    
    def run(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/doublets_split_bam/doublets_split_bam.py '
            f'--doublets_species {self.doublets_species} '
            f'--fj_path {self.fj_path} '
            f'--otsu_min_support_read {self.otsu_min_support_read} '
            f'--name {self.name} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    doublets_species_list = df_mapfile['doublets_species']
    fj_path_list = df_mapfile['fj_path']
    otsu_min_support_read_list = df_mapfile['otsu_min_support_read']
    name_list = df_mapfile['name']
    outdir_list = df_mapfile['outdir']
    return doublets_species_list, fj_path_list, otsu_min_support_read_list, name_list, outdir_list

def run_single(doublets_species, fj_path, otsu_min_support_read, name, outdir):
    runner = Doublets_split_bam_run(doublets_species, fj_path, otsu_min_support_read, name, outdir)
    runner.run()

def main():
    mapfile = sys.argv[1]
    doublets_species_list, fj_path_list, otsu_min_support_read_list, name_list, outdir_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, doublets_species_list, fj_path_list, otsu_min_support_read_list, name_list, outdir_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
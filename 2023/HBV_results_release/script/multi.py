import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Hbv_results_release_run():
    def __init__(self, rds, prefix, fj_path, consensus_path, outdir, name):
        self.rds = rds
        self.prefix = prefix
        self.fj_path = fj_path
        self.consensus_path = consensus_path
        self.outdir = outdir
        self.name = name
    
    def run(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/HBV_results_release/script/analysis.py '
            f'--rds {self.rds} '
            f'--prefix {self.prefix} '
            f'--fj_path {self.fj_path} '
            f'--consensus_path {self.consensus_path} '
            f'--outdir {self.outdir} '
            f'--name {self.name} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    rds_list = df_mapfile['rds']
    prefix_list = df_mapfile['prefix']
    fj_path_list = df_mapfile['fj_path']
    consensus_path_list = df_mapfile['consensus_path']
    outdir_list = df_mapfile['outdir']
    name_list = df_mapfile['name']
    return rds_list, prefix_list, fj_path_list, consensus_path_list, outdir_list, name_list

def run_single(rds, prefix, fj_path, consensus_path, outdir, name):
    runner = Hbv_results_release_run(rds, prefix, fj_path, consensus_path, outdir, name)
    runner.run()

def main():
    mapfile = sys.argv[1]
    rds_list, prefix_list, fj_path_list, consensus_path_list, outdir_list, name_list = parse_mapfile(mapfile)
    with ProcessPoolExecutor(max_workers=1) as executor:
        executor.map(run_single, rds_list, prefix_list, fj_path_list, consensus_path_list, outdir_list, name_list)

if __name__ == '__main__':
    main()
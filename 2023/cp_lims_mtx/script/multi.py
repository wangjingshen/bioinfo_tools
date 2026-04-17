import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Cp_lims_mtx_run():
    def __init__(self, expM, name, mode):
        self.expM = expM
        self.name = name
        self.mode = mode
    
    def run(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/cp_lims_mtx/script/run.py '
            f'--expM {self.expM} '
            f'--name {self.name} '
            f'--mode {self.mode} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    expM_list = df_mapfile['expM']
    name_list = df_mapfile['name']
    mode_list = df_mapfile['mode']
    return expM_list, name_list, mode_list

def run_single(expM, name, mode):
    runner = Cp_lims_mtx_run(expM, name, mode)
    runner.run()

def main():
    mapfile = sys.argv[1]
    expM_list, name_list, mode_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        executor.map(run_single, expM_list, name_list, mode_list)

if __name__ == '__main__':
    main()
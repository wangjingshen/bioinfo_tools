import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Linker_stat_run():
    def __init__(self, R1, mismatch, outdir, name):
        self.R1 = R1
        self.mismatch = mismatch
        self.outdir = outdir
        self.name = name
    
    def run(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/linker_stat/script/analysis.py '
            f'--R1 {self.R1} '
            f'--mismatch {self.mismatch} '
            f'--outdir {self.outdir} '
            f'--name {self.name} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    R1_list = df_mapfile['R1']
    mismatch_list = df_mapfile['mismatch']
    outdir_list = df_mapfile['outdir']
    name_list = df_mapfile['name']
    return R1_list, mismatch_list, outdir_list, name_list

def run_single(R1, mismatch, outdir, name):
    runner = Linker_stat_run(R1, mismatch, outdir, name)
    runner.run()

def main():
    mapfile = sys.argv[1]
    R1_list, mismatch_list, outdir_list, name_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=3) as executor:
        executor.map(run_single, R1_list, mismatch_list, outdir_list, name_list)


if __name__ == '__main__':
    main()
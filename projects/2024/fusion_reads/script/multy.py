import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Fusion_reads_run():
    def __init__(self, path, mod, outdir):
        self.path = path
        self.mod = mod
        self.outdir = outdir
    
    def run(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/fusion_reads/script/fusion_reads.py '
            f'--path {self.path} '
            f'--mod {self.mod} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    path_list = df_mapfile['path']
    mod_list = df_mapfile['mod']
    outdir_list = df_mapfile['outdir']
    return path_list, mod_list, outdir_list

def run_single(path, mod, outdir):
    runner = Fusion_reads_run(path, mod, outdir)
    runner.run()

def main():
    mapfile = sys.argv[1]
    path_list, mod_list, outdir_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, path_list, mod_list, outdir_list):
            print('Done.')


if __name__ == '__main__':
    main()
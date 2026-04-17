import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Fusion_bam_run():
    def __init__(self, path, mod, name):
        self.path = path
        self.mod = mod
        self.name = name
    
    def run(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/fusion_bam/script/fusion_bam.py '
            f'--path {self.path} '
            f'--mod {self.mod} '
            f'--name {self.name}'
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    path_list = df_mapfile['path']
    mod_list = df_mapfile['mod']
    name_list = df_mapfile['name']
    return path_list, mod_list, name_list

def run_single(path, mod, name):
    runner = Fusion_bam_run(path, mod, name)
    runner.run()

def main():
    mapfile = sys.argv[1]
    path_list, mod_list, name_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, path_list, mod_list, name_list):
            print('Done.')


if __name__ == '__main__':
    main()
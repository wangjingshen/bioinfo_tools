import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class R1_polyT_backreads_zl_run():
    def __init__(self, r1_cele, r2_cele, name):
        self.r1_cele = r1_cele
        self.r2_cele = r2_cele
        self.name = name
    
    def run(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/R1_polyT_backreads_zl/script/R1_polyT_backreads_zl.py '
            f'--r1_cele {self.r1_cele} '
            f'--r2_cele {self.r2_cele} '
            f'--name {self.name} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    r1_cele_list = df_mapfile['r1_cele']
    r2_cele_list = df_mapfile['r2_cele']
    name_list = df_mapfile['name']
    return r1_cele_list, r2_cele_list, name_list

def run_single(r1_cele, r2_cele, name):
    runner = R1_polyT_backreads_zl_run(r1_cele, r2_cele, name)
    runner.run()

def main():
    mapfile = sys.argv[1]
    r1_cele_list, r2_cele_list, name_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, r1_cele_list, r2_cele_list, name_list):
            print('Done.')


if __name__ == '__main__':
    main()
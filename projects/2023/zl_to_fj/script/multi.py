# not use


import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    zl_mapfile_list = df_mapfile['zl_mapfile']
    fj_mapfile_list = df_mapfile['fj_mapfile']
    return zl_mapfile_list, fj_mapfile_list

def run_single(zl_mapfile, fj_mapfile):
    cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/zl_to_fj/script/analysis.py '
            f'--zl_mapfile {zl_mapfile} '
            f'--fj_mapfile {fj_mapfile} '
        )
    subprocess.check_call(cmd, shell=True)

def main():
    mapfile = sys.argv[1]
    zl_mapfile_list, fj_mapfile_list = parse_mapfile(mapfile)
    with ProcessPoolExecutor(max_workers=4) as executor:
        executor.map(run_single, zl_mapfile_list, fj_mapfile_list)

if __name__ == '__main__':
    main()
# not use


import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\s+', header=None)
    matrix_10X_list = df_mapfile.iloc[:,0].tolist()
    name_list = df_mapfile.iloc[:,1].tolist()
    return matrix_10X_list, name_list

def run_single(matrix_10X, name):
    cmd = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script/mix_determination/script/mix_determination.R '
            f'--matrix_10X {matrix_10X} '
            f'--name {name}'
        )
    subprocess.check_call(cmd, shell=True)

def main():
    mapfile = sys.argv[1]
    matrix_10X_list, name_list = parse_mapfile(mapfile)
    with ProcessPoolExecutor(max_workers=1) as executor:
        executor.map(run_single, matrix_10X_list, name_list)

if __name__ == '__main__':
    main()
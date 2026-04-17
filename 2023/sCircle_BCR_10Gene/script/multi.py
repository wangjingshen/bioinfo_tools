import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Analysis_run():
    def __init__(self, matrix_10X, fj_UMI, outdir, name):
        self.matrix_10X = matrix_10X
        self.fj_UMI = fj_UMI
        self.outdir = outdir
        self.name = name
    
    def run(self):
        cmd = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script/sCircle_BCR_10Gene/script/analyis.R '
            f'--matrix_10X {self.matrix_10X} '
            f'--fj_UMI {self.fj_UMI} '
            f'--outdir {self.outdir} '
            f'--name {self.name} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    matrix_10X_list = df_mapfile['matrix_10X']
    fj_UMI_list = df_mapfile['fj_UMI']
    outdir_list = df_mapfile['outdir']
    name_list = df_mapfile['name']
    return matrix_10X_list, fj_UMI_list, outdir_list, name_list

def run_single(matrix_10X, fj_UMI, outdir, name):
    runner = Analysis_run(matrix_10X, fj_UMI, outdir, name)
    runner.run()

def main():
    mapfile = sys.argv[1]
    matrix_10X_list, fj_UMI_list, outdir_list, name_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, matrix_10X_list, fj_UMI_list, outdir_list, name_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Fastqc_run():
    def __init__(self, fq_path, outdir):
        self.fq_path = fq_path
        self.outdir = outdir
    
    def run_fastqc(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/fastqc/fastqc.py '
            f'--fq_path {self.fq_path} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_fastqc()

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    fq_path_list = df_mapfile['fq_path']
    outdir_list = df_mapfile['outdir']

    return fq_path_list, outdir_list

def run_single(fq_path, outdir):
    runner = Fastqc_run(fq_path, outdir)
    runner.run()

def main():
    mapfile = sys.argv[1]
    fq_path_list, outdir_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, fq_path_list, outdir_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
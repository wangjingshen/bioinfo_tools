import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Hbv_virus_gene_consensus():
    def __init__(self, fj_path, outdir):
        self.fj_path = fj_path
        self.outdir = outdir
    
    def run(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/HBV_virus_gene_consensus/HBV_virus_gene_consensus.py '
            f'--fj_path {self.fj_path} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    fj_path_list = df_mapfile['fj_path']
    outdir_list = df_mapfile['outdir']
    return fj_path_list, outdir_list

def run_single(fj_path, outdir):
    runner = Hbv_virus_gene_consensus(fj_path, outdir)
    runner.run()

def main():
    mapfile = sys.argv[1]
    fj_path_list, outdir_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, fj_path_list, outdir_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
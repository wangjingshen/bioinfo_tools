import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class SAHMI_step0():
    def __init__(self, rna_path, fj_path, rna_spname, fj_spname):
        self.rna_path = rna_path
        self.fj_path = fj_path
        self.rna_spname = rna_spname
        self.fj_spname = fj_spname

    
    def run(self):
        cmd = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script/SAHMI/script/data_process.R '
            f'--rna_path {self.rna_path} '
            f'--fj_path {self.fj_path} '
            f'--rna_spname {self.rna_spname} '
            f'--fj_spname {self.fj_spname} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    rna_path_list = df_mapfile['rna_path']
    fj_path_list = df_mapfile['fj_path']
    rna_spname_list = df_mapfile['rna_spname']
    fj_spname_list = df_mapfile['fj_spname']
    return rna_path_list, fj_path_list, rna_spname_list, fj_spname_list

def run_single(rna_path, fj_path, rna_spname, fj_spname):
    runner = SAHMI_step0(rna_path, fj_path, rna_spname, fj_spname)
    runner.run()

def main():
    mapfile = sys.argv[1]
    rna_path_list, fj_path_list, rna_spname_list, fj_spname_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, rna_path_list, fj_path_list, rna_spname_list, fj_spname_list):
            print('Done.')


if __name__ == '__main__':
    main()
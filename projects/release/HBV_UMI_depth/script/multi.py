import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd


class Get_virus_bam_multi():
    def __init__(self, code_path, fj_path, otsu_min_support_read, outdir, name):
        self.code_path = code_path
        self.fj_path = fj_path
        self.otsu_min_support_read = otsu_min_support_read
        self.outdir = outdir
        self.name = name

    def run_get_virus_bam(self):
        cmd = (
            f'python {self.code_path}/get_virus_bam.py '
            f'--fj_path {self.fj_path} '
            f'--otsu_min_support_read {self.otsu_min_support_read} '
            f'--outdir {self.outdir} '
            f'--name {self.name} '
        )
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_get_virus_bam()


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\t')
    code_path_list = df_mapfile['code_path']
    fj_path_list = df_mapfile['fj_path']
    otsu_min_support_read_list = df_mapfile['otsu_min_support_read']
    outdir_list = df_mapfile['outdir']
    name_list = df_mapfile['name']
    
    return code_path_list, fj_path_list, otsu_min_support_read_list, outdir_list, name_list

def run_single(code_path, fj_path, otsu_min_support_read, outdir, name):
    runner = Get_virus_bam_multi(code_path, fj_path, otsu_min_support_read, outdir, name)
    runner.run()

def main():
    mapfile = sys.argv[1]
    code_path_list, fj_path_list, otsu_min_support_read_list, outdir_list, name_list = parse_mapfile(mapfile)
    with ProcessPoolExecutor(max_workers = 1) as executor:
        for result in executor.map(run_single, code_path_list, fj_path_list, otsu_min_support_read_list, outdir_list, name_list):
            print(result, 'done')

if __name__ == '__main__':
    main()

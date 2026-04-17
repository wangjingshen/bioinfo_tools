import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Raw_data_check_run():
    def __init__(self, download_path, raw_path, sample_name, outdir):
        self.download_path = download_path
        self.raw_path = raw_path
        self.sample_name = sample_name
        self.outdir = outdir
    
    def run(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/raw_data_check/script/check.py '
            f'--download_path {self.download_path} '
            f'--raw_path {self.raw_path} '
            f'--sample_name {self.sample_name} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    download_path_list = df_mapfile['download_path']
    raw_path_list = df_mapfile['raw_path']
    sample_name_list = df_mapfile['sample_name']
    outdir_list = df_mapfile['outdir']
    return download_path_list, raw_path_list, sample_name_list, outdir_list

def run_single(download_path, raw_path, sample_name, outdir):
    runner = Raw_data_check_run(download_path, raw_path, sample_name, outdir)
    runner.run()

def main():
    mapfile = sys.argv[1]
    download_path_list, raw_path_list, sample_name_list, outdir_list = parse_mapfile(mapfile)
    with ProcessPoolExecutor(max_workers=10) as executor:
        executor.map(run_single, download_path_list, raw_path_list, sample_name_list, outdir_list)

if __name__ == '__main__':
    main()
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd


class Subreads():
    def __init__(self, fastq_path, seed, data_size, outdir):
        self.fastq_path = fastq_path
        self.seed = seed
        self.data_size = data_size
        self.outdir = outdir

    def run_subreads(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/subreads/subreads.py '
            f'--fastq_path {self.fastq_path} '
            f'--seed {self.seed} '
            f'--data_size {self.data_size} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_subreads()


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\t')
    fastq_path_list = df_mapfile['fastq_path']
    seed_list = df_mapfile['seed']
    data_size_list = df_mapfile['data_size']
    outdir_list = df_mapfile['outdir']
    
    return fastq_path_list, seed_list, data_size_list, outdir_list


def run_single(fastq_path, seed, data_size, outdir):
    runner = Subreads(fastq_path, seed, data_size, outdir)
    runner.run()


def main():
    mapfile = sys.argv[1]
    fastq_path_list, seed_list, data_size_list, outdir_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, fastq_path_list, seed_list, data_size_list, outdir_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
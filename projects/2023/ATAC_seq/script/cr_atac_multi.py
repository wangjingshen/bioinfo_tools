import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd


class ATAC():
    def __init__(self, fastq_path, run_path, prefix, species, R2_10X, outdir):
        self.fastq_path = fastq_path
        self.run_path = run_path
        self.prefix = prefix
        self.species = species
        self.R2_10X = R2_10X
        self.outdir = outdir

    def run_atac(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/ATAC_seq/cr_atac.py '
            f'--fastq_path {self.fastq_path} '
            f'--run_path {self.run_path} '
            f'--prefix {self.prefix} '
            f'--species {self.species} '
            f'--R2_10X {self.R2_10X} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_atac()


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\t')
    fastq_path_list = df_mapfile['fastq_path']
    run_path_list = df_mapfile['run_path']
    prefix_list = df_mapfile['prefix']
    species_list = df_mapfile['species']
    R2_10X_list = df_mapfile['R2_10X']
    outdir_list = df_mapfile['outdir']
    
    return fastq_path_list, run_path_list, prefix_list, species_list, R2_10X_list, outdir_list


def run_single(fastq_path, run_path, prefix, species, R2_10X, outdir):
    runner = ATAC(fastq_path, run_path, prefix, species, R2_10X, outdir)
    runner.run()


def main():
    mapfile = sys.argv[1]
    fastq_path_list, run_path_list, prefix_list, species_list, R2_10X_list, outdir_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=3) as executor:
        for result in executor.map(run_single, fastq_path_list, run_path_list, prefix_list, species_list, R2_10X_list, outdir_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
import pysam
import pandas as pd
from collections import defaultdict
import pyranges as pr
from functools import reduce
import glob
import os
import argparse
import subprocess

class Hbv_results_release:
    def __init__(self, args):
        self.rds = args.rds
        self.prefix = args.prefix
        self.virus_tsne = glob.glob(f"{args.fj_path}/06.analysis_capture_virus/*_virus_tsne.tsv")[0]
        self.virus_matrix = glob.glob(f"{args.fj_path}/08.count_capture_virus_mtx/*_virus_matrix/")[0]
        self.HBV_positive_stat = glob.glob(f"{args.fj_path}/06.analysis_capture_virus/stat.txt")[0]
        self.consensus_path = glob.glob(f'{args.consensus_path}/3.consensus')[0]
        self.outdir = args.outdir
        self.name = args.name
    
    def mkdir(self):
        if not os.path.exists(f"{self.outdir}/01/{self.name}"):
            os.system(f"mkdir -p {self.outdir}/01/{self.name}")
        if not os.path.exists(f"{self.outdir}/02/{self.name}"):
            os.system(f"mkdir -p {self.outdir}/02/{self.name}")
        if not os.path.exists(f"{self.outdir}/03/{self.name}"):
            os.system(f"mkdir -p {self.outdir}/03/{self.name}")
        if not os.path.exists(f"{self.outdir}/04/{self.name}"):
            os.system(f"mkdir -p {self.outdir}/04/{self.name}")
        if not os.path.exists(f"{self.outdir}/05/{self.name}"):
            os.system(f"mkdir -p {self.outdir}/05/{self.name}")

    def plot_and_stat(self):
        '''
        02 HBV status plot
        02 HBV genes sum stat
        03 HBV virus genes plot
        04 cell type HBV percent stat 
        '''
        cmd = (f'/SGRNJ/Public/Software/conda_env/r4.1_env/bin/Rscript '
               f'/SGRNJ06/randd/USER/wangjingshen/script/HBV_results_release/script/HBV_gene_plot_stat.R '
               f'--rds {self.rds} '
               f'--prefix {self.prefix} '
               f'--virus_tsne {self.virus_tsne} '
               f'--virus_matrix {self.virus_matrix} '
               f'--HBV_positive_stat {self.HBV_positive_stat} '
               f'--outdir_02 {self.outdir}/02/{self.name} '
               f'--outdir_03 {self.outdir}/03/{self.name} '
               f'--outdir_04 {self.outdir}/04/{self.name} '
               f'--name {self.name}'
        )
        subprocess.check_call(cmd, shell = True)
    
    def cp_virus_mtx_and_consensus(self):
        '''
        virus mtx to 01
        consensus to 05
        '''
        cmd1 = f'cp -r {self.virus_matrix} {self.outdir}/01/{self.name}/'
        subprocess.check_call(cmd1, shell = True)
        cmd2 = f'cp -r {self.consensus_path} {self.outdir}/05/{self.name}/'
        subprocess.check_call(cmd2, shell = True)

    
    def run(self):
        self.mkdir()
        self.plot_and_stat()
        self.cp_virus_mtx_and_consensus()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--rds', help = 'rds', required = True)
    parsers.add_argument('--prefix', help = 'prefix', required = True)
    parsers.add_argument('--fj_path', help = 'fj path', required = True)
    parsers.add_argument('--consensus_path', help = 'consensus path', required = True)
    parsers.add_argument('--outdir', help = 'outdir', required = True)
    parsers.add_argument('--name', help = 'name', required=True)

    args = parsers.parse_args()
    runner = Hbv_results_release(args) 
    runner.run()

if __name__ == '__main__':
    main()
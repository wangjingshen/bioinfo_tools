#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import glob


class Bam_primer():
    def __init__(self, bam, mode, name, top_n):
        '''
        '''
        self.bam = bam
        self.mode = mode
        self.name = name
        self.top_n = int(top_n)


    def mkdir(self):
        cmd1 = f"mkdir -p {self.name}"
        subprocess.check_call(cmd1, shell=True)

    def run_get_top_bp(self):
        cmd2 = (f'samtools view {self.bam} | '
                r'''awk -F "\t" '{print $10}' | '''
                f'cut -c -{self.top_n} > {self.name}/{self.name}_top{self.top_n}bp.tsv')
        subprocess.check_call(cmd2, shell=True)

    def run_get_top_bp_host(self):
        '''
        # https://www.biostars.org/p/56246/
        '''
        cmd3 = f"samtools view -b -F 4 {self.bam} > {self.name}/{self.name}_mapped.bam"
        cmd4 = (f'samtools view {self.name}/{self.name}_mapped.bam | '
                r'''awk -F "\t" '{print $10}' | '''
                f'cut -c -{self.top_n} > {self.name}/{self.name}_top{self.top_n}bp.tsv')
        subprocess.check_call(cmd3, shell=True)
        subprocess.check_call(cmd4, shell=True)


    def run_stat(self):
        cmd1 = (f'Rscript /SGRNJ06/randd/USER/wangjingshen/script/sc16S_bam_primer/script/stat.R '
               f'--df {self.name}/{self.name}_top{self.top_n}bp.tsv '
               f'--name {self.name} '
               f'--top_n {self.top_n} '
               f'--outdir {self.name} ' )
        subprocess.check_call(cmd1, shell=True)

    def run(self):
        self.mkdir()
        if(self.mode == "host"):
            self.run_get_top_bp_host()
        else:
            self.run_get_top_bp()
        self.run_stat()

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\s+', header = 0)
    bam_list = df_mapfile['bam']
    mode_list = df_mapfile['mode']
    name_list = df_mapfile['name']
    top_n_list = df_mapfile['top_n']
    return bam_list, mode_list, name_list, top_n_list


def run_single(bam, mode, name, top_n):
    runner = Bam_primer(bam, mode, name, top_n)
    runner.run()


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='mapfile', required=True)
    # add version
    parser.add_argument('--version', action='version', version='1.0')
    
    args = parser.parse_args()
    bam_list, mode_list, name_list, top_n_list = parse_mapfile(args.mapfile)

    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, bam_list, mode_list, name_list, top_n_list)

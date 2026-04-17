#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
import pandas as pd


class SC16S_analysis_pipeline():
    def __init__(self, args):
        '''
        '''
        self.matrix_10X = args.matrix_10X
        self.species = args.species
        self.rna_spname = args.rna_spname
        self.gname = args.gname
        self.fj_path = args.fj_path
        self.fj_spname = args.fj_spname
        self.step = args.step.strip().split(',')


    def anno(self):
        '''
        '''
        cmd1 = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script/seurat/script/seurat.R '
            f'--matrix_10X {self.matrix_10X} '
            f'--species {self.species} '
            f'--spname {self.rna_spname} '
            f'--gname {self.gname} '
            f'--outdir 00.rds '
        )

        subprocess.check_call(cmd1, shell=True)


    def filter(self):
        '''
        '''
        cmd1 = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/2024/script/pathseq_environment_filter/script/analysis_v2.R '
            f'--fj_path {self.fj_path} '
            f'--sample {self.fj_spname} '
        )
        subprocess.check_call(cmd1, shell=True)


    def analysis_raw(self):
        '''
        raw sc16S UMI analysis
        '''
        cmd1 = (
        # raw
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/2024/script/pathseq_down/script/analysis_v2.R '
            f'--rds 00.rds/seurat.rds '
            f'--df_genus_path 01.raw_mtx/ '
            f'--rna_spname {self.rna_spname} '
            f'--spname {self.fj_spname} '
            f'--barplot_topn 10 '
            f'--barplot_width 10 '
            f'--outdir 03.analysis_raw '
            f'--saveRDS F '     
        )
        subprocess.check_call(cmd1, shell=True)


    def analysis_filter(self):
        '''
        filter sc16S UMI analysis
        '''
        cmd1 = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/2024/script/pathseq_down/script/analysis_v2.R '
            f'--rds 00.rds/seurat.rds '
            f'--df_genus_path 02.filter_mtx/ '
            f'--rna_spname {self.rna_spname} '
            f'--spname {self.fj_spname} '
            f'--barplot_topn 10 '
            f'--barplot_width 10 '
            f'--outdir 04.analysis_filter '
            f'--saveRDS F '     
        )
        subprocess.check_call(cmd1, shell=True)


    def run(self):
        if 'anno' in self.step:
            self.anno()
        if 'filter' in self.step:
            self.filter()
        if 'analysis_raw' in self.step:
            self.analysis_raw()
        if 'analysis_filter' in self.step:
            self.analysis_filter()


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix_10X', help='matrix_10X', required=True)
    parser.add_argument('--species', help='species', required=True)
    parser.add_argument('--rna_spname', help='rna spname', required=True)
    parser.add_argument('--gname', help='gname', required=True)
    parser.add_argument('--fj_path', help='fj path', required=True)
    parser.add_argument('--fj_spname', help='fj spname', required=True)
    parser.add_argument('--step')

    # add version
    parser.add_argument('--version', action='version', version='1.0')
    args = parser.parse_args()

    if(args.step is None): # None run all steps
        args.step = 'anno,filter,analysis_raw,analysis_filter'

    runner = SC16S_analysis_pipeline(args)
    runner.run()
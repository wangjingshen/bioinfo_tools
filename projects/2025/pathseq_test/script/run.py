#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import glob

JAVA_options = '-Xmx60G'
MICROBE_bwa_image = '/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16s_sgr/reference/pathseq_microbe/pathseq_microbe.fa.img'
MICROBE_dict = '/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16s_sgr/reference/pathseq_microbe/pathseq_microbe.dict'
MICROBE_taxonomy_file = '/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16s_sgr/reference/pathseq_microbe/pathseq_taxonomy.db'
MIN_clipped_read_length = 70
TMP_dir = '/SGRNJ06/randd/USER/wangjingshen/tmp/'


class Pathseq_test():
    def __init__(self, fq, sample, host_species):
        self.fq = fq
        self.sample = sample
        self.host_species = host_species
        if self.host_species == "human":
            self.filter_bwa_image = '/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/reference/pathseq_host/pathseq_host.fa.img'
            self.kmer_file = '/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/reference/pathseq_host/pathseq_host.bfi'
        if self.host_species == "mouse":
            self.filter_bwa_image = '/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/reference/pathseq_mm10/mm10.fasta.img'
            self.kmer_file = '/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/reference/pathseq_mm10/mm10.bfi'

    def run_fastqtosam(self):
        cmd1 = f"mkdir -p {self.sample}/"
        cmd2 = (
            f'picard FastqToSam '
            f'F1={self.fq} '
            f'O={self.sample}/{self.sample}.bam '
            f'SM={self.sample} '
            f'RG={self.sample} '
            f'TMP_DIR={TMP_dir} '
        )
        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)

    def run_pathseq(self):
        cmd1 = (
            f'gatk PathSeqPipelineSpark '
            f'--input {self.sample}/{self.sample}.bam ' 
            f'--filter-bwa-image {self.filter_bwa_image} '
            f'--kmer-file {self.kmer_file} '
            f'--min-clipped-read-length {MIN_clipped_read_length} '
            f'--microbe-bwa-image {MICROBE_bwa_image} '
            f'--microbe-dict {MICROBE_dict} '
            f'--taxonomy-file {MICROBE_taxonomy_file} '
            f'--divide-by-genome-length true '
            f'--java-options {JAVA_options} '
            f'--tmp-dir {TMP_dir} '
            f'--is-host-aligned false '
            f'--output {self.sample}/{self.sample}_pathseq.bam '
            f'--scores-output {self.sample}/{self.sample}_pathseq_score.txt '
        )
        subprocess.check_call(cmd1, shell=True)


    def run(self):
        self.run_fastqtosam()
        self.run_pathseq()

def parse_mapfile(mapfile, host_species):
    df_mapfile = pd.read_csv(mapfile, sep='\t', header = None)
    df_mapfile.columns = ["fq", "sample"]
    df_mapfile['host_species'] = host_species

    fq_list = df_mapfile['fq']
    sample_list = df_mapfile['sample']
    host_species_list = df_mapfile['host_species']
    
    return fq_list, sample_list, host_species_list


def run_single(fq, sample, host_species):
    runner = Pathseq_test(fq, sample, host_species)
    runner.run()


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='mapfile', required=True)
    parser.add_argument('--host_species', required=True)

    # add version
    parser.add_argument('--version', action='version', version='1.0')
    
    args = parser.parse_args()
    fq_list, sample_list, host_species_list = parse_mapfile(args.mapfile, args.host_species)    

    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, fq_list, sample_list, host_species_list)
#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import glob


class Unmap_reads():
    def __init__(self, sample, fq_dir, library_id, match_dir, host_species):
        '''
        fq_dir: fastq dir
        library_id: library id
        match_dir: match rna dir
        host_species: host_species
        '''
        self.sample = sample
        self.fq_dir = fq_dir
        self.library_id = library_id
        self.fq1 = glob.glob(f"{fq_dir}/{library_id}*R1*gz")[0]
        self.fq2 = glob.glob(f"{fq_dir}/{library_id}*R2*gz")[0]
        self.match_dir = match_dir
        self.match_bc = f"{match_dir}/outs/filtered/barcodes.tsv.gz"

        self.host_species = host_species
        if self.host_species == "human":
            self.genomeDir = '/SGRNJ06/randd/public/genome/rna/celescope_v2/hs/'

        if self.host_species == "mouse":
            self.genomeDir = '/SGRNJ06/randd/public/genome/rna/celescope_v2/mmu/'
   

    def run_starsolo(self):
        '''
        cmd1: run starsolo; Within for get map and unmap reads
        cmd2: get unmap reads
        cmd3: subset reads
        cmd4: addRG
        cmd5: del
        '''
        cmd1 = (
            f'/SGRNJ06/randd/USER/wangjingshen/script/unmap_fq/script/starsolo.py '
            f'--subparser_assay sc16s '
            f'--outdir {self.sample}/01.starsolo/ '
            f'--sample {self.sample} '
            f'--thread 8 '
            f'--chemistry auto '
            f'--adapter_3p AAAAAAAAAAAA '
            f'--genomeDir {self.genomeDir} '
            f'--outFilterMatchNmin 50 '
            f'--soloCellFilter "EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.001 10000" '
            f'--starMem 32 '
            f'--soloFeatures "Gene GeneFull_Ex50pAS" '
            f'--fq1 {self.fq1} '
            f'--fq2 {self.fq2} '
            f'--STAR_param "--outSAMunmapped Within" '
        )

        cmd2 = f'samtools view -b -h -f 4 {self.sample}/01.starsolo/outs/{self.sample}_Aligned.sortedByCoord.out.bam > {self.sample}/01.starsolo/outs/{self.sample}_unmapped.bam '  # select unmap reads

        cmd3 = (
            f'samtools view {self.sample}/01.starsolo/outs/{self.sample}_unmapped.bam | '
            '''awk -F "\t" '{print $1}' > '''
            f'{self.sample}/01.starsolo/outs/{self.sample}_unmapped_name.txt'
        )

        cmd4 = f"mkdir -p {self.sample}/02.unmapped_fq/"

        cmd5 = f'seqtk subseq {self.fq1} {self.sample}/01.starsolo/outs/{self.sample}_unmapped_name.txt | gzip > {self.sample}/02.unmapped_fq/{self.sample}_unmapped_R1.fq.gz'
        cmd6 = f'seqtk subseq {self.fq2} {self.sample}/01.starsolo/outs/{self.sample}_unmapped_name.txt | gzip > {self.sample}/02.unmapped_fq/{self.sample}_unmapped_R2.fq.gz'

        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)
        subprocess.check_call(cmd3, shell=True)    
        subprocess.check_call(cmd4, shell=True)
        subprocess.check_call(cmd5, shell=True)
        subprocess.check_call(cmd6, shell=True)

    def make_unmap_mapfile(self):
        '''
        mapfile has four columns, which are library id, fq dir, sample, and match RNA dir in order.
        '''
        df_mapfile = pd.read_csv(self.mapfile, sep='\t', header = None)
        df_mapfile.columns = ["library_id", "fq_dir", "sample", "match_dir"]
        
        current_path = os.getcwd()
        df_mapfile['fq_dir'] = f'{current_path}/{self.sample}/02.unmapped_fq/'
        df_mapfile.to_csv(f'{self.sample}/mapfile', sep="\t", index=False, header=False)


    def run(self):
        self.run_starsolo()
        self.make_unmap_mapfile()


def parse_mapfile(mapfile, host_species):
    '''
    mapfile has four columns, which are library id, fq dir, sample, and match RNA dir in order.
    '''
    df_mapfile = pd.read_csv(mapfile, sep='\t', header = None)
    df_mapfile.columns = ["library_id", "fq_dir", "sample", "match_dir"]
    df_mapfile['host_species'] = host_species

    sample_list = df_mapfile['sample']
    fq_dir_list = df_mapfile['fq_dir']
    library_id_list = df_mapfile['library_id']
    match_dir_list = df_mapfile['match_dir']
    host_species_list = df_mapfile['host_species']

    return sample_list, fq_dir_list, library_id_list, match_dir_list, host_species_list


def run_single(sample, fq_dir, library_id, match_dir, host_species):
    runner = Unmap_reads(sample, fq_dir, library_id, match_dir, host_species)
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
    sample_list, fq_dir_list, library_id_list, match_dir_list, host_species_list = parse_mapfile(args.mapfile, args.host_species)    
    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, sample_list, fq_dir_list, library_id_list, match_dir_list, host_species_list)

#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import glob

MICROBE_bwa_image = '/SGRNJ06/randd/USER/wangjingshen/rd_project/invadeseq/reference/best_practice_data/resources/pathseq_microbe/pathseq_microbe.fa.img'
MICROBE_dict = '/SGRNJ06/randd/USER/wangjingshen/rd_project/invadeseq/reference/best_practice_data/resources/pathseq_microbe/pathseq_microbe.dict'
MICROBE_taxonomy_file = '/SGRNJ06/randd/USER/wangjingshen/rd_project/invadeseq/reference/best_practice_data/resources/pathseq_microbe/pathseq_taxonomy.db'
MIN_clipped_read_length = 70
TMP_dir = '/SGRNJ06/randd/USER/wangjingshen/tmp/'
ASSETS_dir = '/SGRNJ06/randd/USER/wangjingshen/rd_project/scsweetseq/assets/'
JAVA_mem = '-Xmx48G'


class SC_16S():
    def __init__(self, sample, fq_dir, library_id, genomeDir, host_species, match_dir):
        '''
        fq_dir: fastq dir
        library_id: library id
        genomeDir: genomeDir
        host_species: host_species
        '''
        self.sample = sample
        self.fq_dir = fq_dir
        self.library_id = library_id
        self.fq1 = glob.glob(f"{fq_dir}/{library_id}*1*gz")[0]
        self.fq2 = glob.glob(f"{fq_dir}/{library_id}*2*gz")[0]
        self.genomeDir = genomeDir
        self.match_bc = f"{match_dir}/outs/filtered/barcodes.tsv.gz")

        self.host_species = host_species
        if self.host_species == "human":
            self.filter_bwa_image = '/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/reference/pathseq_host/pathseq_host.fa.img'
            self.kmer_file = '/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/reference/pathseq_host/pathseq_host.bfi'
        if self.host_species == "mouse":
            self.filter_bwa_image = '/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/reference/pathseq_mm10/mm10.fasta.img'
            self.kmer_file = '/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/reference/pathseq_mm10/mm10.bfi'
    

    def run_barcode(self):
        cmd = (
            f'/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/script/barcode.py '
            f'--sample {self.sample} '
            f'--assay sc_16S '
            f'--protocol auto '
            f'--assets_dir {ASSETS_dir} '
            f'--fq1 {self.fq1} '
            f'--fq2 {self.fq2} '
        )
        subprocess.check_call(cmd, shell=True)

    def run_cutadapt(self):
        cmd = (
            f'/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/script/cutadapt.py '
            f'--sample {self.sample} '
            f'--fq {self.sample}_2.fq '
        )
        subprocess.check_call(cmd, shell=True)

    def run_star(self):
        cmd = (
            f'/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/script/star_mixin.py '
            f'--sample {self.sample} '
            f'--genomeDir {self.genomeDir} '
            f'--fq {self.sample}_clean_2.fq '
            f'--STAR_param "--outSAMunmapped Within" '
        )
        subprocess.check_call(cmd, shell=True)

    def sort_add_group(self):
        cmd1 = f'samtools sort -n {self.sample}_Aligned.out.bam -@ 4 -o {self.sample}_Aligned.out.sort.bam '
        cmd2 = (
            f'picard AddOrReplaceReadGroups '
            f'INPUT={self.sample}_Aligned.out.sort.bam '
            f'OUTPUT={self.sample}_Aligned.out.sort.addRG.bam '
            f'RGID={self.sample}_rgid '
            f'RGLB={self.sample}_rglb '
            f'RGPL=illumina '
            f'RGPU={self.sample}_rgpu '
            f'RGSM={self.sample}_rgsm '
        )
        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)

    def run_pathseq(self):
        cmd = (
            f'gatk PathSeqPipelineSpark '
            f'--input {self.sample}_Aligned.out.sort.addRG.bam ' 
            f'--filter-bwa-image {self.filter_bwa_image} '
            f'--kmer-file {self.kmer_file} '
            f'--min-clipped-read-length {MIN_clipped_read_length} '
            f'--microbe-bwa-image {MICROBE_bwa_image} '
            f'--microbe-dict {MICROBE_dict} '
            f'--taxonomy-file {MICROBE_taxonomy_file} '
            f'--divide-by-genome-length true '
            f'--java-options {JAVA_mem} '
            f'--tmp-dir {TMP_dir} '
            f'--is-host-aligned false '
            f'--output {self.sample}_pathseq.bam '
            f'--scores-output {self.sample}_pathseq_score.txt '
        )
        subprocess.check_call(cmd, shell=True)


    def run_UMI_matrix(self):
        cmd = (
            f'/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/script/UMI_matrix.py '
            f'{self.sample}_Aligned.out.sort.bam '  # here, use sort bam
            f'{self.sample} '
            f'{self.match_bc}'
            f'{self.sample}_pathseq.bam '
            f'{self.sample}_pathseq_score.txt '
            f'{self.sample}.invadeseq.readname '
            f'{self.sample}.invadeseq.unmap_cbub.bam '
            f'{self.sample}.invadeseq.unmap_cbub.fasta '
            f'{self.sample}.invadeseq.list '
            f'{self.sample}.invadeseq.raw.readnamepath '
            f'{self.sample}.invadeseq.genus.cell '
            f'{self.sample}.invadeseq.genus.csv '
            f'{self.sample}.invadeseq.validate.csv '
        )
        subprocess.check_call(cmd, shell=True)


    def run_analysis(self):
        cmd = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script/PathSeq/script/UMI_df.R '
            f'--sample {self.sample} '
            f'--rna_dir {self.rna_dir}  '
            f'--pathseq_dir {self.pathseq_dir} '
        )
        subprocess.check_call(cmd, shell=True)


    def run(self):
        self.run_barcode()
        self.run_cutadapt()
        self.run_star()
        self.unmap_add_group()
        self.run_pathseq()
        self.run_UMI_matrix()
        self.run_analysis()

def parse_mapfile(mapfile, genomeDir, host_species):
    df_mapfile = pd.read_csv(mapfile, sep='\t')
    df_mapfile['genomeDir'] = genomeDir
    df_mapfile['host_species'] = host_species
    sample_list = df_mapfile['sample']
    fq_dir_list = df_mapfile['fq_dir']
    library_id_list = df_mapfile['library_id']
    genomeDir_list = df_mapfile['genomeDir']
    host_species_list = df_mapfile['host_species']
    
    return sample_list, fq_dir_list, library_id_list, genomeDir_list, host_species_list


def run_single(sample, fq_dir, library_id, genomeDir, host_species):
    runner = SC_16S(sample, fq_dir, library_id, genomeDir, host_species)
    runner.run()


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='mapfile', required=True)
    parser.add_argument('--genomeDir', help='genomeDir.', required=True)
    parser.add_argument('--host_species', required=True)

    # add version
    parser.add_argument('--version', action='version', version='1.0')
    
    args = parser.parse_args()
    sample_list, fq_dir_list, library_id_list, genomeDir_list, host_species_list = parse_mapfile(args.mapfile, args.genomeDir, args.host_species)
    #genomeDir_list = list(args.genomeDir * len(sample_list))
    #host_species_list = list(args.host_species * len(sample_list))
    #prin(sample_list)
    #print(genomeDir_list)
    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, sample_list, fq_dir_list, library_id_list, genomeDir_list, host_species_list)
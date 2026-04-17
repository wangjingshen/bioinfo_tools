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


class Pathseq_run():
    def __init__(self, sample, fj_dir, pathseq_dir, host_species):
        '''
        fj_dir: fj run celescope rna with --STAR_param "--outSAMunmapped Within"
        rna_dir: rna run celescope rna
        pathseq_dir: pathseq dir
        sample: sample name
        '''
        self.sample = sample
        self.fj_dir = fj_dir
        self.fj_map_bam = glob.glob(f"{fj_dir}/outs/*Aligned.sortedByCoord.out.bam")[0]
        self.pathseq_dir = pathseq_dir
        self.host_species = host_species
        if self.host_species == "human":
            self.filter_bwa_image = '/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/reference/pathseq_host/pathseq_host.fa.img'
            self.kmer_file = '/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/reference/pathseq_host/pathseq_host.bfi'
        if self.host_species == "mouse":
            self.filter_bwa_image = '/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/reference/pathseq_mm10/mm10.fasta.img'
            self.kmer_file = '/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/reference/pathseq_mm10/mm10.bfi'

    def add_group(self):
        cmd = (
            f'picard AddOrReplaceReadGroups '
            f'INPUT = {self.fj_map_bam} '
            f'OUTPUT = {self.pathseq_dir}/{self.sample}_addRG.bam '
            f'RGID={self.sample}_rgid '
            f'RGLB={self.sample}_rglb '
            f'RGPL=illumina '
            f'RGPU={self.sample}_rgpu '
            f'RGSM={self.sample}_rgsm '
        )
        subprocess.check_call(cmd, shell=True)

    
    def run_pathseq(self):
        cmd = (
            f'gatk PathSeqPipelineSpark '
            f'--input {self.pathseq_dir}/{self.sample} ' 
            f'--filter-bwa-image /SGRNJ06/randd/USER/wangjingshen/rd_project/invadeseq/reference/mm10/pathseq_mm10/mm10.fasta.img '
            f'--kmer-file /SGRNJ06/randd/USER/wangjingshen/rd_project/invadeseq/reference/mm10/pathseq_mm10/mm10.bfi '
            f'--min-clipped-read-length {MIN_clipped_read_length} '
            f'--microbe-bwa-image {MICROBE_bwa_image} '
            f'--microbe-dict {MICROBE_dict} '
            f'--taxonomy-file {MICROBE_taxonomy_file} '
            f'--divide-by-genome-length true '
            f'--java-options "-Xmx48G" '
            f'--tmp-dir {TMP_dir} '
            f'--is-host-aligned false '
            f'--output{self.pathseq_dir}/{self.sample}_pathseq.bam '
            f'--scores-output {self.pathseq_dir}/{self.sample}_pathseq_score.txt '
        )
        subprocess.check_call(cmd, shell=True)

   def run(self):
        self.add_group()
        self.run_pathseq()


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\t')
    sample_list = df_mapfile['sample']
    fj_dir_list = df_mapfile['fj_dir']
    rna_dir_list = df_mapfile['rna_dir']
    pathseq_dir_list = df_mapfile['pathseq_dir']
    
    return sample_list, fj_dir_list, rna_dir_list, pathseq_dir_list


def run_single(sample, fj_dir, rna_dir, pathseq_dir):
    runner = UMI_matrix_run(sample, fj_dir, rna_dir, pathseq_dir)
    runner.run()


def main():
    mapfile = sys.argv[1]
    sample_list, fj_dir_list, rna_dir_list, pathseq_dir_list = parse_mapfile(mapfile)
    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, sample_list, fj_dir_list, rna_dir_list, pathseq_dir_list)
        #for result in executor.map(run_single, sample_list, fj_dir_list, rna_dir_list, pathseq_dir_list):
        #    print('Done.')

if __name__ == '__main__':
    main()
from pathlib import Path
import os
import subprocess
import pysam
from collections import defaultdict
import pandas as pd
import glob
import argparse
import sys
sys.path.append('/SGRNJ06/randd/USER/wangjingshen/script/tools')
import utils as utils


def bam2fq(bam, outdir, name):
    cmd1 = f'bedtools bamtofastq -i {bam} -fq {outdir}/{name}.fq'
    subprocess.check_call(cmd1, shell=True)

def fq2fa(fq, outdir, name):
    cmd1 = f'/SGRNJ/Public/Software/conda_env/cellrank/bin/seqtk seq -A {fq} > {outdir}/{name}.fa'
    subprocess.check_call(cmd1, shell=True)

def star(fq, genome, outdir, name):
    '''
    outdir: sample_name/05.fusion/
    '''
    cmd = (
        f'/SGRNJ/Public/Software/conda_env/celescope1.5.1b0/bin/STAR --runThreadN 4 '
        f'--genomeDir {genome} '
        f'--readFilesIn {fq} '
        f'--outFilterMultimapNmax 1 '
        f'--outFileNamePrefix {outdir}/{name}_ '
        f'--outSAMtype BAM Unsorted ' 
    )
    subprocess.check_call(cmd, shell=True)
    #f'--outFilterMatchNmin 50 '  # 原版是放进去的

class Fusion_bam:
    def __init__(self, args):
        self.bam = glob.glob(f"{args.path}/04.star_virus/*_Aligned.sortedByCoord.out.bam")[0]
        self.mod = args.mod
        self.name = args.name
        self.outdir = f'{args.name}_fusion'

        if(self.mod == "human"):
            self.fusion_human_name = f'{self.name}_fusion_human_HBV'
        if(self.mod == "homo_mus"):
            self.fusion_human_name = f'{self.name}_fusion_human_HBV'
            self.fusion_mouse_name = f'{self.name}_fusion_mouse_HBV'

        if not os.path.exists(self.outdir):
            os.system(f"mkdir -p {self.outdir}")
    
    def bam2fq_run(self):
        bam2fq(self.bam, self.outdir, name = f'{self.name}_HBV_align')
        self.HBV_align_fq = f'{self.outdir}/{self.name}_HBV_align.fq'


    @utils.add_log
    def star_run(self):
        if(self.mod == "human"):
            #Fusion_bam.star_run.logger.info(f'{self.name} star begin.')
            star(self.HBV_align_fq, genome = "/SGRNJ06/randd/public/genome/rna/hs/hs_ensembl_99", outdir=self.outdir,  name = f'{self.fusion_human_name}')
            #Fusion_bam.star_run.logger.info(f'{self.name} star done.')
        if(self.mod == "homo_mus"):
            #Fusion_bam.star_run.logger.info(f'{self.name} star begin.')
            star(self.HBV_align_fq, genome = "/SGRNJ06/randd/public/genome/rna/hs/hs_ensembl_99", outdir=self.outdir,  name = f'{self.fusion_human_name}')
            star(self.HBV_align_fq, genome = "/SGRNJ06/randd/public/genome/rna/mmu/mmu_ensembl_99", outdir=self.outdir,  name = f'{self.fusion_mouse_name}')
            #Fusion_bam.star_run.logger.info(f'{self.name} star done.')

    def bam2fa_run(self):
        if(self.mod == "human"):
            bam2fq(f'{self.outdir}/{self.fusion_human_name}_Aligned.out.bam', outdir=self.outdir, name = f'{self.fusion_human_name}')
            fq2fa(f'{self.outdir}/{self.fusion_human_name}.fq', outdir=self.outdir, name = f'{self.fusion_human_name}')

        if(self.mod == "homo_mus"):
            # human
            bam2fq(f'{self.outdir}/{self.fusion_human_name}_Aligned.out.bam', outdir=self.outdir, name = f'{self.fusion_human_name}')
            fq2fa(f'{self.outdir}/{self.fusion_human_name}.fq', outdir=self.outdir, name = f'{self.fusion_human_name}')
            # mouse
            bam2fq(f'{self.outdir}/{self.fusion_mouse_name}_Aligned.out.bam', outdir=self.outdir, name = f'{self.fusion_mouse_name}')
            fq2fa(f'{self.outdir}/{self.fusion_mouse_name}.fq', outdir=self.outdir, name = f'{self.fusion_mouse_name}')
    
    def run(self):
        self.bam2fq_run()
        self.star_run()
        self.bam2fa_run()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--path', help='path', required=True)
    parsers.add_argument('--mod', help='mod', required=True)
    parsers.add_argument('--name', help='name', required=True)

    args = parsers.parse_args()
    runner = Fusion_bam(args) 
    runner.run()

if __name__ == '__main__':
    main()
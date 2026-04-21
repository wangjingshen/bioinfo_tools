from pathlib import Path
import os
import subprocess
import pysam
from collections import defaultdict
import pandas as pd
import glob
import argparse


def get_align_reads(bam, fq, outdir):
    cmd1 = f'samtools view {bam}|cut -f 1|sort|uniq > {outdir}/align.name'
    cmd2 = f'/SGRNJ/Public/Software/conda_env/cellrank/bin/seqtk seq -a {fq} > {outdir}/all.fa'
    cmd3 = (f"/SGRNJ/Public/Software/conda_env/cellrank/bin/seqtk subseq {outdir}/all.fa {outdir}/align.name > {outdir}/align.fa")
    subprocess.check_call(cmd1, shell=True)
    subprocess.check_call(cmd2, shell=True)
    subprocess.check_call(cmd3, shell=True)

def get_fusion_name(bam, mod):
    hbv_name = set()
    human_name = set()
    if(mod == "homo_mus"):
        mouse_name = set()

    if(mod == "human"):
        with pysam.AlignmentFile(bam) as bam, open(f'{outdir}/human_HBV_fusion_name','w') as fh:
            for i in bam:
                if(i.reference_name == "HBV_NC_003977"):
                    hbv_name.add(i.qname)
                else:
                    human_name.add(i.qname)
            huamn_hbv_fusion_name =  human_name & hbv_name
            for i in huamn_hbv_fusion_name:
                fh.write(f'{i}\n')
            #return(huamn_hbv_fusion_name)

    if(mod == "homo_mus"):
        with pysam.AlignmentFile(bam) as bam, open(f'{outdir}/human_HBV_fusion_name','w') as fh1, open(f'{outdir}/mouse_HBV_fusion_name','w') as fh2:
            for i in bam:
                if(i.reference_name == "HBV_NC_003977"):
                    hbv_name.add(i.qname)
                if(i.reference_name.startswith("homo")):
                    human_name.add(i.qname)
                if(i.reference_name.startswith("mus")):
                    mouse_name.add(i.qname)
            human_hbv_fusion_name =  human_name & hbv_name
            mouse_hbv_fusion_name =  mouse_name & hbv_name
            for i in human_hbv_fusion_name:
                fh1.write(f'{i}\n')
            for i in mouse_hbv_fusion_name:
                fh2.write(f'{i}\n')
            #return(human_hbv_fusion_name, mouse_hbv_fusion_name)


class Fusion_reads:
    def __init__(self, args):
        self.fq = glob.glob(f"{args.path}/02.cutadapt/*_clean_2.fq")[0]
        self.bam = glob.glob(f"{args.path}/03.star/*_Aligned.sortedByCoord.out.bam")[0]
        self.mod = args.mod
        self.outdir = args.outdir
        if not os.path.exists(self.outdir):
            os.system(f"mkdir -p {self.outdir}")
    
    def get_align_reads_run(self):
        get_align_reads(self.bam, self.fq, self.outdir)
        self.align_reads = f'{self.outdir}/align.fa'

    def get_fusion_name_run(self):
        if(self.mod == "human"):
            get_fusion_name(self.bam, self.mod)
            self.human_hbv_fusion_name = f'{self.outdir}/human_hbv_fusion_name'
        if(self.mod == "homo_mus"):
            get_fusion_name(self.bam, self.mod)
            self.human_hbv_fusion_name = f'{self.outdir}/human_hbv_fusion_name'
            self.mouse_hbv_fusion_name = f'{self.outdir}/mouse_hbv_fusion_name'
            
    def subseq_run(self):
        if(self.mod == "human"):
            cmd = (f"/SGRNJ/Public/Software/conda_env/cellrank/bin/seqtk subseq {self.align_reads} {self.human_hbv_fusion_name} > {outdir}/human_hbv_fusion.fa ")
            subprocess.check_call(cmd1, shell = True)

        if(self.mod == "homo_mus"):
            cmd1 = (f"/SGRNJ/Public/Software/conda_env/cellrank/bin/seqtk subseq {self.align_reads} {self.human_hbv_fusion_name} > {outdir}/human_hbv_fusion.fa ")
            cmd2 = (f"/SGRNJ/Public/Software/conda_env/cellrank/bin/seqtk subseq {self.align_reads} {self.mouse_hbv_fusion_name} > {outdir}/mouse_hbv_fusion.fa ")
            subprocess.check_call(cmd1, shell = True)
            subprocess.check_call(cmd2, shell = True)
    
    def run(self):
        self.get_align_reads_run()
        self.get_fusion_name_run()
        self.subseq_run()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--path', help='path', required=True)
    parsers.add_argument('--mod', help='mod', required=True)
    parsers.add_argument('--outdir', help='dir to save result file', required=True)

    args = parsers.parse_args()
    runner = Fusion_reads(args) 
    runner.run()

if __name__ == '__main__':
    main()
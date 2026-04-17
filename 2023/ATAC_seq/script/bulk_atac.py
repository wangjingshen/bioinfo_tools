import os
import pysam
import pandas as pd
import subprocess
import argparse
from collections import defaultdict
from itertools import zip_longest
import pyfastx

class Bulk_atac():
    '''
    use cr-atac run bulk atac(R1,R2)
    '''
    def __init__(self, fastq_path, prefix, species, outdir):
        self.fastq_path = fastq_path
        self.prefix = prefix
        self.species = species
        self.R2_10X_file = "/SGRNJ06/randd/USER/wangjingshen/script/ATAC_seq/data/atac_8_pbmc_1k_nextgem_S1_L001_R2_001.fastq.gz"
        self.outdir = outdir
    

    def cp_fastq(self):
        '''
        download R1:  10X input R1
        download R2:  10X input R3
        '''
        if not os.path.exists(self.prefix):
            os.system(f"mkdir -p {self.prefix}")
        
        cmd1 = f"cp {self.fastq_path}/{self.prefix}*R1*fastq.gz {self.prefix}/{self.prefix}_S1_L001_R1_001.fastq.gz"
        cmd2 = f"cp {self.fastq_path}/{self.prefix}*R2*fastq.gz {self.prefix}/{self.prefix}_S1_L001_R3_001.fastq.gz"
        #cmd1 = f"cp {self.fastq_path}/{self.prefix}_R1.fastq.gz {self.prefix}/{self.prefix}_S1_L001_R1_001.fastq.gz"
        #cmd2 = f"cp {self.fastq_path}/{self.prefix}_R2.fastq.gz {self.prefix}/{self.prefix}_S1_L001_R3_001.fastq.gz"
        cmd3 = f"zcat {self.prefix}/{self.prefix}_S1_L001_R3_001.fastq.gz | sed -e 's/2:N:0/3:N:0/g' > {self.prefix}/{self.prefix}_S1_L001_R3_001.fastq"
        cmd4 = f"gzip -f {self.prefix}/{self.prefix}_S1_L001_R3_001.fastq"        
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)


    def get_R2(self):
        '''
        get R2 from 10X R2
        '''
        #R1 = pyfastx.Fastq(f"{self.prefix}/{self.prefix}_S1_L001_R1_001.fastq.gz")
        #n_reads = len(R1)
        #print(n_reads)
        with pysam.FastxFile(f"{self.prefix}/{self.prefix}_S1_L001_R1_001.fastq.gz") as R2_analysis, pysam.FastxFile(f"{self.R2_10X_file}") as R2_10X, open(f"{self.prefix}/{self.prefix}_S1_L001_R2_001.fastq", mode='w') as fout:
            for (i_analysis,i_10X) in zip_longest(R2_analysis, R2_10X):
                if(i_10X == None):
                    #zcat /SGRNJ03/randd/user/wangjingshen/bio_soft/cellranger-atac-2.0.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz | head -n 1
                    i_analysis.sequence = "GGATAGGGTCATAGCT"      # set the extra as the first barcode in 737K-cratac-v1.txt.gz    
                    i_analysis.quality = i_analysis.quality[0:16]
                elif( i_analysis != None):
                    i_analysis.sequence = i_10X.sequence
                    i_analysis.quality = i_analysis.quality[0:16]
                else:
                    break
                fout.write(str(i_analysis) + '\n')

        cmd1 = f"sed -i 's/1:N:0/2:N:0/g' {self.prefix}/{self.prefix}_S1_L001_R2_001.fastq"
        cmd2 = f"gzip {self.prefix}/{self.prefix}_S1_L001_R2_001.fastq"
        #cmd1 = f"seqkit head -n {n_reads} {self.R2_10X_file} > {self.prefix}/{self.prefix}_S1_L001_R2_001.fastq"
        #cmd3 = f" rm {self.prefix}/{self.prefix}_S1_L001_R1_001.fastq.gz.fxi"

        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        #subprocess.check_call(cmd3, shell = True)
        
        
    def run_cratac(self):
        if(self.species == "human"):
            reference = "/SGRNJ06/randd/USER/wangjingshen/share/10x_atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/"
        if(self.species == "mouse"):
            reference = "/SGRNJ06/randd/USER/wangjingshen/share/10x_atac/refdata-cellranger-arc-mm10-2020-A-2.0.0/"
        if(self.species == "human_mouse"):
            reference = "/SGRNJ06/randd/USER/wangjingshen/share/10x_atac/refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0/"
        if(self.species == "mustela_putorius"):
            reference = "/SGRNJ06/randd/USER/cjj/ref/atac_mustela/mustela_putorius/"
        if(self.species == "mustela_putorius_v2"):
            reference = "/SGRNJ06/randd/USER/wangjingshen/share/genome/mustela_putorius_v2/mustela_putorius/"
        
        cmd = f"cellranger-atac count \
                    --id={self.outdir} \
                    --reference={reference} \
                    --fastqs={self.prefix} \
                    --localcores=8 --localmem=64"
        subprocess.check_call(cmd, shell = True)
    
    def run(self):
        self.cp_fastq()
        self.get_R2()
        self.run_cratac()

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--fastq_path', help='path of fastq file', required=True)
    parsers.add_argument('--prefix', help='sample name', required=True)
    parsers.add_argument('--species', help='species', required=True)
    parsers.add_argument('--outdir', help='dir of output', required=True)

    args = parsers.parse_args()
    runner = Bulk_atac(fastq_path, prefix, species, outdir) 
    runner.run()

if __name__ == '__main__':
    main()

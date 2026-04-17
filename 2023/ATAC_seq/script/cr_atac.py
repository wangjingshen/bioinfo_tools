import os
import pysam
import pandas as pd
import subprocess
import argparse
from collections import defaultdict
from itertools import zip_longest


class Atac():
    def __init__(self, args):
        self.fastq_path = args.fastq_path
        self.run_path = args.run_path
        #self.input_count = int(args.input_count)
        self.prefix = args.prefix
        self.species = args.species
        self.R2_10X = args.R2_10X
        self.outdir = args.outdir
    

    def cp_fastq(self):
        if not os.path.exists(self.run_path):
            os.system(f"mkdir -p {self.run_path}")
        
        cmd1 = f"cp {self.fastq_path}/{self.prefix}_R1.fastq.gz {self.run_path}/{self.prefix}_S1_L001_R1_001.fastq.gz"
        cmd2 = f"cp {self.fastq_path}/{self.prefix}_R2.fastq.gz {self.run_path}/{self.prefix}_R2.fastq.gz"
        cmd3 = f"cp {self.fastq_path}/{self.prefix}_R3.fastq.gz {self.run_path}/{self.prefix}_S1_L001_R3_001.fastq.gz"
        
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)

        I1_file = f"{self.fastq_path}/{self.prefix}_I1.fastq.gz"
        if(os.path.isfile(I1_file)):
            print("I1 exist!")
            cmd4 = f"cp {self.fastq_path}/{self.prefix}_I1.fastq.gz {self.run_path}/{self.prefix}_S1_L001_I1_001.fastq.gz"
            subprocess.check_call(cmd4, shell = True)
        else:
            print("I1 not exist!")



    def process_R2(self):
        cmd1 = f"zcat {self.run_path}/{self.prefix}_R2.fastq.gz | sed -n '$=' > {self.run_path}/tmp"
        subprocess.check_call(cmd1, shell = True)
        n_lines = int(pd.read_csv(f"{self.run_path}/tmp").columns[0])
        print(n_lines)

        ## pbmc
        if( self.R2_10X == "pbmc"):
            self.R2_10X_file = "/SGRNJ06/randd/USER/wangjingshen/project/10X_ATAC/data/10X_barcode_R2_file/atac_2_pbmc_1k_nextgem_S1_L001_R2_001.fastq.gz"
            with pysam.FastxFile(f"{self.run_path}/{self.prefix}_R2.fastq.gz") as R2_analysis, pysam.FastxFile(f"{self.R2_10X_file}") as R2_10X, open(f"{self.run_path}/{self.prefix}_S1_L001_R2_001.fastq", mode='w') as fout:
                if(n_lines > 209972328):
                    for (i_analysis,i_10X) in zip_longest(R2_analysis, R2_10X):
                        if(i_10X == None):
                            #zcat /SGRNJ03/randd/user/wangjingshen/bio_soft/cellranger-atac-2.0.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz | head -n 1
                            i_analysis.sequence = "GGATAGGGTCATAGCT"      # set the extra as the first barcode in 737K-cratac-v1.txt.gz    
                        else:
                            i_analysis.sequence = i_10X.sequence
                        i_analysis.quality = i_analysis.quality[0:16]
                        fout.write(str(i_analysis) + '\n')

                else:
                    for (i_analysis,i_10X) in zip(R2_analysis, R2_10X):
                        i_analysis.sequence = i_10X.sequence
                        i_analysis.quality = i_analysis.quality[0:16] 
                        fout.write(str(i_analysis) + '\n')
            
            cmd2 = f"rm {self.run_path}/tmp "        
            cmd3 = f" gzip {self.run_path}/{self.prefix}_S1_L001_R2_001.fastq"
            cmd4 = f" rm {self.run_path}/{self.prefix}_R2.fastq.gz"
            subprocess.check_call(cmd2, shell = True)
            subprocess.check_call(cmd3, shell = True)
            subprocess.check_call(cmd4, shell = True)
        
        ## tiny
        if( self.R2_10X == "tiny"):
            self.R2_10X_file = "/SGRNJ06/randd/USER/wangjingshen/project/10X_ATAC/data/10X_barcode_R2_file/tiny_2_S0_L001_R2_001.fastq.gz"
            with pysam.FastxFile(f"{self.run_path}/{self.prefix}_R2.fastq.gz") as R2_analysis, pysam.FastxFile(f"{self.run_path}/10X_R2.fastq") as R2_10X, open(f"{self.run_path}/{self.prefix}_S1_L001_R2_001.fastq", mode='w') as fout:
                for (i_analysis,i_10X) in zip(R2_analysis, R2_10X):
                    i_analysis.sequence = i_10X.sequence
                    i_analysis.quality = i_analysis.quality[0:16]
                    fout.write(str(i_analysis) + '\n')
            
            cmd2 = f"rm {self.run_path}/tmp "        
            cmd3 = f" gzip {self.run_path}/{self.prefix}_S1_L001_R2_001.fastq"
            cmd4 = f" rm {self.run_path}/{self.prefix}_R2.fastq.gz"
            subprocess.check_call(cmd2, shell = True)
            subprocess.check_call(cmd3, shell = True)
            subprocess.check_call(cmd4, shell = True)


        
        ## sc, need test
        if( self.R2_10X == "sc"):
            barcode_whitelist = pd.read_csv("/SGRNJ03/randd/user/wangjingshen/bio_soft/cellranger-atac-2.0.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz",header=None)
            barcode_whitelist.columns = ["barcode"]
            ref_dict = {x:y for x,y in zip((set(record_dict_df.sequence)),set(barcode_whitelist.barcode))}
            with pysam.FastxFile(f"{self.run_path}/{self.prefix}_R2.fastq.gz") as R2_analysis, open(f"{self.run_path}/{self.prefix}_S1_L001_R2_001.fastq", mode='w') as fout:
                for i_analysis in R2_analysis:
                    i_analysis.sequence = ref_dict[i_analysis.sequence]
                    i_analysis.quality = i_analysis.quality[0:16]
                    fout.write(str(i_analysis) + '\n')

            cmd2 = f"rm {self.run_path}/tmp "        
            cmd3 = f" gzip {self.run_path}/{self.prefix}_S1_L001_R2_001.fastq"
            cmd4 = f" rm {self.run_path}/{self.prefix}_R2.fastq.gz"
            subprocess.check_call(cmd2, shell = True)
            subprocess.check_call(cmd3, shell = True)
            subprocess.check_call(cmd4, shell = True)

        # 10X   
        if( self.R2_10X == "10X"):
            cmd1 = f"mv {self.run_path}/{self.prefix}_R2.fastq.gz {self.run_path}/{self.prefix}_S1_L001_R2_001.fastq.gz"
            subprocess.check_call(cmd1, shell = True)
            
        
    def run_cratac(self):
        if(self.species == "human"):
            reference = "/SGRNJ06/randd/USER/wangjingshen/share/10x_atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/"
        if(self.species == "mouse"):
            reference = "/SGRNJ06/randd/USER/wangjingshen/share/10x_atac/refdata-cellranger-arc-mm10-2020-A-2.0.0/"
        if(self.species == "human_mouse"):
            reference = "/SGRNJ06/randd/USER/wangjingshen/share/10x_atac/refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0/"
        
        cmd = f"cellranger-atac count \
                    --id={self.outdir} \
                    --reference={reference} \
                    --fastqs={self.run_path} \
                    --localcores=8 --localmem=64"
        subprocess.check_call(cmd, shell = True)
    
    def run(self):
        self.cp_fastq()
        self.process_R2()
        self.run_cratac()

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--fastq_path', help='path of fastq file', required=True)
    parsers.add_argument('--run_path', help='run path', required=True)
    parsers.add_argument('--prefix', help='sample name', required=True)
    parsers.add_argument('--species', help='species', required=True)
    parsers.add_argument('--R2_10X', help='R2 from 10X, pbmc(large) or tiny(small)', required=True)
    parsers.add_argument('--outdir', help='dir of output', required=True)

    args = parsers.parse_args()
    runner = Atac(args) 
    runner.run()

if __name__ == '__main__':
    main()

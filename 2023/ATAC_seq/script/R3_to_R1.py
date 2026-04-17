import os
import pysam
import pandas as pd
import subprocess
import argparse
from collections import defaultdict
import glob
#from generate_mapfile_function import make_mapfile

class R3_to_R1():
    '''
    '''
    def __init__(self, args):
        self.I1 = glob.glob(os.path.join(args.fastq_path,"*I1.fastq.gz"))[0]
        self.R1 = glob.glob(os.path.join(args.fastq_path,"*R1.fastq.gz"))[0]
        self.R2 = glob.glob(os.path.join(args.fastq_path,"*R2.fastq.gz"))[0]
        self.R3 = glob.glob(os.path.join(args.fastq_path,"*R3.fastq.gz"))[0]
        self.prefix = args.prefix

        if not os.path.exists("_tmp/"):
            os.system(f"mkdir -p _tmp")

        self.replace_path = f"_tmp/{self.prefix}_replace"

        if not os.path.exists(f"{self.replace_path}"):
            os.system(f"mkdir -p {self.replace_path}")


    def replace(self):
        record_dict = defaultdict(list)
        with pysam.FastxFile(f"{self.R1}") as R1, pysam.FastxFile(f"{self.R3}") as R3,\
            open(f"{self.replace_path}/{self.prefix}_R1.fastq", mode='w') as R1_out:
            for (R1_i, R3_i) in zip(R1, R3):
                R1_i.sequence = R3_i.sequence
                R1_i.quality = R3_i.quality
                R1_out.write(str(R1_i) + '\n')
        
        cmd1 = f" gzip {self.replace_path}/{self.prefix}_R1.fastq"
        cmd2 = f" cp {self.I1} {self.replace_path}/"
        cmd3 = f" cp {self.R2} {self.replace_path}/"
        cmd4 = f" cp {self.R3} {self.replace_path}/"
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)
    
    def run(self):
        self.replace()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--fastq_path', help='path of fastq file', required=True)
    parsers.add_argument('--prefix', help='sample name', required=True)

    args = parsers.parse_args()
    runner = R3_to_R1(args) 
    runner.run()

if __name__ == '__main__':
    main()

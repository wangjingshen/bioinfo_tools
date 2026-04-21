import os
import pysam
import pandas as pd
import subprocess
import argparse
from collections import defaultdict
import glob
#from generate_mapfile_function import make_mapfile

class Filter_poor_reads():
    '''
    insert + filter
    '''
    def __init__(self, args):
        self.I1 = glob.glob(os.path.join(args.fastq_path,"*I1.fastq.gz"))[0]
        self.R1 = glob.glob(os.path.join(args.fastq_path,"*R1.fastq.gz"))[0]
        self.R2 = glob.glob(os.path.join(args.fastq_path,"*R2.fastq.gz"))[0]
        self.R3 = glob.glob(os.path.join(args.fastq_path,"*R3.fastq.gz"))[0]
        self.prefix = args.prefix

        if not os.path.exists("_tmp/"):
            os.system(f"mkdir -p _tmp")

        self.cutadapt_path = f"_tmp/{self.prefix}_cutadat"
        self.filter_path = f"_tmp/{self.prefix}_filter"

        if not os.path.exists(f"{self.cutadapt_path}"):
            os.system(f"mkdir -p {self.cutadapt_path}")
        if not os.path.exists(f"{self.filter_path}"):
            os.system(f"mkdir -p {self.filter_path}")


    def run_cutadapt(self):
        '''
        # R1 adapter R2 sequence:  CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
        # R2 adapter R1 sequence:  CTGTCTCTTATACACATCTGACGCTGCCGACGA
        '''
        cmd1 = f'cutadapt -j 8 --minimum-length 0 -a "GGGGGGGGGG" -a "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC" -o {self.cutadapt_path}/{self.prefix}_clean_R1.fastq.gz {self.R1}'
        cmd2 = f'cutadapt -j 8 --minimum-length 0 -a "GGGGGGGGGG" -a "CTGTCTCTTATACACATCTGACGCTGCCGACGA" -o {self.cutadapt_path}/{self.prefix}_clean_R3.fastq.gz {self.R3}'
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)


    def run_get_insert(self):
        def get_insert(name):
            record_dict = defaultdict(list)
            with pysam.FastxFile(f"{self.cutadapt_path}/{self.prefix}_clean_{name}.fastq.gz") as f:
                for i in f:
                    record_dict["name"].append(i.name)
                    record_dict["sequence"].append(i.sequence)
            record_df = pd.DataFrame(record_dict)
            record_df['length'] = record_df.sequence.apply(len)

            sdict = {}
            sdict.update({
                "sample_prefix": self.prefix,
                "total_reads": len(record_df),
                "mean_insert_size": record_df['length'].mean(),
                "median_insert_size": record_df['length'].median()
            })
        
            with open(f"{self.cutadapt_path}/{self.prefix}_{name}_insert_size_summarize.txt", "w") as fd:
                for k, v in sdict.items():
                    fd.write(f"{k}: {v}\n")

            #df_res = f"{self.outdir}/{self.prefix}_insert_size_df.tsv"
            #record_df.to_csv(df_res, sep="\t")
        get_insert("R1")
        get_insert("R3")   


    def filter_poor_reads(self):
        record_dict = defaultdict(list)
        with pysam.FastxFile(f"{self.I1}") as I1, pysam.FastxFile(f"{self.R1}") as R1, pysam.FastxFile(f"{self.R2}") as R2,pysam.FastxFile(f"{self.R3}") as R3,\
            pysam.FastxFile(f"{self.cutadapt_path}/{self.prefix}_clean_R1.fastq.gz") as R1_cut, pysam.FastxFile(f"{self.cutadapt_path}/{self.prefix}_clean_R3.fastq.gz") as R3_cut,\
            open(f"{self.filter_path}/{self.prefix}_I1.fastq", mode='w') as I1_out, open(f"{self.filter_path}/{self.prefix}_R1.fastq", mode='w') as R1_out,\
            open(f"{self.filter_path}/{self.prefix}_R2.fastq", mode='w') as R2_out, open(f"{self.filter_path}/{self.prefix}_R3.fastq", mode='w') as R3_out:

            for (I1_i, R1_i, R2_i, R3_i, R1_cut_i, R3_cut_i) in zip(I1, R1, R2, R3, R1_cut, R3_cut):
                if((len(R1_cut_i.sequence) >= 30) & (len(R3_cut_i.sequence) >= 30)):
                    I1_out.write(str(I1_i) + '\n')
                    R1_out.write(str(R1_i) + '\n')
                    R2_out.write(str(R2_i) + '\n')
                    R3_out.write(str(R3_i) + '\n')
        
        cmd1 = f" gzip {self.filter_path}/{self.prefix}_I1.fastq"
        cmd2 = f" gzip {self.filter_path}/{self.prefix}_R1.fastq"
        cmd3 = f" gzip {self.filter_path}/{self.prefix}_R2.fastq"
        cmd4 = f" gzip {self.filter_path}/{self.prefix}_R3.fastq"
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)

    def run(self):
        self.run_cutadapt()
        self.run_get_insert()
        self.filter_poor_reads()

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--fastq_path', help='path of fastq file', required=True)
    parsers.add_argument('--prefix', help='sample name', required=True)

    args = parsers.parse_args()
    runner = Filter_poor_reads(args) 
    runner.run()

if __name__ == '__main__':
    main()

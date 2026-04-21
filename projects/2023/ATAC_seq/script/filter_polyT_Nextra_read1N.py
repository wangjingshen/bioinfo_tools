import os
import pysam
import pandas as pd
import subprocess
import argparse
from collections import defaultdict
import glob
#from generate_mapfile_function import make_mapfile

class Filter_polyT_Nextra_read1N():
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
        # polyT:  CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
        # Nextra read1N:  TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
        # Retain the following sequences(-g)
        '''
        cmd = f'cutadapt -j 8 --minimum-length 0 -g "TTTTTTTTTT" -g "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" -o {self.cutadapt_path}/{self.prefix}_clean_R1.fastq.gz {self.R1}'
        subprocess.check_call(cmd, shell = True)


    def filter_poor_reads(self):
        record_dict = defaultdict(list)
        with pysam.FastxFile(f"{self.I1}") as I1, pysam.FastxFile(f"{self.R1}") as R1, pysam.FastxFile(f"{self.R2}") as R2,pysam.FastxFile(f"{self.R3}") as R3,\
            pysam.FastxFile(f"{self.cutadapt_path}/{self.prefix}_clean_R1.fastq.gz") as R1_cut,\
            open(f"{self.filter_path}/{self.prefix}_I1.fastq", mode='w') as I1_out, open(f"{self.filter_path}/{self.prefix}_R1.fastq", mode='w') as R1_out,\
            open(f"{self.filter_path}/{self.prefix}_R2.fastq", mode='w') as R2_out, open(f"{self.filter_path}/{self.prefix}_R3.fastq", mode='w') as R3_out:

            for (I1_i, R1_i, R2_i, R3_i, R1_cut_i) in zip(I1, R1, R2, R3, R1_cut):
                if(len(R1_cut_i.sequence) >= 141):
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
    
    def get_summary(self):
        # raw R1
        record_dict1 = defaultdict(list)
        with pysam.FastxFile(self.R1) as f:
            for i in f:
                record_dict1["name"].append(i.name)
        record_df1 = pd.DataFrame(record_dict1)

        # update R1
        record_dict2 = defaultdict(list)
        with pysam.FastxFile(f"{self.filter_path}/{self.prefix}_R1.fastq.gz") as f:
            for i in f:
                record_dict2["name"].append(i.name)
        record_df2 = pd.DataFrame(record_dict2)

        sdict = {}
        sdict.update({
            "sample_prefix": self.prefix,
            "raw_reads": len(record_df1),
            "filter_reads": len(record_df2),
            "valid percent": str(round(100*len(record_df2)/len(record_df1),3))+" %"
        })
        
        with open(f"{self.cutadapt_path}/{self.prefix}_summary.txt", "w") as fd:
            for k, v in sdict.items():
                fd.write(f"{k}: {v}\n")

        #df_res = f"{self.outdir}/{self.prefix}_insert_size_df.tsv"
        #record_df.to_csv(df_res, sep="\t")

    def run(self):
        self.run_cutadapt()
        self.filter_poor_reads()
        self.get_summary()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--fastq_path', help='path of fastq file', required=True)
    parsers.add_argument('--prefix', help='sample name', required=True)

    args = parsers.parse_args()
    runner = Filter_polyT_Nextra_read1N(args) 
    runner.run()

if __name__ == '__main__':
    main()

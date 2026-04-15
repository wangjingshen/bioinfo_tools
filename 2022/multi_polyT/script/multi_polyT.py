import pysam
import pandas as pd
from collections import defaultdict
import pyranges as pr
from functools import reduce
import glob
import os
import argparse
import subprocess
#import faulthandler
#faulthandler.enable()


class Multi_polyT:
    def __init__(self, args):
        self.mapfile = pd.read_csv(f"{args.mapfile}", sep="\t", header=None)
        self.R1 = self.mapfile.iloc[0,1] + "/" + self.mapfile.iloc[0,0] + "_R1.fastq.gz"
        self.R2 = self.mapfile.iloc[0,1] + "/" + self.mapfile.iloc[0,0] + "_R2.fastq.gz"      
    
    def get_multi_polyT_R1(self):
        '''
        --no-group-separator  # for rm --
        '''
        cmd1 = f"mkdir -p raw_fastq/"
        cmd2 = f"zcat {self.R1} | grep -E -A 2 -B 1 --no-group-separator 'TTTTTTTTTT[ACG]+[ATCG]*TTTTTTTTTT' > raw_fastq/multi_polyT_R1.fastq"
        cmd3 = f"gzip raw_fastq/multi_polyT_R1.fastq"
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)

    def get_multi_polyT_R2(self):
        ref_record_dict1 = defaultdict(list)
        with pysam.FastxFile(f"raw_fastq/multi_polyT_R1.fastq.gz") as R1:
            for i in R1:
                ref_record_dict1['name'].append(i.name)
                ref_record_dict1['name1'].append(i.name)
        ref_record_df1 = pd.DataFrame(ref_record_dict1)

        ref_record_dict2 = defaultdict(list)
        with pysam.FastxFile(f"{self.R2}") as R2, open(f"raw_fastq/multi_polyT_R2.fastq", mode='w') as out:
            for i in R2:
                ref_record_dict2['name'].append(i.name)
                ref_record_dict2['record'].append(i)
            ref_record_df2 = pd.DataFrame(ref_record_dict2)
            ref_record_df2_filter = ref_record_df2.loc[ref_record_df2.name.isin(ref_record_df1.name),]

            for _, row in ref_record_df2_filter.iterrows():
                out.write(str(row["record"]) + '\n')  

        cmd1 = f"gzip raw_fastq/multi_polyT_R2.fastq"
        subprocess.check_call(cmd1, shell = True)

    def make_mapfile(self):
        mapfile = self.mapfile
        mapfile.iloc[0,0] = "multi_polyT"
        mapfile.iloc[0,1] = "raw_fastq/"
        mapfile.iloc[0,2] = mapfile.iloc[0,2]+"_multi_polyT"
        mapfile.to_csv("multi_polyT.mapfile", sep="\t", header=None, index=None)

    def run_celescope(self):
        cmd1 = f"/SGRNJ/Public/Software/conda_env/celescope1.5.1b0/bin/multi_capture_virus \
            --mapfile multi_polyT.mapfile \
            --virus_genome /SGRNJ03/randd/test_rd/dxh/data/HBV_genome/ \
            --thread 4\
            --not_consensus \
            --outFilterMatchNmin 80 \
            --umi_threshold otsu \
            --allowNoPolyT \
            --gtf /SGRNJ06/randd/USER/wangjingshen/project/huashan_HBV/data/2022-10-18_gtf/HBV.gtf"
        subprocess.check_call(cmd1, shell = True)
    
    def run(self):
        self.get_multi_polyT_R1()
        self.get_multi_polyT_R2()
        self.make_mapfile()
        self.run_celescope()
    

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--mapfile', help='mapfile', required=True)

    args = parsers.parse_args()
    runner = Multi_polyT(args) 
    runner.run()

if __name__ == '__main__':
    main()
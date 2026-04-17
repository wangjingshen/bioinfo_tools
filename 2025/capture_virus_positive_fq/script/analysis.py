import collections
import subprocess
import pandas as pd
import pysam
import argparse
import os


class FeatureCounts_capture_virus():
    def __init__(self, args):
        self.args = args
        self.bam = args.bam
        self.filter_umi_file = args.filter_umi_file 
        self.filter_read_count_json = args.filter_read_count_json
        self.outdir = args.outdir
        self.sample = args.sample

        self.filter_bam = f'{self.outdir}/{self.sample}_filter.bam'
        self.filter_fq = f'{self.outdir}/{self.sample}_filter.fq'

        if not os.path.exists(f'{self.outdir}'):
            os.system(f'mkdir -p {self.outdir}')
    
    def get_valid_barcodes(filter_umi_file):     ## UMI min support reads.    barcode UMI not otsu
        df = pd.read_table(filter_umi_file)
        df = df[df['UMI'] > 0]
        valid_barcodes = set(df['barcode'])
        return valid_barcodes    
    
    def get_valid_umis(filter_read_count_df, valid_barcodes):
        df = pd.read_table(filter_read_count_df)
        filter_bc_df = df[df["barcode"].isin(valid_barcodes)][["barcode","UMI"]].drop_duplicates()
        barcode_umis = collections.defaultdict(set)
        all_bc = filter_bc_df["barcode"].to_list()
        for bc in all_bc:
            umi_set = set(filter_bc_df[filter_bc_df["barcode"]==bc]["UMI"].to_list())
            barcode_umis[bc] = barcode_umis[bc] | umi_set
        return barcode_umis

    def run_filter(self):
        valid_barcodes = FeatureCounts_capture_virus.get_valid_barcodes(self.args.filter_umi_file)
        barcode_umis = FeatureCounts_capture_virus.get_valid_umis(self.args.filter_read_count_json, valid_barcodes)

        with pysam.AlignmentFile(self.args.bam, "rb") as raw_bam:
            header = raw_bam.header
            with pysam.AlignmentFile(self.filter_bam, "wb", header=header) as filter_bam:
                for read in raw_bam:
                    attr = read.query_name.split('_')
                    barcode = attr[0]
                    umi = attr[1]
                    if barcode in valid_barcodes and umi in barcode_umis[barcode]:
                        filter_bam.write(read)

    def bam2fq(self):
        cmd1 = f'samtools index {self.filter_bam}'
        cmd2 = f'bedtools bamtofastq -i {self.filter_bam} -fq {self.filter_fq}'    
        cmd3 = f'gzip {self.filter_fq}'    
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)


    def run(self):
        self.run_filter()
        self.bam2fq()

    
def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--bam', help='input bam file', required=True)
    parsers.add_argument('--filter_umi_file', help='filter umi file', required=True)
    parsers.add_argument('--filter_read_count_json', help='filter_read_count_json file', required=True)
    parsers.add_argument('--outdir', help='outdir', required=True)
    parsers.add_argument('--sample', help='sample', required=True)

    args = parsers.parse_args()
    runner = FeatureCounts_capture_virus(args) 
    runner.run()

if __name__ == '__main__':
    main()

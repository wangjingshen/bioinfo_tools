import subprocess
import pandas as pd
import pysam
import argparse
import os
import sys
from pathlib import Path
from typing import Dict, Set
import time
import json

dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from utils.utils import mkdir, logger, execute_cmd

class VirusPositiveFastq():
    def __init__(self, bam:str, filter_umi_file:str, filter_read_count_json:str, outdir:str, sample:str, threads:int):
        self.bam = Path(bam)
        self.filter_umi_file = Path(filter_umi_file) 
        self.filter_read_count_json = Path(filter_read_count_json)
        self.outdir = Path(outdir)
        mkdir(self.outdir)
        self.sample = sample
        self.threads = threads
        self.out_bam = f'{self.outdir}/{self.sample}_filter.bam'
        self.out_fq = f'{self.outdir}/{self.sample}_filter.fq'

        self._valid_barcodes: Set[str] | None = None
        self._valid_barcode_umis: Set[str] | None = None

    def _get_valid_barcodes(self) -> Set[str]:
        df = pd.read_csv(self.filter_umi_file, usecols=["barcode", "sum_UMI"])
        return set(df.loc[df['sum_UMI'] > 0, 'barcode'])   
    
    def _get_valid_barcode_umis(self, valid_barcodes: Set[str]) -> Dict[str, Set[str]]:
        data = json.load(self.filter_read_count_json.open())
        records = set()   # records = [] slow
        for barcode, umi_dict in data.items():
            if barcode not in valid_barcodes:
                continue
            for ref, umi_counts in umi_dict.items():
                for umi, count in umi_counts.items():
                    if count > 0:
                        records.add(f'{barcode}_{umi}')   # records.append(f'{barcode}_{umi}')  slow
        return records

    def filter_bam(self):
        logger.info("Running filter bam...")
        if self._valid_barcodes is None:
            self._valid_barcodes = self._get_valid_barcodes()
        if self._valid_barcode_umis is None:
            self._valid_barcode_umis = self._get_valid_barcode_umis(self._valid_barcodes)

        with pysam.AlignmentFile(self.bam, "rb") as raw_bam:
            with pysam.AlignmentFile(self.out_bam, "wb", header = raw_bam.header) as out_bam:
                for read in raw_bam:
                    try:
                        barcode, umi, _ = read.query_name.split('_', 2)
                    except ValueError:
                        continue

                    if f'{barcode}_{umi}' in self._valid_barcode_umis:
                        out_bam.write(read)

    def bam2fq(self) -> None:
        logger.info("Running bam2fq...")
        cmds = [
            f'samtools sort -@ {self.threads} -o {self.out_bam}.sorted {self.out_bam}',
            f'samtools index {self.out_bam}.sorted',
            f'bedtools bamtofastq -i {self.out_bam}.sorted -fq {self.out_fq}',    
            f'pigz -p {self.threads} {self.out_fq}'  #gzip {self.out_fq}
        ]
        for cmd in cmds:
            execute_cmd(cmd)

    def run(self) -> None:
        t0 = time.time()
        self.filter_bam()
        self.bam2fq()
        logger.info(f"Done.({time.time() - t0:.1f}s)")
    

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--bam', help='input bam file', required=True)
    parsers.add_argument('--filter_umi_file', help='filter umi file', required=True)
    parsers.add_argument('--filter_read_count_json', help='filter_read_count_json file', required=True)
    parsers.add_argument('--outdir', help='outdir', required=True)
    parsers.add_argument('--sample', help='sample', required=True)
    parsers.add_argument('--threads', default = 4, help='threads')

    args = parsers.parse_args()
    runner = VirusPositiveFastq(args.bam, args.filter_umi_file, args.filter_read_count_json, args.outdir, args.sample, args.threads) 
    runner.run()

if __name__ == '__main__':
    main()

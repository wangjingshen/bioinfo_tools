import os
import subprocess
import argparse
import gzip
import pandas as pd
import glob


def file_len(fname):
    with gzip.open(fname, 'rb') as f:
        for i,l in enumerate(f):
            pass
    return( (i+1)/4 )

class Check_data():
    def __init__(self, args):
        self.download_path = args.download_path
        self.raw_path = args.raw_path
        self.sample_name = args.sample_name
        self.outdir = args.outdir
        if not os.path.exists(f"{self.outdir}"):
            os.system(f"mkdir -p {self.outdir}")

    def compare_reads(self):
        if(self.download_path.find(',') != -1):
            download_reads = 0
            for i in self.download_path.split(','):
                #print(i)
                download_reads += file_len(glob.glob(f"{i}/*R1.fastq.gz")[0])
        else:
            download_reads = file_len(glob.glob(f"{self.download_path}/*R1.fastq.gz")[0])
        
        raw = pd.read_csv(f"{self.raw_path}/01.barcode/stat.txt", header=None, sep=":")
        raw_reads = int(raw.iloc[0,1].replace(',',''))   #  raw_lines = 4 * int(raw.iloc[0,1].replace(',',''))

        print(f"{self.sample_name}")
        print(f"reads of download file: {download_reads}")
        print(f"reads of raw file: {raw_reads}")
        if(download_reads == raw_reads):
            check_res = f"{self.sample_name}: reads are equal"
        else:
            check_res = f"{self.sample_name}: reads are not equal"
        print(check_res)
        print('------')
        print('')


        # save res --
        sdict = {}
        sdict.update({
            "sample name": self.sample_name,
            "raw_path": self.raw_path,
            "download_path": self.download_path,
            "reads of raw file": raw_reads,
            "reads of download file": download_reads,
            "check res": check_res
            })
        with open(f"{self.outdir}/{self.sample_name}_summarize.txt", "w") as fd:
            for k, v in sdict.items():
                fd.write(f"{k}: {v}\n")


    def run(self):
        self.compare_reads()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--download_path', help='download path', required=True)
    parsers.add_argument('--raw_path', help='raw_path')
    parsers.add_argument('--sample_name', help='sample_name')
    parsers.add_argument('--outdir', help='outdir')

    args = parsers.parse_args()
    runner = Check_data(args) 
    runner.run()

if __name__ == '__main__':
    main()
        
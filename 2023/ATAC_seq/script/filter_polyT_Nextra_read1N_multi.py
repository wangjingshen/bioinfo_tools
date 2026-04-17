import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd


class Filter_polyT_Nextra_read1N_r():
    def __init__(self, fastq_path, prefix):
        self.fastq_path = fastq_path
        self.prefix = prefix

    def run_cmd(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/ATAC_seq/filter_polyT_Nextra_read1N.py '
            f'--fastq_path {self.fastq_path} '
            f'--prefix {self.prefix} '
        )
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_cmd()


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\t',header=None)
    fastq_path_list = df_mapfile.iloc[:,1]
    prefix_list = df_mapfile.iloc[:,0]
    return fastq_path_list, prefix_list


def run_single(fastq_path, prefix):
    runner = Filter_polyT_Nextra_read1N_r(fastq_path, prefix)
    runner.run()


def main():
    mapfile = sys.argv[1]
    species = sys.argv[2]
    R2_10X = sys.argv[3]
    fastq_path_list, prefix_list = parse_mapfile(mapfile)
    print(fastq_path_list)
    print(prefix_list)
    with ProcessPoolExecutor(max_workers=4) as executor:
        for result in executor.map(run_single, fastq_path_list, prefix_list):
            print(result, 'done')
    
    # make mapfile
    raw_mapfile = pd.read_csv(mapfile,sep="\t",header=None)
    new_mapfile = raw_mapfile.copy()
    new_mapfile.columns = ["prefix", "fastq_path", "run_path"]
    new_mapfile['fastq_path'] = "_tmp/" + new_mapfile['prefix']+ "_filter"
    new_mapfile['species'] = species
    new_mapfile['outdir'] = new_mapfile["run_path"]+"_outdir"
    new_mapfile['R2_10X'] = R2_10X
    new_mapfile['run_path'] = "_tmp/" + new_mapfile['run_path']

    new_mapfile.to_csv("mapfile",sep="\t",index=False)

if __name__ == '__main__':
    main()
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Hbv_virus_gene_consensus_run():
    def __init__(self, fj_path, max_n_reads, small_threshold, large_threshold, min_consensus_read, outdir):
        self.fj_path = fj_path
        self.max_n_reads = max_n_reads
        self.small_threshold = small_threshold
        self.large_threshold = large_threshold
        self.min_consensus_read = min_consensus_read
        self.outdir = outdir
    
    def run(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/2023/HBV_virus_gene_consensus/analysis.py '
            f'--fj_path {self.fj_path} '
            f'--max_n_reads {self.max_n_reads} '
            f'--small_threshold {self.small_threshold} '
            f'--large_threshold {self.large_threshold} '
            f'--min_consensus_read {self.min_consensus_read} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    fj_path_list = df_mapfile['fj_path']
    max_n_reads_list = df_mapfile['max_n_reads']
    small_threshold_list = df_mapfile['small_threshold']
    large_threshold_list = df_mapfile['large_threshold']
    min_consensus_read_list = df_mapfile['min_consensus_read']
    outdir_list = df_mapfile['outdir']
    return fj_path_list, max_n_reads_list, small_threshold_list, large_threshold_list, min_consensus_read_list, outdir_list

def run_single(fj_path, max_n_reads, small_threshold, large_threshold, min_consensus_read, outdir):
    runner = Hbv_virus_gene_consensus_run(fj_path, max_n_reads, small_threshold, large_threshold, min_consensus_read, outdir)
    runner.run()

def main():
    mapfile = sys.argv[1]
    fj_path_list, max_n_reads_list, small_threshold_list, large_threshold_list, min_consensus_read_list, outdir_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, fj_path_list, max_n_reads_list, small_threshold_list, large_threshold_list, min_consensus_read_list, outdir_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from bulk_atac import Bulk_atac


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\t')
    fastq_path_list = df_mapfile['fastq_path']
    prefix_list = df_mapfile['prefix']
    species_list = df_mapfile['species']
    outdir_list = df_mapfile['outdir']
    
    return fastq_path_list, prefix_list, species_list, outdir_list


def run_single(fastq_path, prefix, species, outdir):
    runner = Bulk_atac(fastq_path, prefix, species, outdir)
    runner.run()


def main():
    mapfile = sys.argv[1]
    fastq_path_list, prefix_list, species_list, outdir_list = parse_mapfile(mapfile)
    if(len(outdir_list.unique()) == len(outdir_list)):
        with ProcessPoolExecutor(max_workers=11) as executor:
            executor.map(run_single, fastq_path_list, prefix_list, species_list, outdir_list)
    else:
        print("sample duplicate, exit. \n \
               please run /SGRNJ06/randd/USER/wangjingshen/script_dev/merge_fq_ATAC/test/merge.sh \n \
               then use merge_mapfile run bulk_atac_multi.")
if __name__ == '__main__':
    main()
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from space_tag import SpaceTag


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\t')
    space_dir_list = df_mapfile['space_dir']
    tag_dir_list = df_mapfile['tag_dir']
    sample_list = df_mapfile['sample']
    
    return space_dir_list, tag_dir_list, sample_list


def run_single(space_dir, tag_dir, sample):
    runner = SpaceTag(space_dir, tag_dir, sample)
    runner.run()


def main():
    mapfile = sys.argv[1]
    space_dir_list, tag_umi_list, sample_list = parse_mapfile(mapfile)
    with ProcessPoolExecutor(max_workers=1) as executor:
        executor.map(run_single, space_dir_list, tag_umi_list, sample_list)

if __name__ == '__main__':
    main()
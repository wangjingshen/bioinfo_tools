import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Top_reads_umi_stat_run():
    def __init__(self, cele_path, species, name):
        self.cele_path = cele_path
        self.species = species
        self.name = name
    
    def run_top_reads_umi_stat(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/top_reads_umi_stat/script/top_reads_umi_stat.py '
            f'--cele_path {self.cele_path} '
            f'--species {self.species} '
            f'--name {self.name} '
        )
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_top_reads_umi_stat()

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    cele_path_list = df_mapfile['cele_path']
    species_list = df_mapfile['species']
    name_list = df_mapfile['name']

    return cele_path_list, species_list, name_list

def run_single(cele_path, species, name):
    runner = Top_reads_umi_stat_run(cele_path, species, name)
    runner.run()

def main():
    mapfile = sys.argv[1]
    cele_path_list, species_list, name_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, cele_path_list, species_list, name_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
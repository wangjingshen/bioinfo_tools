import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Rename():
    def __init__(self, celescope_path, new_path, old_name, new_name):
        self.celescope_path = celescope_path
        self.new_path = new_path
        self.old_name = old_name
        self.new_name = new_name
    
    def run_rename(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/release_celescope/release_celescope.py '
            f'--celescope_path {self.celescope_path} '
            f'--new_path {self.new_path} '
            f'--old_name {self.old_name} '
            f'--new_name {self.new_name} '
        )
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_rename()

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    celescope_path_list = df_mapfile['celescope_path']
    new_path_list = df_mapfile['new_path']
    old_name_list = df_mapfile['old_name']
    new_name_list = df_mapfile['new_name']

    return celescope_path_list, new_path_list, old_name_list, new_name_list

def run_single(celescope_path, new_path, old_name, new_name):
    runner = Rename(celescope_path, new_path, old_name, new_name)
    runner.run()

def main():
    mapfile = sys.argv[1]
    celescope_path_list, new_path_list, old_name_list, new_name_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, celescope_path_list, new_path_list, old_name_list, new_name_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
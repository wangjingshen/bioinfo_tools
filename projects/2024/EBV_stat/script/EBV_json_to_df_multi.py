import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Json_to_df_single():
    def __init__(self, cele_path, prefix, mod, outdir):
        self.cele_path = cele_path
        self.prefix = prefix
        self.mod = mod
        self.outdir = outdir
    
    def run_json_to_df(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/EBV_json_to_df.py '
            f'--cele_path {self.cele_path} '
            f'--prefix {self.prefix} '
            f'--mod {self.mod} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_json_to_df()

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    cele_path_list = df_mapfile['cele_path']
    prefix_list = df_mapfile['prefix']
    mod_list = df_mapfile['mod']
    outdir_list = df_mapfile['outdir']

    return cele_path_list, prefix_list, mod_list, outdir_list

def run_single(cele_path, prefix, mod, outdir):
    runner = Json_to_df_single(cele_path, prefix, mod, outdir)
    runner.run()

def main():
    mapfile = sys.argv[1]
    cele_path_list, prefix_list, mod_list, outdir_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, cele_path_list, prefix_list, mod_list, outdir_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
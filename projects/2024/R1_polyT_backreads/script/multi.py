import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class R1_polyT_backreads_run():
    def __init__(self, probe_res, name, outdir):
        self.probe_res = probe_res
        self.name = name
        self.outdir = outdir
    
    def run(self):
        cmd = (
            f'python /SGRNJ06/randd/USER/wangjingshen/script/R1_polyT_backreads/script/R1_polyT_backreads.py '
            f'--probe_res {self.probe_res} '
            f'--name {self.name} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    probe_res_list = df_mapfile['probe_res']
    name_list = df_mapfile['name']
    outdir_list = df_mapfile['outdir']
    return probe_res_list, name_list, outdir_list

def run_single(probe_res, name, outdir):
    runner = R1_polyT_backreads_run(probe_res, name, outdir)
    runner.run()

def main():
    mapfile = sys.argv[1]
    probe_res_list, name_list, outdir_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, probe_res_list, name_list, outdir_list):
            print('Done.')


if __name__ == '__main__':
    main()
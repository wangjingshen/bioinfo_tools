import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Signac_run():
    def __init__(self, analysis_dir, mod, sname, resolution):
        self.analysis_dir = analysis_dir
        self.mod = mod
        #self.peak_h5 = peak_h5
        #self.metadata = metadata
        #self.fragments = fragments
        self.sname = sname
        self.resolution = resolution
    
    def run(self):
        cmd = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script/Signac/script/analysis.R '
            f'--analysis_dir {self.analysis_dir} '
            f'--mod {self.mod} '
            f'--sname {self.sname} '
            f'--resolution {self.resolution} '
            f'--peak_h5 None '
            f'--sname None '
            f'--resolution None '            
        )
        subprocess.check_call(cmd, shell=True)


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    analysis_dir_list = df_mapfile['analysis_dir']
    mod_list = df_mapfile['mod']
    sname_list = df_mapfile['sname']
    resolution_list = df_mapfile['resolution']

    return analysis_dir_list, mod_list,  sname_list, resolution_list

def run_single(analysis_dir, mod, sname, resolution):
    runner = Signac_run(analysis_dir, mod, sname, resolution)
    runner.run()

def main():
    mapfile = sys.argv[1]
    analysis_dir_list, mod_list, sname_list, resolution_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, analysis_dir_list, mod_list, sname_list, resolution_list):
            print('Done')


if __name__ == '__main__':
    main()
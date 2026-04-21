import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Tss_frag_run():
    def __init__(self, fragment, spname, gname, species, outdir):
        self.fragment = fragment
        self.spname = spname
        self.gname = gname
        self.species = species
        self.outdir = outdir
    
    def run(self):
        cmd = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script/ATAC_QC/script/Tss_frag.R '
            f'--fragment {self.fragment} '
            f'--spname {self.spname} '
            f'--gname {self.gname} '
            f'--species {self.species} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep=" ")
    fragment_list = df_mapfile['fragment']
    spname_list = df_mapfile['spname']
    gname_list = df_mapfile['gname']
    species_list = df_mapfile['species']
    outdir_list = df_mapfile['outdir']
    return fragment_list, spname_list, gname_list, species_list, outdir_list

def run_single(fragment, spname, gname, species, outdir):
    runner = Tss_frag_run(fragment, spname, gname, species, outdir)
    runner.run()

def main():
    mapfile = sys.argv[1]
    fragment_list, spname_list, gname_list, species_list, outdir_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        executor.map(run_single, fragment_list, spname_list, gname_list, species_list, outdir_list)


if __name__ == '__main__':
    main()
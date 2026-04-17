import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd


class CelliD():
    def __init__(self, rds, resolution, species, mode, organ, ref, nfeatures, outdir):
        self.rds = rds
        self.resolution = resolution
        self.species = species
        self.mode = mode
        self.organ = organ
        self.ref = ref
        self.nfeatures = nfeatures
        self.outdir = outdir

    def run_CelliD(self):
        cmd = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/project/CelliD/script/run_CelliD.R '
            f'--rds {self.rds} '
            f'--resolution {self.resolution} '
            f'--species {self.species} '
            f'--mode {self.mode} '
            f'--organ {self.organ} '
            f'--ref {self.ref} '
            f'--nfeatures {self.nfeatures} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_CelliD()


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\t')
    rds_list = df_mapfile['rds']
    resolution_list = df_mapfile['resolution']
    species_list = df_mapfile['species']
    mode_list = df_mapfile['mode']
    organ_list = df_mapfile['organ']
    ref_list = df_mapfile['ref']
    nfeatures_list = df_mapfile['nfeatures']
    outdir_list = df_mapfile['outdir']
    
    return rds_list, resolution_list, species_list, mode_list, organ_list, ref_list, nfeatures_list, outdir_list


def run_single(rds, resolution, species, mode, organ, ref, nfeatures, outdir):
    runner = CelliD(rds, resolution, species, mode, organ, ref, nfeatures, outdir)
    runner.run()


def main():
    mapfile = sys.argv[1]
    rds_list, resolution_list, species_list, mode_list, organ_list, ref_list, nfeatures_list, outdir_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=2) as executor:
        for result in executor.map(run_single, rds_list, resolution_list, species_list, mode_list, organ_list, ref_list, nfeatures_list, outdir_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
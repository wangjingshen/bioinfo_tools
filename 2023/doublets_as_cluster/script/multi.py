import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Doublets_as_cluster_run():
    def __init__(self, rds, mix_species, outdir, name):
        self.rds = rds
        self.mix_species = mix_species
        self.outdir = outdir
        self.name = name
    
    def run(self):
        cmd = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script/doublets_as_cluster/script/analysis.R '
            f'--rds {self.rds} '
            f'--mix_species {self.mix_species} '
            f'--outdir {self.outdir} '
            f'--name {self.name} '
        )
        #activate_environment('r4.1_env')
        subprocess.check_call(cmd, shell=True)


#def activate_environment(env_name):
#    activate_cmd = f'source activate {env_name}'
#    subprocess.check_call(activate_cmd, shell=True)

#def activate_conda_environment(env_name):
#    activate_cmd = f"conda activate {env_name}"
#    subprocess.call(['bash', '-c', activate_cmd])

# 切换到名为 "liutao_r" 的 Conda 环境
#activate_conda_environment("liutao_r")


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\t")
    rds_list = df_mapfile['rds']
    mix_species_list = df_mapfile['mix_species']
    outdir_list = df_mapfile['outdir']
    name_list = df_mapfile['name']
    return rds_list, mix_species_list, outdir_list, name_list

def run_single(rds, mix_species, outdir, name):
    runner = Doublets_as_cluster_run(rds, mix_species, outdir, name)
    runner.run()

def main():
    mapfile = sys.argv[1]
    rds_list, mix_species_list, outdir_list, name_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        for result in executor.map(run_single, rds_list, mix_species_list, outdir_list, name_list):
            print(result, 'done')


if __name__ == '__main__':
    main()
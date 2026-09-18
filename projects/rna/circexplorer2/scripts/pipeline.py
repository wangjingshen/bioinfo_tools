import subprocess
import pandas as pd
import argparse
import os
import sys
from pathlib import Path
import time


sys.path.append("/SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/")
from utils.utils import find_file, mkdir, logger, execute_cmd, make_space_input
pipeline_root = Path(__file__).resolve().parents[1]
from merge_circ_annotation import merge_circ_annotation
from build_circ_matrix import build_circ_matrix

class Circexplorer2():
    def __init__(self, celescope_dir:str, name:str, refFlat:str, reference_fa:str, match_window:int):
        self.name = name
        self.celescope_chimeric = f'{celescope_dir}/01.starsolo/{self.name}_Chimeric.out.junction'
        self.barcode = f'{celescope_dir}/outs/filtered/barcodes.tsv.gz'
        self.refFlat = refFlat
        self.reference_fa = reference_fa
        self.match_window = match_window


    def get_bc_umi_chimeric(self) -> None:
        mkdir(f'circexplorer2/{self.name}/01.chimeric/')
        cmd = [
            f'tail -n +2 {self.celescope_chimeric} | awk "NF==22" > circexplorer2/{self.name}/01.chimeric/{self.name}_valid_Chimeric.out.junction'
        ]
        execute_cmd(cmd)
    
    def circexplorer_parse(self) -> None:
        mkdir(f'circexplorer2/{self.name}/02.parse/')
        cmd = (f'/SGRNJ06/randd/USER/wangjingshen/soft/miniforge3/envs/circexplorer2_env/bin/CIRCexplorer2 parse '
               f'-t STAR '
               f'-b circexplorer2/{self.name}/02.parse/{self.name}_back_spliced_junction.bed '
               f'circexplorer2/{self.name}/01.chimeric/{self.name}_valid_Chimeric.out.junction ')
        execute_cmd(cmd)

    def circexplorer_anno(self) -> None:
        mkdir(f'circexplorer2/{self.name}/03.annotate/')
        cmd = (f'/SGRNJ06/randd/USER/wangjingshen/soft/miniforge3/envs/circexplorer2_env/bin/CIRCexplorer2 annotate '
               f'-r {self.refFlat} '
               f'-g {self.reference_fa} '
               f'-b circexplorer2/{self.name}/02.parse/{self.name}_back_spliced_junction.bed '
               f'-o circexplorer2/{self.name}/03.annotate/{self.name}_circularRNA_known.txt '
               f'--no-fix ')
        execute_cmd(cmd)

    def build_mtx(self) -> None:
        mkdir(f'circexplorer2/{self.name}/04.matrix/')
        merge_circ_annotation(self.name, self.match_window)
        build_circ_matrix(self.name, self.barcode)

    def run(self) -> None:
        self.get_bc_umi_chimeric()
        self.circexplorer_parse()
        self.circexplorer_anno()
        self.build_mtx()

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--celescope_dir', help='celescope_dir', required=True)
    parsers.add_argument('--name', help='name', required=True)
    parsers.add_argument('--refFlat', help='refFlat', required=True)
    parsers.add_argument('--reference_fa', help='reference_fa', required=True)
    parsers.add_argument('--match_window', default = 10, type =int, help='match_window', required=True)

    args = parsers.parse_args()
    runner = Circexplorer2(args.celescope_dir, args.name, args.refFlat, args.reference_fa, args.match_window) 
    runner.run()

if __name__ == '__main__':
    main()

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


class SpaceVirus():
    def __init__(self, space_dir:str, virus_df:str, name:str):
        self.space_dir = space_dir
        self.virus_df = virus_df
        self.name = Path(name)

    def run_space_virus(self) -> None:
        self.space_input = make_space_input(self.space_dir, self.name)
        cmd = (f'Rscript {pipeline_root}/scripts/space_virus.R '
              f'--space_input {self.space_input} '
              f'--virus_df {self.virus_df} '
              f'--name {self.name} ')
        execute_cmd(cmd)
    
    def rm_file(self) -> None:
        cmd = f'rm -rf {self.space_input}'
        execute_cmd(cmd)

    def run(self) -> None:
        self.run_space_virus()
        self.rm_file()

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--space_dir', help='space_dir', required=True)
    parsers.add_argument('--virus_df', help='virus_df', required=True)
    parsers.add_argument('--name', help='name', required=True)

    args = parsers.parse_args()
    runner = SpaceVirus(args.space_dir, args.virus_df, args.name) 
    runner.run()

if __name__ == '__main__':
    main()

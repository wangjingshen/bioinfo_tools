import pandas as pd
import argparse
import os
import subprocess


class Doublets_split_mtx:
    def __init__(self, args):
        self.mtx = args.mtx
        self.doublets_species = args.doublets_species
        self.prefix = args.prefix
        self.outdir = args.outdir
    
    def run_analysis(self):
        cmd1 = (
            f"Rscript /SGRNJ06/randd/USER/wangjingshen/script/doublets_split_mtx/script/analysis.R "
            f"--mtx {self.mtx} "
            f"--doublets_species {self.doublets_species} "
            f"--outdir {self.outdir} "
            f"--prefix {self.prefix} "
        )
        subprocess.check_call(cmd1, shell = True)

    def tar(self):
        cmd1 = f"tar -cf {self.outdir}/{self.prefix}_human_cellranger3_matrix_10X.tar --directory {self.outdir}/{self.prefix}_human ."
        cmd2 = f"tar -cf {self.outdir}/{self.prefix}_mouse_cellranger3_matrix_10X.tar --directory {self.outdir}/{self.prefix}_mouse ."
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
    
    def run(self):
        self.run_analysis()
        self.tar()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--mtx', help='fj path', required=True)
    parsers.add_argument('--doublets_species', help='doublets species', required=True)
    parsers.add_argument('--outdir', help='dir to save result file', required=True)
    parsers.add_argument('--prefix', help='prefix', required=True)


    args = parsers.parse_args()
    runner = Doublets_split_mtx(args) 
    runner.run()

if __name__ == '__main__':
    main()
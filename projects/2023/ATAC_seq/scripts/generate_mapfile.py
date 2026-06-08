import os
import pysam
import pandas as pd
import subprocess
import argparse
from collections import defaultdict
import glob

class Generate_mapfile():
    def __init__(self, args):
        self.raw_mapfile = args.raw_mapfile
        self.species = args.species
        self.R2_10X = args.R2_10X

    def make_mapfile(self):
        raw_mapfile = pd.read_csv(self.raw_mapfile,sep="\t",header=None)
        new_mapfile = raw_mapfile.copy()
        new_mapfile.columns = ["prefix", "fastq_path", "run_path"]
        new_mapfile['species'] = self.species
        new_mapfile['outdir'] = new_mapfile["run_path"]+"_outdir"
        new_mapfile['R2_10X'] = self.R2_10X
        new_mapfile.to_csv("mapfile",sep="\t",index=False)
    
    def run(self):
        self.make_mapfile()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--raw_mapfile', help='raw mapfile', required=True)
    parsers.add_argument('--species', help='species', required=True)
    parsers.add_argument('--R2_10X', help='R2_10X, pbmc,tiny,10X,sc', required=True)

    args = parsers.parse_args()
    runner = Generate_mapfile(args) 
    runner.run()

if __name__ == '__main__':
    main()

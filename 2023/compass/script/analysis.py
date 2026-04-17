import pysam
import pandas as pd
from collections import defaultdict
import pyranges as pr
from functools import reduce
import glob
import os
import argparse
import subprocess
import faulthandler
faulthandler.enable()  # debug
from Bio.Align.Applications import MuscleCommandline
from io import StringIO
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from xopen import xopen
from itertools import groupby
import sys
sys.path.append('/SGRNJ06/randd/USER/wangjingshen/script/tools')

import utils as utils



class Compass:
    def __init__(self, rds, h5ad, mode, threads, outdir, name):
        self.rds = rds
        self.h5ad = h5ad
        self.mode = mode
        self.subset = subset
        self.subset_var = subset_var
        self.threads = int(threads)
        self.outdir = outdir
        self.name = name
    
    def compass_input(self):
        if(self.mode == "rds"):
            cmd = f'Rscript /SGRNJ06/randd/USER/wangjingshen/script/compass/script/get_input.R '
                        f'--rds {self.rds} '
                        f'--outdir {self.outdir} '
                        f'--name {self.name} '
            subprocess.check_call(cmd, shell = True)
        
        #if(self.mode == 'h5ad'):


    def compass_run(self):
        cmd = f'compass --data {self.outdir}/compass_input{self.name}.tsv '
                f'--num-processes {self.threads}'
                f'--species {self.species}'
        subprocess.check_call(cmd, shell = True)
        

    def compass_plot(self):



def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--rds', help = 'rds')
    parsers.add_argument('--h5ad', help = 'h5ad')
    parsers.add_argument('--mode', help = 'mode, rds or h5ad', default = 'rds')
    parsers.add_argument('--species', help = 'species, homo_sapiens or mouse', default = 'homo_sapiens')
    parsers.add_argument('--outdir', help = 'outdir', default = './')
    parsers.add_argument('--name', help = 'name', default = '')
    args = parsers.parse_args()

    if not os.path.exists(f"{args.outdir}/"):
        os.system(f"mkdir -p {args.outdir}/")
    
    runner = Compass(args.rds, args.h5ad, args.mode, args.species, args.outdir, args.name) 
    runner.run()

if __name__ == '__main__':
    main()
import os
import subprocess
import argparse
import pandas as pd

def run(raw_mapfile, outdir):
    raw_mapfile = pd.read_csv(raw_mapfile, sep="\t")
    raw_mapfile['outdir'] = outdir
    raw_mapfile.to_csv('mapfile', index = False, sep="\t")

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--raw_mapfile', help='mapfile', required=True)
    parsers.add_argument('--outdir', help='outdir')

    args = parsers.parse_args()
    run(args.raw_mapfile, args.outdir)

if __name__ == '__main__':
    main()        
import os
import subprocess
import argparse
import pandas as pd


def merge_mapfile(mapfile, name):
    df_mapfile = pd.read_csv(mapfile, header= None, sep="\t")
    mapfile_list = df_mapfile.iloc[:,0].to_list()
    df = pd.DataFrame()
    for filename in mapfile_list:
        df = df.append(pd.read_csv(filename, header=None, sep="\t"), ignore_index=True)
    
    df.rename(columns={0:'library_id',1:'path',2:'sample_name'}, inplace=True)
    df.sort_values(by=['sample_name', 'path'], inplace=True)
    df.to_csv(f'{name}.mapfile', sep="\t", index=False, header=False)

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--mapfile', help='mapfile', required=True)
    parsers.add_argument('--name', help='name', required=True)

    args = parsers.parse_args()
    merge_mapfile(args.mapfile, args.name) 

if __name__ == '__main__':
    main()
        
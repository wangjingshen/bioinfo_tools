import argparse
import os
import subprocess
import pandas as pd
import numpy as np


def get_reaction_consistencies(compass_reaction_penalties, min_range=1e-3):
    """
    Converts the raw penalties outputs of compass into scores per reactions where higher numbers indicate more activity
    """
    df = -np.log(compass_reaction_penalties + 1)
    df = df[df.max(axis=1) - df.min(axis=1) >= min_range]   # filter reactions
    df = df - df.min().min()
    return df

def run(reactions, metadata, outdir, name):
    reaction_penalties = pd.read_csv(reactions, sep="\t", index_col = 0)
    reaction_consistencies = get_reaction_consistencies(reaction_penalties)
    # metadata
    reaction_metadata = pd.read_csv(metadata, index_col = 0)

    test =reaction_consistencies.reset_index().reset_index().iloc[:,0:2]
    test['reaction'] = test['index'].apply(lambda x: x[:-4])
    test1 = test.merge(reaction_metadata, how='left', left_on='reaction', right_index=True, validate='m:1')
    # check method 1
    #reaction_consistencies.index.to_list() == test1['index'].to_list()
    # check method 2
    #sum([(reaction_consistencies.index[i]) == (test1['index'][i]) for i in range(6490)])
    reaction_consistencies.index = reaction_consistencies.index + "___" + test1.subsystem + "___" + test1.reaction_name    
    reaction_consistencies.to_csv(f'{outdir}/reaction_consistencies{name}.tsv')


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--reactions', help='reactions', required=True)
    parsers.add_argument('--metadata', help='metadata')
    parsers.add_argument('--outdir', help='outdir')
    parsers.add_argument('--name', help='name')

    args = parsers.parse_args()
    if args.metadata == None:
        args.metadata="/SGRNJ06/randd/USER/wangjingshen/script/compass/data/reaction_metadata.csv"
    if args.outdir == None:
        args.outdir="./"
    if not os.path.exists(args.outdir):
        os.system(f'mkdir -p {args.outdir}')
        #cmd1 = f'mkdir -p {outdir}'
        #subprocess.check_call(cmd1, shell = True)
    if args.name == None:
        args.name=""
    if args.name!='':
        args.name = "_" + args.name
    run(args.reactions, args.metadata, args.outdir, args.name)

if __name__ == '__main__':
    main()

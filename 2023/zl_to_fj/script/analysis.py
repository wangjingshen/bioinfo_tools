import os
import subprocess
import argparse
import pandas as pd


def Zl_to_fj(zl_mapfile, fj_mapfile):

    zl_mapfile = pd.read_csv(zl_mapfile, sep="\t", header=None)
    zl_mapfile.rename(columns={0:'library_id', 1:'path', 2:'sample_name'}, inplace=True)
    zl_mapfile.sort_values(by=['sample_name', 'path'])

    fj_mapfile = pd.read_csv(fj_mapfile, sep="\t", header=None)
    fj_mapfile.rename(columns={0:'library_id',1:'path',2:'sample_name',3:'zl_path'}, inplace=True)
    fj_mapfile.sort_values(by=['sample_name','path'], inplace=True)
    fj_mapfile.drop_duplicates(['zl_path'], inplace=True)
    fj_mapfile['zl_sample_name'] = fj_mapfile.zl_path.map( lambda x : x.split('/')[-1])

    mapfile_end = zl_mapfile.merge(fj_mapfile.iloc[:,3:5], left_on='sample_name', right_on='zl_sample_name', how='left').iloc[:,0:4]  # del zl_sample_name
    #zl_mapfile.dropna()  # del no fj
    mapfile_end.sort_values(by=['zl_path'], inplace=True)
    mapfile_end.to_csv('zl_to_fj.mapfile', sep="\t", index=False, header=False)

    cmd1 = f"multi_capture_virus \
                --mapfile zl_to_fj.mapfile \
                --virus_genome /SGRNJ03/randd/test_rd/dxh/data/HBV_genome/ \
                --thread 4 \
                --not_consensus \
                --outFilterMatchNmin 80 \
                --umi_threshold otsu \
                --allowNoPolyT \
                --gtf /SGRNJ06/randd/USER/wangjingshen/project/huashan_HBV/data/2022-10-18_gtf/HBV.gtf"
    cmd2 = f'sjm log/sjm.job'
    subprocess.check_call(cmd1, shell = True)
    subprocess.check_call(cmd2, shell = True)



def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--zl_mapfile', help='zl_mapfile, split by ,', required=True)
    parsers.add_argument('--fj_mapfile', help='fj_mapfile, split by ,', required=True)

    args = parsers.parse_args()
    Zl_to_fj(args.zl_mapfile, args.fj_mapfile) 

if __name__ == '__main__':
    main()
        
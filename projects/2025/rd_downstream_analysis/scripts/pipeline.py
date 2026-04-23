import pandas as pd
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor
import json

def qc_plot(csv, spname = None):
    '''
    01 qc plot
    '''
    cmd1 = (
        f'/SGRNJ06/randd/USER/wangjingshen/soft/miniforge3/envs/r4.1_env_bak/bin/Rscript '
        f'/SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/2025/rd_downstream_analysis/scripts/qc_plot.R '
        f'--csv {csv} '
        f'--spname {spname} '
        f'--outdir 01.qc_plot/ '
    )
    subprocess.check_call(cmd1, shell=True)

def getrds(matrix_10X, spname, gname, species):
    '''
    02 getrds
    '''
    cmd1 = (
        f'/SGRNJ01/Public/Software/conda_env/r4.1_env/bin/Rscript '
        f'/SGRNJ06/randd/USER/wangjingshen/bioinfo_tools/projects/2022/seurat/scripts/seurat.R '
        f'--matrix_10X {matrix_10X} '
        f'--spname {spname} '
        f'--gname {gname} '
        f'--species {species} '
        f'--outdir 02.rds/ '
    )
    print(cmd1)
    subprocess.check_call(cmd1, shell=True)


def getrds_cjj(matrix_10X, spname, gname, species):
    '''
    02 getrds
    '''
    cmd1 = (
        f'/SGRNJ/Public/Software/conda_env/CJJV3/bin/Rscript '
        f' /SGRNJ06/randd/USER/cjj/Script/Self-Used-Scripts/seurat/Seurat3_V2.R '
        f'--expM {matrix_10X} '
        f'--spname {spname} '
        f'--ugname {gname} '
        f'--prefix rds_ '
        f'--outdir 02.rds/ '
        f'--step 1,2 '
        f'--is_10X 10X'
    )
    subprocess.check_call(cmd1, shell=True)


def probe(rna_mapfile, genomeDir, probe_file, probe_file_mode):
    '''
    03 probe
    just run sample,barcode
    '''
    cmd1 = (
        f'multi_rna '
        f'--mapfile {rna_mapfile} '
        f'--genomeDir {genomeDir} '
        f'--mod sjm '
        f'--probe_file {probe_file} '
        f'--probe_file_mode {probe_file_mode} '
        f'--steps_run sample,barcode'

    )
    subprocess.check_call(cmd1, shell=True)

def probe_run():
    '''
    03 probe run sjm
    '''
    cmd1 = f'sjm sjm/sjm.job'
    subprocess.check_call(cmd1, shell=True)


def parse_mapfile(mapfile):
    '''
    mapfile has 4 columns, which are dir, spname, gname, species in order.
    '''
    df_mapfile = pd.read_csv(mapfile, sep='\s+', header = None)
    dir_list = df_mapfile.iloc[:,0]
    spname_list = df_mapfile.iloc[:,1]
    gname_list = df_mapfile.iloc[:,2]

    return dir_list, spname_list, gname_list

def main():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help="mapfile", required=True)
    parser.add_argument('--rna_mapfile', help="rna_mapfile", required=True)
    parser.add_argument('--species', help="species", required=True)
    #parser.add_argument('--version', help="version", required=True)
    parser.add_argument('--genomeDir', help="genomeDir", default='/SGRNJ06/randd/public/genome/rna/hs/hs_ensembl_99')  # probe run barcode, not use, set default value
    parser.add_argument('--probe_file', help="probe_file", required=True)
    parser.add_argument('--probe_file_mode', help="probe_file_mode", required=True)

    args = parser.parse_args()
    dir_list, spname_list, gname_list = parse_mapfile(args.mapfile)
    matrix_10X = ",".join(dir_list + "/" + spname_list + "/outs/filtered/")  
    spname = ",".join(spname_list)
    gname = ",".join(gname_list)
    species = args.species
    #version = args.version
    rna_mapfile = args.rna_mapfile.split(",")

    if(dir_list.nunique() > 1):
        csv = ",".join(pd.Series(dir_list.unique()).str.split(",").apply(lambda x: ','.join(x) + "/merge.csv"))
    else:
        csv = dir_list[0] + "/merge.csv"

    if(len(rna_mapfile)>1):
        rna_mapfile = pd.concat([pd.read_csv(i,sep="\t",header=None) for i in rna_mapfile])
        rna_mapfile = rna_mapfile[ rna_mapfile.iloc[:,2].isin(spname.split(","))]
        rna_mapfile.to_csv("rna_mapfile", sep="\t", index = False, header = False)
        rna_mapfile = f'rna_mapfile'

    if(len(rna_mapfile)==1):
        rna_mapfile = rna_mapfile[0]

    genomeDir = args.genomeDir
    probe_file = args.probe_file
    probe_file_mode = args.probe_file_mode

    # run
    qc_plot(csv, spname)
    getrds(matrix_10X, spname, gname, species)
    probe(rna_mapfile, genomeDir, probe_file, probe_file_mode)


if __name__ == '__main__':
    main()
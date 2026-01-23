import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import argparse
from pathlib import Path
root = Path(__file__).resolve().parents[1]

def parse_mapfile(mapfile, gene_list, database):
    df_mapfile = pd.read_csv(mapfile, sep='\s+')
    df_mapfile['gene_list'] = gene_list
    df_mapfile['database'] = database

    vcf_list = df_mapfile['vcf']
    sample_list = df_mapfile['sample']
    match_dir_list = df_mapfile['match_dir']
    gene_list_list = df_mapfile['gene_list']
    database_list = df_mapfile['database']

    return vcf_list, sample_list, match_dir_list, gene_list_list, database_list


def run_single(vcf, sample, match_dir, gene_list, database):
    cmd1 = (
        f'python {root}/script/filter_snp.py '
        f'--vcf {vcf} '
        f'--outdir {sample}/04.filter_snp/ '
        f'--sample {sample} '
        f'--bcftools_filter "QUAL>=100" '
        f'--ref_threshold_method otsu '
        f'--alt_threshold_method otsu '
        f'--ref_min_support_read 2 '
        f'--alt_min_support_read 2 '
        f'--vaf 0.2 '
        f'--subparser_assay snp '
        f'--thread 10 '
    )
    subprocess.check_call(cmd1, shell=True)

    sample_RNA = sample.replace("Target", "RNA")
    print(sample_RNA)
    cmd2 = (
        f'python {root}/script/analysis_snp.py ' 
        f'--outdir {sample}/05.analysis_snp '
        f'--sample {sample} '
        f'--thread 10 '
        f'--gene_list {gene_list} '
        f'--database {database} ' 
        f'--plot_top_n 20 '
        f'--match_dir {match_dir} '
        f'--vcf {sample}/04.filter_snp/{sample}_filtered.vcf '
        f'--subparser_assay snp '
    )

    subprocess.check_call(cmd2, shell=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapfile", help="mapfile")
    parser.add_argument("--database", help="database")
    parser.add_argument("--gene_list", help="gene_list")
    args = parser.parse_args()

    vcf_list, sample_list, match_dir_list, gene_list_list, database_list = parse_mapfile(args.mapfile, args.gene_list, args.database)
    with ProcessPoolExecutor(max_workers=1) as executor:
        executor.map(run_single, vcf_list, sample_list, match_dir_list, gene_list_list, database_list)

if __name__ == '__main__':
    main()
import subprocess
import argparse
from pathlib import Path
import sys

dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from utils.utils import find_file, mkdir, logger, execute_cmd

microbe_bwa_image = '/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_microbe/pathseq_microbe.fa.img'
microbe_dict = '/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_microbe/pathseq_microbe.dict'
microbe_taxonomy_file = '/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_microbe/pathseq_taxonomy.db'


def write_cmd(mapfile, species, sample):
    if(species == "mouse"):
        genome = "/SGRNJ06/randd/public/genome/rna/mmu/mmu_ensembl_110_nofilter"
        filter_bwa_image = "/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_mm39/Mus_musculus.GRCm39.dna.primary_assembly.fa.img"
        kmer_file = "/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_mm39/Mus_musculus.GRCm39.dna.primary_assembly.bfi"

    if(species == "human"):
        genome = "/SGRNJ06/randd/public/genome/rna/hs/hs_ensembl_110_nofilter/"
        filter_bwa_image = "/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_host/pathseq_host.fa.img"
        kmer_file = "/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_host/pathseq_host.fa.fai"
    

    cmd1 = (f'/SGRNJ/Public/Software/conda_env/celescope3.0.0/bin/multi_space '
            f' --mapfile {mapfile} '
            f' --genomeDir {genome} '
            f' --outFilterMatchNmin 15 '
            f' --mod shell '
            f' --step sample,starsolo ')
    execute_cmd(cmd1)

    cmds = [f'head -n 4 shell/{sample}.sh > shell/cmd_space.sh',
            f'sed -i "s#celescope#/SGRNJ/Public/Software/conda_env/celescope3.0.0/bin/celescope#" shell/cmd_space.sh',
            f'''sed -i "s/space starsolo/space starsolo --STAR_param '--outSAMunmapped Within'/" shell/cmd_space.sh''']
    for i in cmds:
        execute_cmd(i)

    cmd2 = (f'multi_pathseq '
            f'--mapfile {mapfile} '
            f'--genomeDir {genome} '
            f'--filter_bwa_image {filter_bwa_image} '
            f'--kmer_file {kmer_file} '
            f'--microbe_bwa_image {microbe_bwa_image} '
            f'--microbe_dict {microbe_dict} '
            f'--microbe_taxonomy_file {microbe_taxonomy_file} '
            f'--downsample_reads 10000000 '
            f'--mod shell '
            f'--step pathseq,count_pathseq,analysis_pathseq')
    execute_cmd(cmd2)

    cmds = [f'tail -n 5 shell/{sample}.sh > shell/cmd_pathseq.sh ',
            f'sed -i "s#celescope#/SGRNJ/Public/Software/conda_env/celescope2.2.0/bin/celescope#" shell/cmd_pathseq.sh ',
            f'cat shell/cmd_space.sh  shell/cmd_pathseq.sh > shell/{sample}.sh ',
            f'rm shell/cmd_space.sh  shell/cmd_pathseq.sh']
    for i in cmds:
        execute_cmd(i)


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--mapfile', help='mapfile', required=True)
    parsers.add_argument('--species', help='species', required=True)
    parsers.add_argument('--sample', help='sample', required=True)
    args = parsers.parse_args()

    write_cmd(args.mapfile, args.species, args.sample)

if __name__ == '__main__':
    main()
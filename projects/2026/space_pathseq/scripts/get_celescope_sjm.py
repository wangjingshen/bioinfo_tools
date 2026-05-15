import subprocess
import argparse
from pathlib import Path
import sys

def add_root():
    root = Path(__file__).resolve().parent
    while not (root / "utils").exists():
        parent = root.parent
        if parent == root:
            raise FileNotFoundError("utils not found！")
        root = parent
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    return root

bioinfo_root = add_root()
from utils.utils import find_file, mkdir, logger, execute_cmd

microbe_bwa_image = '/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_microbe/pathseq_microbe.fa.img'
microbe_dict = '/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_microbe/pathseq_microbe.dict'
microbe_taxonomy_file = '/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_microbe/pathseq_taxonomy.db'


def write_cmd(mapfile_space, mapfile_pathseq, species, sample):
    if(species == "mouse"):
        genome = "/SGRNJ06/randd/public/genome/rna/mmu/mmu_ensembl_110_nofilter"
        filter_bwa_image = "/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_mm39/Mus_musculus.GRCm39.dna.primary_assembly.fa.img"
        kmer_file = "/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_mm39/Mus_musculus.GRCm39.dna.primary_assembly.bfi"

    if(species == "human"):
        genome = "/SGRNJ06/randd/public/genome/rna/hs/hs_ensembl_110_nofilter/"
        filter_bwa_image = "/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_host/pathseq_host.fa.img"
        kmer_file = "/SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_host/pathseq_host.fa.fai"
    

    cmd1 = (f'/SGRNJ/Public/Software/conda_env/celescope3.0.0/bin/multi_space '
            f' --mapfile {mapfile_space} '
            f' --genomeDir {genome} '
            f' --outFilterMatchNmin 15 '
            f' --mod sjm '
            f' --step sample,starsolo ')
    execute_cmd(cmd1)

    cmds = [f'cp sjm/sjm.job sjm/sjm_space.job',
            f'''sed -i "s/space starsolo/space starsolo --STAR_param '--outSAMunmapped Within'/" sjm/sjm_space.job''']
    for i in cmds:
        execute_cmd(i)

    cmd2 = (f'/SGRNJ/Public/Software/conda_env/celescope2.2.0/bin/multi_pathseq '
            f'--mapfile {mapfile_pathseq} '
            f'--genomeDir {genome} '
            f'--filter_bwa_image {filter_bwa_image} '
            f'--kmer_file {kmer_file} '
            f'--microbe_bwa_image {microbe_bwa_image} '
            f'--microbe_dict {microbe_dict} '
            f'--microbe_taxonomy_file {microbe_taxonomy_file} '
            f'--downsample_reads 10000000 '
            f'--mod sjm '
            f'--step pathseq,count_pathseq')
    execute_cmd(cmd2)

    cmds = [f'mv sjm/sjm.job sjm/sjm_pathseq.job ',
            f'sed -i "s#celescope3.0.0#celescope2.2.0#" sjm/sjm_pathseq.job',
            f'sed -i "s#01.starsolo#/outs#" sjm/sjm_pathseq.job']
    for i in cmds:
        execute_cmd(i)


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--mapfile_space', help='space mapfile', required=True)
    parsers.add_argument('--mapfile_pathseq', help='pathseq mapfile', required=True)
    parsers.add_argument('--species', help='species', required=True)
    parsers.add_argument('--sample', help='sample', required=True)
    args = parsers.parse_args()

    write_cmd(args.mapfile_space, args.mapfile_pathseq, args.species, args.sample)

if __name__ == '__main__':
    main()
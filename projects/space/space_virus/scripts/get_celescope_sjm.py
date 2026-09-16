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

miRNA_pattern = 'L22U10C8L30C8L30'

def write_cmd(space_mapfile, virus_mapfile, product, virus_genome, mirna_genome, outFilterMatchNmin, skip_intron):
    base_cmd = (
        f'/SGRNJ/Public/Software/conda_env/celescope3.0.0/bin/multi_space '
        f'--mapfile {space_mapfile} '
        f'--genomeDir {virus_genome} '
        f'--outFilterMatchNmin {outFilterMatchNmin} '
        f'--mod sjm '
        f'--step sample,starsolo '
    )

    cmd1 = base_cmd
    star_param = ""

    if product == "rna":
        if skip_intron:
            star_param = '--STAR_param " --alignIntronMax 1 --alignMatesGapMax 1 --outFilterIntronMotifs RemoveNoncanonical" '

    elif product == "miRNA":
        cmd1 += (
            f'--chemistry customized '
            f'--pattern {miRNA_pattern} '
            f'--whitelist "/SGRNJ06/randd/PROJECT/R25030501_Spatial_FFPE_tgx/20260601_FAN/DBiT96barcode.txt /SGRNJ06/randd/PROJECT/R25030501_Spatial_FFPE_tgx/20260601_FAN/DBiT96barcode.txt" '
            f'--mirna_genomeDir {mirna_genome} '
        )
        star_mirna = "--clip5pNbases 30"
        if skip_intron:
            star_mirna += " --alignIntronMax 1 --alignMatesGapMax 1 --outFilterIntronMotifs RemoveNoncanonical"
        star_param = f'--STAR_param " {star_mirna}" '

    cmd1 += star_param
    execute_cmd(cmd1)


    cmd2 = [f'mv sjm/sjm.job sjm/sjm_space.job']
    execute_cmd(cmd2)

    cmd3 = (f'/SGRNJ/Public/Software/conda_env/celescope3.0.0/bin/multi_capture_virus '
            f'--mapfile {virus_mapfile} '
            f'--virus_genomeDir {virus_genome} '
            f'--thread 4 '
            f'--outFilterMatchNmin {outFilterMatchNmin} '
            f'--not_consensus '
            f'--umi_threshold_method otsu '
            f'--mod sjm '
            f'--step count_virus,filter_virus')
    execute_cmd(cmd3)

    cmd4 = [f'mv sjm/sjm.job sjm/sjm_virus.job ']
            #f'sed -i "s#04.star_virus#/outs#" sjm/sjm_virus.job']
    execute_cmd(cmd4)


def main():
    parsers = argparse.ArgumentParser(description="CeleScope spatial virus & miRNA analysis pipeline")

    parsers.add_argument('--space_mapfile', help='Spatial sample mapfile', required=True)
    parsers.add_argument('--virus_mapfile', help='Virus reference mapfile', required=True)
    parsers.add_argument('--product', help='Product type: rna or miRNA', required=True, choices=["rna", "miRNA"])
    parsers.add_argument('--virus_genome', help='Virus STAR genome index directory', required=True)
    parsers.add_argument('--mirna_genome', help='miRNA STAR genome index directory (required when product=miRNA)')
    parsers.add_argument('--outFilterMatchNmin', default=50, type=int, help='STAR outFilterMatchNmin threshold, default: 50')
    parsers.add_argument('--skip_intron', action="store_true", help='Enable intron skipping STAR parameters (--alignIntronMax 1 etc.)')

    args = parsers.parse_args()

    if args.product == "miRNA" and args.mirna_genome is None:
        parsers.error("Error: --mirna_genome must be provided when product=miRNA")

    write_cmd(args.space_mapfile, args.virus_mapfile, args.product, args.virus_genome, args.mirna_genome, args.outFilterMatchNmin, args.skip_intron)

if __name__ == '__main__':
    main()
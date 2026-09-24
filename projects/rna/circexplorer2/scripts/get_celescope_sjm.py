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


def write_cmd(mapfile, product, genome, chimMultimapNmax):
    '''
    Parameter reference: Genome-wide analysis of drosophila circular RNAs reveals their structural and sequence properties and age-dependent neural accumulation. Cell Rep. 2014 Dec 11;9(5):1966-1980.
    chimOutType & chimJunctionOverhangMin use default value
    '''
    cmd = (
        f'/SGRNJ/Public/Software/conda_env/celescope3.0.0/bin/multi_{product} ' \
        f'--mapfile {mapfile} '
        f'--genomeDir {genome} '
        f'--mod sjm '
        f'--STAR_param "--chimSegmentMin 20 \
--chimScoreMin 1 \
--alignIntronMax 1000000 \
--outFilterMismatchNmax 4 \
--alignTranscriptsPerReadNmax 10000 \
--outFilterMultimapNmax 2 \
--chimMultimapNmax {chimMultimapNmax} \
--chimOutType Junctions \
--chimJunctionOverhangMin 20" ')

    execute_cmd(cmd)


def main():
    parsers = argparse.ArgumentParser(description="CeleScope ")

    parsers.add_argument('--mapfile', help='mapfile', required=True)
    parsers.add_argument('--product', help='Product type: rna or ffpe', required=True, choices=["rna", "ffpe"])
    parsers.add_argument('--genome', help='STAR genome index directory', required=True)
    parsers.add_argument('--chimMultimapNmax', default = 0, help='chimMultimapNmax')
    args = parsers.parse_args()

    write_cmd(args.mapfile, args.product, args.genome, args.chimMultimapNmax)

if __name__ == '__main__':
    main()
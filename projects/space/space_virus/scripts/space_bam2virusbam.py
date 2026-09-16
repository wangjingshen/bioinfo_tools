import pysam
import sys
import subprocess

from bam2virusbam import reformat_bam


def space_reformat_bam(name):
    cmd1 = f'mv {name}/outs/  {name}/04.star_virus/'
    subprocess.check_call(cmd1, shell=True)

    reformat_bam(f'{name}/04.star_virus/{name}_Aligned.sortedByCoord.out.bam', f'{name}/04.star_virus/{name}_virus_Aligned.out.bam')

if __name__ == "__main__":
    name = sys.argv[1]
    space_reformat_bam(name)
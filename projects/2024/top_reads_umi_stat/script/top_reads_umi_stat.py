import os
import argparse
import subprocess
from featureCounts import FeatureCounts

class Top_reads_umi_stat:
    def __init__(self, args):
        self.cele_path = args.cele_path
        self.bam = f"{args.cele_path}/{args.name}/outs/{args.name}_Aligned.sortedByCoord.out.bam"
        self.species = args.species
        if(args.species == "human"):
            self.gtf = "/SGRNJ06/randd/public/genome/rna/celescope_v2/hs/Homo_sapiens.GRCh38.99.gtf"
        if(args.species == "mouse"):
            self.gtf = "/SGRNJ06/randd/public/genome/rna/celescope_v2/mmu/Mus_musculus.GRCm38.99.gtf"
        self.name = args.name
    
    def run_FeatureCounts(self):
        cmd1 = (f"/SGRNJ/Public/Software/conda_env/celescope1.5.1/bin/featureCounts -s 1 "
               f"-a {self.gtf} "
               f"-o {self.name} "
               f"-R BAM "
               f"-T 4 "
               f"-t exon "
               f"{self.bam} ")
        subprocess.check_call(cmd1, shell = True)

    def run_stat(self):
        cmd1 = (f"Rscript /SGRNJ06/randd/USER/wangjingshen/script/top_reads_umi_stat/script/stat.R "
               f"--cele_path {self.cele_path} "
               f"--species {self.species} "
               f"--name {self.name} "
               f"--outdir stat")
        subprocess.check_call(cmd1, shell = True)


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--cele_path', help='cele_path', required=True)
    parsers.add_argument('--species', help='species', required=True)
    parsers.add_argument('--name', help='name', required=True)

    args = parsers.parse_args()
    runner = Top_reads_umi_stat(args) 
    runner.run_FeatureCounts()
    runner.run_stat()

if __name__ == '__main__':
    main()
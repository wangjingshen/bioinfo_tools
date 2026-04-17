import os
import argparse
import subprocess

class FeatureCounts:
    def __init__(self, args):
        self.bam = f"{args.cele_path}/outs/{args.name}_Aligned.sortedByCoord.out.bam"
        if(args.species == "human"):
            self.gtf = "/SGRNJ06/randd/public/genome/rna/celescope_v2/hs/Homo_sapiens.GRCh38.99.filter.gtf"
            #self.gtf = "/SGRNJ06/randd/public/genome/rna/celescope_v2/hs/Homo_sapiens.GRCh38.99.gtf"
        if(args.species == "mouse"):
            self.gtf = "/SGRNJ06/randd/public/genome/rna/celescope_v2/mmu/Mus_musculus.GRCm38.99.filtered.gtf"
            #self.gtf = "/SGRNJ06/randd/public/genome/rna/celescope_v2/mmu/Mus_musculus.GRCm38.99.gtf"
        self.name = args.name

    def run(self):
        cmd1 = (f"/SGRNJ/Public/Software/conda_env/celescope1.5.1/bin/featureCounts -s 1 "
               f"-a {self.gtf} "
               f"-o featureCounts/{self.name} "
               f"-R BAM "
               f"-T 4 "
               f"-t exon "
               f"{self.bam} ")
        subprocess.check_call(cmd1, shell = True)


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--cele_path', help='cele_path', required=True)
    parsers.add_argument('--species', help='species', required=True)
    parsers.add_argument('--name', help='name', required=True)

    args = parsers.parse_args()
    runner = FeatureCounts(args) 
    runner.run()

if __name__ == '__main__':
    main()

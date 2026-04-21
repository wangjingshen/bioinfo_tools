import subprocess
import argparse
import os

class Fastqc:
    def __init__(self, args):
        self.fq_path = args.fq_path
        self.outdir = args.outdir
        if not os.path.exists(f"{self.outdir}"):
            os.system(f"mkdir -p {self.outdir}")
    
    def run_fastqc(self):
        cmd1 = f"/Public/Software/miniconda2/envs/old/bin/fastqc {self.fq_path} --outdir {self.outdir}"
        subprocess.check_call(cmd1, shell = True)

    def run(self):
        self.run_fastqc()

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--fq_path', help='fq path', required=True)
    parsers.add_argument('--outdir', help='outdir', default=1)

    args = parsers.parse_args()
    runner = Fastqc(args) 
    runner.run()


if __name__ == '__main__':
    main()
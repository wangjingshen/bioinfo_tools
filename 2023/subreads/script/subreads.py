# js_env

import os
import subprocess
import argparse
import glob

class Subreads():
    def __init__(self, args):
        self.fq1 = glob.glob(os.path.join(args.fastq_path,'*R1*gz'))[0]
        self.fq2 = glob.glob(os.path.join(args.fastq_path,'*R2*gz'))[0]
        self.seed = args.seed
        self.data_size = args.data_size
        self.reads_count = int(int(args.data_size)*10**9/150/2)
        #self.percent = args.percent
        self.outdir = args.outdir
        print(self.fq1)
        print(self.fq2)
        print(self.reads_count)
        print("------")
        
        if not os.path.exists(self.outdir):
            os.system(f"mkdir -p {self.outdir}")

    def subseq(self):
        '''
        https://zhuanlan.zhihu.com/p/477002661
        '''
        cmd1 = f"seqtk sample -s {self.seed} {self.fq1} {self.reads_count} > {self.outdir}/subseq_{self.data_size}G_R1.fastq" 
        cmd2 = f"seqtk sample -s {self.seed} {self.fq2} {self.reads_count} > {self.outdir}/subseq_{self.data_size}G_R2.fastq"     
        cmd3 = f"gzip {self.outdir}/subseq_{self.data_size}G_R1.fastq"
        cmd4 = f"gzip {self.outdir}/subseq_{self.data_size}G_R2.fastq"

        #f"seqtk sample -s {self.seed} {self.fq1} {self.percent} > {self.outdir}/subseq_{self.percent}_R1.fastq"
        #f"seqtk sample -s {self.seed} {self.fq2} {self.percent} > {self.outdir}/subseq_{self.percent}_R2.fastq"
        #cmd3 = f" gzip {self.outdir}/subseq_{self.percent}_R1.fastq"
        #cmd4 = f" gzip {self.outdir}/subseq_{self.percent}_R2.fastq"

        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)

    def run(self):
        self.subseq()

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--fastq_path', help='path of fastq file', required=True)
    parsers.add_argument('--seed', help='seed', required=True)
    parsers.add_argument('--data_size', help='data size, xG', required=True)
    parsers.add_argument('--outdir', help='dir of output', required=True)

    args = parsers.parse_args()
    runner = Subreads(args) 
    runner.run()

if __name__ == '__main__':
    main()

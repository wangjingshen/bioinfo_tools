import os
import argparse
import subprocess
import sys
sys.path.append('/SGRNJ06/randd/USER/wangjingshen/script/tools')
import utils as utils

class R1_polyT_backreads:
    def __init__(self, args):
        self.probe_res = args.probe_res
        self.name = args.name
        self.outdir = args.outdir

        if not os.path.exists(self.outdir):
            os.system(f'mkdir -p {self.outdir}/map/')

    def get_no_probe_with_polyT_reads(self):
        cmd1 =  (r'''awk -F "\t" '{if($10 == "polyT" && $9 == "no_probe"){print $0}}' ''' 
                f'{self.probe_res} > {self.outdir}/{self.name}_no_probe_polyT.txt')
        cmd2 =  (r'''awk -F "\t" -v OFS="" '{print(">",$2,"\n",$6)}' '''
                 f'{self.outdir}/{self.name}_no_probe_polyT.txt > {self.outdir}/{self.name}_no_probe_polyT.fa')
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        # awk -F "\t" '{print $9}' H001_no_probe_polyT.txt | sort | uniq -c   24709463 no_probe 
        # awk -F "\t" '{print $10}' H001_no_probe_polyT.txt | sort | uniq -c  24709463 polyT

    @utils.add_log
    def cutadapt_polyT(self):
        polyT = "TTTTTTTTTT"
        cmd1 = f'cutadapt -j 8 -g {polyT} -e 0.1 -o {self.outdir}/{self.name}_cut1.fa {self.outdir}/{self.name}_no_probe_polyT.fa --minimum-length 5'
        cmd2 = f'''sed -i 's/^T\+//' {self.outdir}/{self.name}_cut1.fa'''
        cmd3 = f'cutadapt -j 8 -g {polyT} -e 0.1 -o {self.outdir}/{self.name}_cut2.fa {self.outdir}/{self.name}_cut1.fa --minimum-length 5'
        cmd4 = f'''sed -i 's/^T\+//' {self.outdir}/{self.name}_cut2.fa'''
        cmd5 = f'cutadapt -j 8 -g {polyT} -e 0.1 -o {self.outdir}/{self.name}_cut3.fa {self.outdir}/{self.name}_cut2.fa --minimum-length 5'
        cmd6 = f'''sed -i 's/^T\+//' {self.outdir}/{self.name}_cut3.fa'''
        cmd7 = f'cutadapt -j 8 -g {polyT} -e 0.1 -o {self.outdir}/{self.name}_cut4.fa {self.outdir}/{self.name}_cut3.fa --minimum-length 5'
        cmd8 = f'''sed -i 's/^T\+//' {self.outdir}/{self.name}_cut4.fa'''
        cmd9 = f'cutadapt -j 8 -g {polyT} -e 0.1 -o {self.outdir}/{self.name}_cut5.fa {self.outdir}/{self.name}_cut4.fa --minimum-length 5'
        cmd10 = f'''sed -i 's/^T\+//' {self.outdir}/{self.name}_cut5.fa'''
        cmd11 = f'cutadapt -j 8 -g {polyT} -e 0.1 -o {self.outdir}/{self.name}_cut6.fa {self.outdir}/{self.name}_cut5.fa --minimum-length 5'
        cmd12 = f'''sed -i 's/^T\+//' {self.outdir}/{self.name}_cut6.fa'''
        cmd13 = f'seqtk seq -L 5 {self.outdir}/{self.name}_cut6.fa > {self.outdir}/{self.name}_no_probe_cut_polyT.fa'
        cmd14 = f'rm {self.outdir}/{self.name}_no_probe_polyT.fa {self.outdir}/{self.name}_cut1.fa {self.outdir}/{self.name}_cut2.fa {self.outdir}/{self.name}_cut3.fa {self.outdir}/{self.name}_cut4.fa {self.outdir}/{self.name}_cut5.fa {self.outdir}/{self.name}_cut6.fa'

        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)
        subprocess.check_call(cmd5, shell = True)
        subprocess.check_call(cmd6, shell = True)
        subprocess.check_call(cmd7, shell = True)
        subprocess.check_call(cmd8, shell = True)
        subprocess.check_call(cmd9, shell = True)
        subprocess.check_call(cmd10, shell = True)
        subprocess.check_call(cmd11, shell = True)
        subprocess.check_call(cmd12, shell = True)
        subprocess.check_call(cmd13, shell = True)
        subprocess.check_call(cmd14, shell = True)

    # error
    #def fa2fq(self):
    #    cmd = f'seqtk seq -a {self.outdir}/{self.name}_no_probe_cut_polyT.fa > {self.outdir}/{self.name}_no_probe_cut_polyT.fq'
    #    subprocess.check_call(cmd, shell = True)

    def map2HBV(self):
        cmd = (f'/SGRNJ/Public/Software/conda_env/celescope1.5.1b0/bin/STAR --runThreadN 4 '
               f'--genomeDir /SGRNJ03/randd/test_rd/dxh/data/HBV_genome/ '
               f'--readFilesIn {self.outdir}/{self.name}_no_probe_cut_polyT.fa '
               f'--outFilterMultimapNmax 1 '
               f'--outFileNamePrefix {self.outdir}/map/{self.name}_no_probe_cut_polyT_ '
               f'--outSAMtype BAM Unsorted ' 
            )
        subprocess.check_call(cmd, shell = True)
    
    def bam2fa(self):
        cmd1 = f'bedtools bamtofastq -i {self.outdir}/map/{self.name}_no_probe_cut_polyT_Aligned.out.bam -fq {self.outdir}/{self.name}_no_probe_cut_polyT_HBV_Aligned.fq'
        cmd2 = f'/SGRNJ/Public/Software/conda_env/cellrank/bin/seqtk seq -A {self.outdir}/{self.name}_no_probe_cut_polyT_HBV_Aligned.fq > {self.outdir}/{self.name}_no_probe_cut_polyT_HBV_Aligned.fa'
        cmd3 = f'rm {self.outdir}/{self.name}_no_probe_cut_polyT_HBV_Aligned.fq'

        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)
        subprocess.check_call(cmd3, shell=True)

    def run(self):
        self.get_no_probe_with_polyT_reads()
        self.cutadapt_polyT()
        #self.fa2fq()
        self.map2HBV()
        self.bam2fa()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--probe_res', help='probe_res', required=True)
    parsers.add_argument('--name', help='name', required=True)
    parsers.add_argument('--outdir', help='outdir', required=True)

    args = parsers.parse_args()
    runner = R1_polyT_backreads(args) 
    runner.run()

if __name__ == '__main__':
    main()

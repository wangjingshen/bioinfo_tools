import os
import argparse
import subprocess
import sys
sys.path.append('/SGRNJ06/randd/USER/wangjingshen/script/tools')
import utils as utils
import glob

class R1_polyT_backreads_zl:
    def __init__(self, args):
        '''
        r1_fq: r1 run barcode fq
        r2_bam: r2 map to HBV
        '''
        self.r1_fq = glob.glob(f"{args.r1_cele}/*barcode*/*.fq")[0]
        self.r2_bam = glob.glob(f"{args.r2_cele}/*star_virus*/*Aligned.sortedByCoord.out.bam")[0]
        self.name = args.name

        if not os.path.exists(self.name):
            os.system(f'mkdir -p {self.name}/map/')

    def get_r2_map_HBV_R1(self):
        cmd1 = (f'samtools view {self.r2_bam} | '
                r'''awk -F "\t" '{print $1}' '''
                f'> {self.name}/{self.name}_r2_map_HBV.bc')
        cmd2 = f"seqtk subseq {self.r1_fq} {self.name}/{self.name}_r2_map_HBV.bc > {self.name}/{self.name}_r2_map_HBV_R1.fq"
        cmd3 = f'seqtk seq -A {self.name}/{self.name}_r2_map_HBV_R1.fq > {self.name}/{self.name}_r2_map_HBV_R1.fa'
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)


    @utils.add_log
    def cutadapt_polyT(self):
        polyT = "TTTTTTTTTT"
        cmd1 = f'cutadapt -j 8 -g {polyT} -e 0.1 -o {self.name}/{self.name}_cut1.fa {self.name}/{self.name}_r2_map_HBV_R1.fa --minimum-length 5'
        cmd2 = f'''sed -i 's/^T\+//' {self.name}/{self.name}_cut1.fa'''
        cmd3 = f'cutadapt -j 8 -g {polyT} -e 0.1 -o {self.name}/{self.name}_cut2.fa {self.name}/{self.name}_cut1.fa --minimum-length 5'
        cmd4 = f'''sed -i 's/^T\+//' {self.name}/{self.name}_cut2.fa'''
        cmd5 = f'cutadapt -j 8 -g {polyT} -e 0.1 -o {self.name}/{self.name}_cut3.fa {self.name}/{self.name}_cut2.fa --minimum-length 5'
        cmd6 = f'''sed -i 's/^T\+//' {self.name}/{self.name}_cut3.fa'''
        cmd7 = f'cutadapt -j 8 -g {polyT} -e 0.1 -o {self.name}/{self.name}_cut4.fa {self.name}/{self.name}_cut3.fa --minimum-length 5'
        cmd8 = f'''sed -i 's/^T\+//' {self.name}/{self.name}_cut4.fa'''
        cmd9 = f'cutadapt -j 8 -g {polyT} -e 0.1 -o {self.name}/{self.name}_cut5.fa {self.name}/{self.name}_cut4.fa --minimum-length 5'
        cmd10 = f'''sed -i 's/^T\+//' {self.name}/{self.name}_cut5.fa'''
        cmd11 = f'cutadapt -j 8 -g {polyT} -e 0.1 -o {self.name}/{self.name}_cut6.fa {self.name}/{self.name}_cut5.fa --minimum-length 5'
        cmd12 = f'''sed -i 's/^T\+//' {self.name}/{self.name}_cut6.fa'''
        cmd13 = f'seqtk seq -L 5 {self.name}/{self.name}_cut6.fa > {self.name}/{self.name}_r2_map_HBV_R1_cut_polyT.fa'  # save length >5 seq
        cmd14 = f'rm {self.name}/{self.name}_cut1.fa {self.name}/{self.name}_cut2.fa {self.name}/{self.name}_cut3.fa {self.name}/{self.name}_cut4.fa {self.name}/{self.name}_cut5.fa {self.name}/{self.name}_cut6.fa {self.name}/{self.name}_r2_map_HBV_R1.fa'

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


    def map2HBV(self):
        cmd = (f'/SGRNJ/Public/Software/conda_env/celescope1.5.1b0/bin/STAR --runThreadN 4 '
               f'--genomeDir /SGRNJ03/randd/test_rd/dxh/data/HBV_genome/ '
               f'--readFilesIn {self.name}/{self.name}_r2_map_HBV_R1_cut_polyT.fa '
               f'--outFilterMultimapNmax 1 '
               f'--outFileNamePrefix {self.name}/map/{self.name}_r2_map_HBV_R1_cut_polyT_ '
               f'--outSAMtype BAM Unsorted ' 
            )
        subprocess.check_call(cmd, shell = True)
    
    def bam2fa(self):
        cmd1 = f'bedtools bamtofastq -i {self.name}/map/{self.name}_r2_map_HBV_R1_cut_polyT_Aligned.out.bam -fq {self.name}/{self.name}_r2_map_HBV_R1_cut_polyT_HBV_Aligned.fq'
        cmd2 = f'/SGRNJ/Public/Software/conda_env/cellrank/bin/seqtk seq -A {self.name}/{self.name}_r2_map_HBV_R1_cut_polyT_HBV_Aligned.fq > {self.name}/{self.name}_r2_map_HBV_R1_cut_polyT_HBV_Aligned.fa'
        cmd3 = f'rm {self.name}/{self.name}_r2_map_HBV_R1_cut_polyT_HBV_Aligned.fq'

        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)
        subprocess.check_call(cmd3, shell=True)

    def run(self):
        self.get_r2_map_HBV_R1()
        self.cutadapt_polyT()
        self.map2HBV()
        self.bam2fa()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--r1_cele', help='r1_cele', required=True)
    parsers.add_argument('--r2_cele', help='r2_cele', required=True)
    parsers.add_argument('--name', help='name', required=True)

    args = parsers.parse_args()
    runner = R1_polyT_backreads_zl(args)
    runner.run()

if __name__ == '__main__':
    main()

import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import glob

class UMI_matrix_run():
    def __init__(self, samplename, fj_dir, rna_dir, pathseq_dir):
        '''
        fj_dir: fj run celescope rna with --STAR_param "--outSAMunmapped Within"
        rna_dir: rna run celescope rna
        pathseq_dir: pathseq dir
        samplename: sample name
        '''
        self.samplename = samplename
        self.fj_dir = fj_dir
        self.fj_map_bam = glob.glob(f"{fj_dir}/outs/*Aligned.sortedByCoord.out.bam")[0]
        self.rna_dir =rna_dir
        self.rna_bc = glob.glob(f"{rna_dir}/outs/filtered/barcodes.tsv.gz")[0]
        #self.rna_bc = os.path.join(rna_dir, "outs/filtered/barcodes.tsv.gz")
        self.pathseq_dir = pathseq_dir


    def run_UMI_matrix(self):
        cmd = (
            f'/SGRNJ06/randd/USER/wangjingshen/script/PathSeq/script/UMI_matrix.py '
            f' {self.fj_map_bam} '
            f' {self.samplename} '
            f' {self.rna_bc}'
            f' {self.pathseq_dir}/{self.samplename}_pathseq.bam '
            f' {self.pathseq_dir}/{self.samplename}_pathseq_score.txt '
            f' {self.pathseq_dir}/{self.samplename}.invadeseq.readname '
            f' {self.pathseq_dir}/{self.samplename}.invadeseq.unmap_cbub.bam '
            f' {self.pathseq_dir}/{self.samplename}.invadeseq.unmap_cbub.fasta '
            f' {self.pathseq_dir}/{self.samplename}.invadeseq.list '
            f' {self.pathseq_dir}/{self.samplename}.invadeseq.raw.readnamepath '
            f' {self.pathseq_dir}/{self.samplename}.invadeseq.genus.cell '
            f' {self.pathseq_dir}/{self.samplename}.invadeseq.genus.csv '
            f' {self.pathseq_dir}/{self.samplename}.invadeseq.validate.csv '
        )
        subprocess.check_call(cmd, shell=True)

    def run_UMI_df(self):
        cmd = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script/PathSeq/script/UMI_df.R '
            f'--samplename {self.samplename} '
            f'--rna_dir {self.rna_dir}  '
            f'--pathseq_dir {self.pathseq_dir} '
        )
        subprocess.check_call(cmd, shell=True)

   
    def run(self):
        self.run_UMI_matrix()
        self.run_UMI_df()


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\t')
    samplename_list = df_mapfile['samplename']
    fj_dir_list = df_mapfile['fj_dir']
    rna_dir_list = df_mapfile['rna_dir']
    pathseq_dir_list = df_mapfile['pathseq_dir']
    
    return samplename_list, fj_dir_list, rna_dir_list, pathseq_dir_list


def run_single(samplename, fj_dir, rna_dir, pathseq_dir):
    runner = UMI_matrix_run(samplename, fj_dir, rna_dir, pathseq_dir)
    runner.run()


def main():
    mapfile = sys.argv[1]
    samplename_list, fj_dir_list, rna_dir_list, pathseq_dir_list = parse_mapfile(mapfile)
    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, samplename_list, fj_dir_list, rna_dir_list, pathseq_dir_list)
        #for result in executor.map(run_single, samplename_list, fj_dir_list, rna_dir_list, pathseq_dir_list):
        #    print('Done.')

if __name__ == '__main__':
    main()
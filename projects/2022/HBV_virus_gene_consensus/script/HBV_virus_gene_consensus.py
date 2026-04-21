import pysam
import pandas as pd
from collections import defaultdict
import pyranges as pr
from functools import reduce
import glob
import os
import argparse
import subprocess
import faulthandler
faulthandler.enable()


class Virus_gene_consensus:
    def __init__(self, args):
        self.fj_path = args.fj_path
        self.outdir = args.outdir

    
    def mkdir(self):
        if not os.path.exists(f"{self.outdir}/1.bam/"):
            os.system(f"mkdir -p {self.outdir}/1.bam/")
        if not os.path.exists(f"{self.outdir}/2.fastq/"):
            os.system(f"mkdir -p {self.outdir}/2.fastq/")
        if not os.path.exists(f"{self.outdir}/3.consensus/"):
            os.system(f"mkdir -p {self.outdir}/3.consensus/")
        

    def get_virus_barcode_umi(self):
        '''
        count_detail: 08
        virus_tsne: 05
        '''
        count_detail = pd.read_csv(glob.glob(f"{self.fj_path}/08.count_capture_virus_mtx/*count_detail.txt")[0], sep="\t")    # all umis support min_support_reads
        virus_tsne = pd.read_csv(glob.glob(f"{self.fj_path}/*analysis_capture_virus*/*virus_tsne.tsv")[0], sep="\t", index_col=0).fillna(0)
        virus_barcode = list(virus_tsne.barcode[virus_tsne.UMI!=0])
        count_detail_virus = count_detail.loc[ count_detail.Barcode.isin(virus_barcode), ].copy()    # filter negative barcodes  # copy for deal with "A value is trying to be set on a copy of a slice from a DataFrame"
        count_detail_virus['barcode_umi'] = count_detail_virus.Barcode + "_" + count_detail_virus.UMI
        # split gene
        S_list =  list(count_detail_virus.loc[count_detail_virus['geneID'] == 'forS', 'barcode_umi'])
        X_list =  list(count_detail_virus.loc[count_detail_virus['geneID'] == 'forX', 'barcode_umi'])
        pgRNA_list =  list(count_detail_virus.loc[count_detail_virus['geneID'] == 'forpg', 'barcode_umi'])
        rcDNA_list =  list(count_detail_virus.loc[count_detail_virus['geneID'] == 'forrcDNA', 'barcode_umi']) 
        cccDNA_list =  list(count_detail_virus.loc[count_detail_virus['geneID'] == 'forcccDNA', 'barcode_umi'])

        # ref record df
        input_bam = pysam.AlignmentFile(glob.glob(f"{self.fj_path}/*star_virus*/*Aligned.sortedByCoord.out.bam")[0], "rb")
        ref_record_dict = defaultdict(list)
        for i in input_bam:
            ref_record_dict["qname"].append(i.qname)
            ref_record_dict["record"].append(i)
        ref_record_df = pd.DataFrame(ref_record_dict)
        ref_record_df['barcode_umi'] = ref_record_df['qname'].str.split("_").str[0] + "_" + ref_record_df['qname'].str.split("_").str[1]
        input_bam.close()
   
        # out bam
        def probe_bam_out(gene_barcode_umi_list, out_bam):
            input_bam = pysam.AlignmentFile(glob.glob(f"{self.fj_path}/*star_virus*/*Aligned.sortedByCoord.out.bam")[0], "rb")
            out_bam = pysam.AlignmentFile(out_bam, "wb", template = input_bam)

            ref_record_df_filter = ref_record_df[(ref_record_df.barcode_umi.isin(gene_barcode_umi_list))]
            for _, row in ref_record_df_filter.iterrows():
                out_bam.write(row["record"])  

            input_bam.close()
            out_bam.close()

        probe_bam_out(S_list, f'{self.outdir}/1.bam/S.bam')
        probe_bam_out(X_list, f'{self.outdir}/1.bam/X.bam')
        probe_bam_out(pgRNA_list, f'{self.outdir}/1.bam/pgRNA.bam')
        probe_bam_out(rcDNA_list, f'{self.outdir}/1.bam/rcDNA.bam')
        probe_bam_out(cccDNA_list, f'{self.outdir}/1.bam/cccDNA.bam')


    def bam_to_fq(self):
        '''
        bam to fastq for consensus
        '''
        cmd1 = f"/SGRNJ/Public/Software/conda_env/js_env/bin/bedtools bamtofastq -i {self.outdir}/1.bam/S.bam -fq {self.outdir}/2.fastq/S.fq"
        cmd2 = f"/SGRNJ/Public/Software/conda_env/js_env/bin/bedtools bamtofastq -i {self.outdir}/1.bam/X.bam -fq {self.outdir}/2.fastq/X.fq"
        cmd3 = f"/SGRNJ/Public/Software/conda_env/js_env/bin/bedtools bamtofastq -i {self.outdir}/1.bam/pgRNA.bam -fq {self.outdir}/2.fastq/pgRNA.fq"
        cmd4 = f"/SGRNJ/Public/Software/conda_env/js_env/bin/bedtools bamtofastq -i {self.outdir}/1.bam/rcDNA.bam -fq {self.outdir}/2.fastq/rcDNA.fq"
        cmd5 = f"/SGRNJ/Public/Software/conda_env/js_env/bin/bedtools bamtofastq -i {self.outdir}/1.bam/cccDNA.bam -fq {self.outdir}/2.fastq/cccDNA.fq"
        
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)
        subprocess.check_call(cmd5, shell = True)


    def run_consensus(self):
        def consensus(sample,fq):
            cmd1 = f"/SGRNJ/Public/Software/conda_env/celescope1.5.1b0/bin/celescope capture_virus consensus \
                --outdir {self.outdir}/3.consensus/{sample}/ \
                --sample {sample} \
                --assay capture_virus \
                --thread 4 \
                --threshold 0.5 \
                --min_consensus_read 1  \
                --fq {fq}"
            subprocess.check_call(cmd1, shell = True)
        consensus("S", f"{self.outdir}/2.fastq/S.fq")
        consensus("X", f"{self.outdir}/2.fastq/X.fq")
        consensus("pgRNA", f"{self.outdir}/2.fastq/pgRNA.fq")
        consensus("rcDNA", f"{self.outdir}/2.fastq/rcDNA.fq")
        consensus("cccDNA", f"{self.outdir}/2.fastq/cccDNA.fq")
    
    def rm_and_gzip(self):
        cmd1 = f"rm {self.outdir}/3.consensus/*report.html"
        cmd2 = f"rm -rf {self.outdir}/3.consensus/*/*tmp"
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)

        cmd3 = f"gzip {self.outdir}/2.fastq/*.fq"
        cmd4 = f"gzip {self.outdir}/3.consensus/*/*.fq"
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)

    
    def run(self):
        self.mkdir()
        self.get_virus_barcode_umi()
        self.bam_to_fq()
        self.run_consensus()
        self.rm_and_gzip()
    

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--fj_path', help='fj path', required=True)
    parsers.add_argument('--outdir', help='outdir', required=True)


    args = parsers.parse_args()
    runner = Virus_gene_consensus(args) 
    runner.run()

if __name__ == '__main__':
    main()
import pysam
import pandas as pd
from collections import defaultdict
import argparse
import os
import glob
import subprocess


class Bam_virus:
    def __init__(self, args):
        self.fj_path = args.fj_path
        input_bam_path = glob.glob(f"{self.fj_path}/*star*/*Aligned.sortedByCoord.out.bam")[0]
        virus_tsne_path = glob.glob(f"{self.fj_path}/*analysis_capture_virus*/*virus_tsne.tsv")[0]
        virus_read_count_path = glob.glob(f"{self.fj_path}/*count_capture_virus*/*virus_read_count.tsv")[0]

        self.input_bam = pysam.AlignmentFile(input_bam_path, "rb")

        self.virus_tsne = pd.read_csv(virus_tsne_path, sep="\t", index_col=0)
        self.virus_tsne = self.virus_tsne.fillna(0)   # fillna 0
        self.virus_barcode = list(self.virus_tsne['barcode'][self.virus_tsne.UMI!=0])
        self.otsu_min_support_read = int(args.otsu_min_support_read)

        self.virus_read_count = pd.read_csv(virus_read_count_path, sep="\t")
        self.virus_read_count = self.virus_read_count[(self.virus_read_count.read_count >= self.otsu_min_support_read) & (self.virus_read_count.barcode.isin(self.virus_barcode)) ]
        self.virus_barcode_umi = list(self.virus_read_count.barcode + "_" + self.virus_read_count.UMI)
        self.outdir = args.outdir
        self.name = args.name

        if not os.path.exists(self.outdir):
            os.system(f"mkdir -p {self.outdir}")
    
    def virus_bam(self):
        def virus_bam_out(input_bam, out_bam):
            ref_record_dict = defaultdict(list)
            for i in input_bam:
                #if(i.qname.split("_")[0]+"_"+i.qname.split("_")[1] in self.virus_barcode_umi):
                ref_record_dict["qname"].append(i.qname)
                ref_record_dict["record"].append(i)
            ref_record_df = pd.DataFrame(ref_record_dict)
            ref_record_df['barcode_umi'] = ref_record_df['qname'].str.split("_").str[0] + "_" + ref_record_df['qname'].str.split("_").str[1]
            ref_record_df_filter = ref_record_df[ref_record_df.barcode_umi.isin(self.virus_barcode_umi)]
            for _, row in ref_record_df_filter.iterrows():
                out_bam.write(row["record"])
            input_bam.close()
            out_bam.close()

        bam_out = f'{self.outdir}/{self.name}_virus.bam'
        virus_bam_out(self.input_bam, pysam.AlignmentFile(bam_out, "wb", template = self.input_bam))
    
    def bam2tsv(self):
        cmd = (f'samtools view {self.outdir}/{self.name}_virus.bam | '
               r'''awk -F "\t" '{print $1}' | awk -F "_" '{print $1"_"$2"\t"$3}' > '''
               f'{self.outdir}/{self.name}_virus_id.tsv'
              )
        subprocess.check_call(cmd, shell = True)
    
    def run(self):
        self.virus_bam()
        self.bam2tsv()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--fj_path', help='fj path', required=True)
    parsers.add_argument('--otsu_min_support_read', help='otsu_min_support_read', required=True)
    parsers.add_argument('--outdir', help='dir to save result file', required=True)
    parsers.add_argument('--name', help='name of result file', required=True)

    args = parsers.parse_args()
    runner = Bam_virus(args) 
    runner.run()

if __name__ == '__main__':
    main()
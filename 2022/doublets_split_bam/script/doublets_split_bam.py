import pysam
import pandas as pd
from collections import defaultdict
import argparse
import glob
import os
import subprocess

class Doublets_split_bam:
    def __init__(self, args):
        self.doublets_species = args.doublets_species
        self.input_bam = glob.glob(f"{args.fj_path}/04.star_virus/*_Aligned.sortedByCoord.out.bam")[0]
        self.virus_tsne = glob.glob(f"{args.fj_path}/06.analysis_capture_virus/*_virus_tsne.tsv")[0]
        self.read_count = glob.glob(f"{args.fj_path}/05.count_capture_virus/*_read_count.tsv")[0]
        self.otsu_min_support_read = int(args.otsu_min_support_read)
        self.name = args.name
        self.outdir = args.outdir
        if not os.path.exists(self.outdir):
            os.system(f"mkdir -p {self.outdir}")
    
    def get_species_barcode(self):
        doublets_species = pd.read_csv(self.doublets_species, sep="\t")
        self.hs_barcode = list(doublets_species.cell_id[doublets_species.Species == "human"])
        self.mm_barcode = list(doublets_species.cell_id[doublets_species.Species == "mouse"])

    def get_virus_barcode_umi(self):
        virus_tsne = pd.read_csv(self.virus_tsne, sep="\t", index_col=0).fillna(0) 
        virus_barcode = list(virus_tsne['barcode'][virus_tsne.UMI!=0]) # virus barcode 
        read_count = pd.read_csv(self.read_count, sep="\t")
        read_count = read_count[(read_count.read_count >= self.otsu_min_support_read) & (read_count.barcode.isin(virus_barcode)) ]  # otsu and virus barcode
        self.virus_barcode_umi = list(read_count.barcode + "_" + read_count.UMI) # virus barcode umi
    
    def get_bam(self):
        ref_record_dict = defaultdict(list)
        input_bam = pysam.AlignmentFile(self.input_bam, "rb")
        for i in input_bam:
            ref_record_dict["qname"].append(i.qname)
            ref_record_dict["record"].append(i)
        ref_record_df = pd.DataFrame(ref_record_dict)
        ref_record_df['barcode'] = ref_record_df['qname'].str.split("_").str[0]
        ref_record_df['barcode_umi'] = ref_record_df['qname'].str.split("_").str[0] + "_" + ref_record_df['qname'].str.split("_").str[1]

        # species
        def species_bam_out(species_barcode, out_bam):
            input_bam = pysam.AlignmentFile(self.input_bam, "rb")
            out_bam = pysam.AlignmentFile(out_bam, "wb", template = input_bam)
            ref_record_df_filter = ref_record_df[ref_record_df.barcode.isin(species_barcode)]
            for _, row in ref_record_df_filter.iterrows():
                out_bam.write(row["record"])
            input_bam.close()
            out_bam.close()
        species_bam_out(self.mm_barcode, f'{self.outdir}/{self.name}_mm.bam')
        species_bam_out(self.hs_barcode, f'{self.outdir}/{self.name}_hs.bam')

        # species + virus 
        def species_virus_bam_out(species_barcode, out_bam):
            input_bam = pysam.AlignmentFile(self.input_bam, "rb")
            out_bam = pysam.AlignmentFile(out_bam, "wb", template = input_bam)
            ref_record_df_filter = ref_record_df[ref_record_df.barcode.isin(species_barcode) & ref_record_df.barcode_umi.isin(self.virus_barcode_umi)]
            for _, row in ref_record_df_filter.iterrows():
                out_bam.write(row["record"])
            input_bam.close()
            out_bam.close()
        species_virus_bam_out(self.mm_barcode, f'{self.outdir}/{self.name}_mm_virus.bam')
        species_virus_bam_out(self.hs_barcode, f'{self.outdir}/{self.name}_hs_virus.bam')


    def sort_bam(self):
        cmd1 = f"samtools sort -@8 {self.outdir}/{self.name}_hs.bam -o {self.outdir}/{self.name}_hs.sort.bam"
        cmd2 = f"samtools sort -@8 {self.outdir}/{self.name}_hs_virus.bam -o {self.outdir}/{self.name}_hs_virus.sort.bam"
        cmd3 = f"samtools sort -@8 {self.outdir}/{self.name}_mm.bam -o {self.outdir}/{self.name}_mm.sort.bam"
        cmd4 = f"samtools sort -@8 {self.outdir}/{self.name}_mm_virus.bam -o {self.outdir}/{self.name}_mm_virus.sort.bam"
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)

    def index_bam(self):
        cmd1 = f"samtools index {self.outdir}/{self.name}_hs.sort.bam"
        cmd2 = f"samtools index {self.outdir}/{self.name}_hs_virus.sort.bam"
        cmd3 = f"samtools index {self.outdir}/{self.name}_mm.sort.bam"
        cmd4 = f"samtools index {self.outdir}/{self.name}_mm_virus.sort.bam"
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)

    
    def run(self):
        self.get_species_barcode()
        self.get_virus_barcode_umi()
        self.get_bam()
        self.sort_bam()
        self.index_bam()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--doublets_species', help='doublets species', required=True)
    parsers.add_argument('--fj_path', help='fj path', required=True)
    parsers.add_argument('--otsu_min_support_read', help='otsu min support read', required=True)
    parsers.add_argument('--name', help='name', required=True)
    parsers.add_argument('--outdir', help='dir to save result file', required=True)

    args = parsers.parse_args()
    runner = Doublets_split_bam(args) 
    runner.run()

if __name__ == '__main__':
    main()
import pysam
import pandas as pd
from collections import defaultdict
import pyranges as pr
from functools import reduce
import glob
import os
import argparse
import subprocess


class Probe_split:
    def __init__(self, args):
        self.fj_path = args.fj_path
        self.probe_anno = args.probe_anno

    def mkdir(self):
        if not os.path.exists("outdir/1.bam/"):
            os.system(f"mkdir -p outdir/1.bam/")
        if not os.path.exists("outdir/2.fastq/"):
            os.system(f"mkdir -p outdir/2.fastq/")
        if not os.path.exists("outdir/3.consensus/"):
            os.system(f"mkdir -p outdir/3.consensus/")


    def split_bam(self):
        df_merged = pd.read_csv(self.probe_anno, sep="\t", header = 0)

        # set probe type
        df_merged['probe_type'] = "tmp"
        df_merged.loc[(df_merged['probe_name'] != "no_probe") & (df_merged['is_polyT'] != "polyT") , 'probe_type'] = "probe_no_polyT"  # reads_with_probe_no_poly_T
        df_merged.loc[(df_merged['probe_name'] != "no_probe") & (df_merged['is_polyT'] == "polyT") , 'probe_type'] = "probe_with_polyT"  # reads_with_probe_with_poly_T
        df_merged.loc[(df_merged['probe_name'] == "no_probe") & (df_merged['is_polyT'] == "polyT") , 'probe_type'] = "polyT" # reads_no_probe_with_poly_T
        df_merged.loc[(df_merged['probe_name'] == "no_probe") & (df_merged['is_polyT'] != "polyT") , 'probe_type'] = "others"   # reads_no_probe_no_poly_T

        # list to split R2 bam
        S_list =  list(df_merged.loc[df_merged['probe_name'] == 'for-S-probe-3', 'query_name'])
        X_list =  list(df_merged.loc[df_merged['probe_name'] == 'for-X-probe-1', 'query_name'])
        pgRNA_list =  list(df_merged.loc[df_merged['probe_name'] == 'for-pg-probe-2', 'query_name'])
        rcDNA_list =  list(df_merged.loc[df_merged['probe_name'] == 'for-rcdna-probe-1', 'query_name']) 
        cccDNA_list =  list(df_merged.loc[df_merged['probe_name'] == 'for-cccdna_probe', 'query_name'])
        polyT_list =  list(df_merged.loc[df_merged['probe_type'] == 'polyT', 'query_name']) 
        others_list =  list(df_merged.loc[df_merged['probe_type'] == 'others', 'query_name'])

        sdict = {}
        sdict.update({
            "mismatch": 2,
            "total_reads": len(df_merged),
            "reads_with_probe": sum( (df_merged.probe_type == "reads_with_probe_no_polyT") | (df_merged.probe_type == "reads_with_probe_with_polyT")),
            "reads_with_polyT": sum( (df_merged.probe_type == "polyT") | (df_merged.probe_type == "reads_with_probe_with_polyT")),
            "reads_with_others": sum(df_merged.probe_type == "others"),
            "reads_with_probe_no_polyT": sum(df_merged.probe_type == "reads_with_probe_no_polyT"),
            "reads_with_probe_with_polyT":sum(df_merged.probe_type == "reads_with_probe_with_polyT"),
            "reads_no_probe_with_polyT":sum(df_merged.probe_type == "polyT"),
            "reads_no_probe_no_polyT":sum(df_merged.probe_type == "others"),
            "S_reads":len(S_list),
            "X_reads":len(X_list),
            "pgRNA_reads":len(pgRNA_list),
            "rcDNA_reads":len(rcDNA_list),
            "cccDNA_reads":len(cccDNA_list),
            "polyT_reads":len(polyT_list),
            "others_reads":len(others_list)
            })
        with open(f"outdir/1.bam/summarize.txt", "w") as fd:
            for k, v in sdict.items():
                fd.write(f"{k}: {v}\n")

        # ref record df
        input_bam = pysam.AlignmentFile(glob.glob(f"{self.fj_path}/*star_virus*/*Aligned.sortedByCoord.out.bam")[0], "rb")
        ref_record_dict = defaultdict(list)
        for i in input_bam:
            ref_record_dict["qname"].append(i.qname)
            ref_record_dict["record"].append(i)
        ref_record_df = pd.DataFrame(ref_record_dict)
            #ref_record_df_filter = ref_record_df[ref_record_df.qname.isin(probe_list)]    # raw
        ref_record_df['barcode_umi'] = ref_record_df['qname'].str.split("_").str[0] + "_" + ref_record_df['qname'].str.split("_").str[1]
        input_bam.close()

        # out bam
        def probe_bam_out(probe_list, out_bam_name):
            input_bam = pysam.AlignmentFile(glob.glob(f"{self.fj_path}/*star_virus*/*Aligned.sortedByCoord.out.bam")[0], "rb")
            out_bam = pysam.AlignmentFile(out_bam_name, "wb", template = input_bam)

            ref_record_df_filter = ref_record_df[ref_record_df.qname.isin(probe_list)]
            for _, row in ref_record_df_filter.iterrows():
                out_bam.write(row["record"])  

            input_bam.close()
            out_bam.close()

        S_out_bam = f'outdir/1.bam/S.bam'
        X_out_bam = f'outdir/1.bam/X.bam'
        pgRNA_out_bam = f'outdir/1.bam/pgRNA.bam'
        rcDNA_out_bam = f'outdir/1.bam/rcDNA.bam'
        cccDNA_out_bam = f'outdir/1.bam/cccDNA.bam'
        polyT_out_bam = f'outdir/1.bam/polyT.bam'
        others_out_bam = f'outdir/1.bam/others.bam'

        probe_bam_out(S_list, S_out_bam)
        probe_bam_out(X_list, X_out_bam)
        probe_bam_out(pgRNA_list, pgRNA_out_bam)
        probe_bam_out(rcDNA_list, rcDNA_out_bam)
        probe_bam_out(cccDNA_list, cccDNA_out_bam)
        probe_bam_out(polyT_list, polyT_out_bam)
        probe_bam_out(others_list, others_out_bam)
  
        # done
        print("----")
        print("probe split bam done.")


    def bam_to_fq(self):
        cmd1 = f"/SGRNJ/Public/Software/conda_env/js_env/bin/bedtools bamtofastq -i outdir/1.bam/S.bam -fq outdir/2.fastq/S.fq"
        cmd2 = f"/SGRNJ/Public/Software/conda_env/js_env/bin/bedtools bamtofastq -i outdir/1.bam/X.bam -fq outdir/2.fastq/X.fq"
        cmd3 = f"/SGRNJ/Public/Software/conda_env/js_env/bin/bedtools bamtofastq -i outdir/1.bam/pgRNA.bam -fq outdir/2.fastq/pgRNA.fq"
        cmd4 = f"/SGRNJ/Public/Software/conda_env/js_env/bin/bedtools bamtofastq -i outdir/1.bam/rcDNA.bam -fq outdir/2.fastq/rcDNA.fq"
        cmd5 = f"/SGRNJ/Public/Software/conda_env/js_env/bin/bedtools bamtofastq -i outdir/1.bam/cccDNA.bam -fq outdir/2.fastq/cccDNA.fq"
        
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)
        subprocess.check_call(cmd5, shell = True)

    def run(self):
        self.mkdir()
        self.split_bam()
        self.bam_to_fq()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--fj_path', help='fj path', required=True)
    parsers.add_argument('--probe_anno', help='probe anno', required=True)

    args = parsers.parse_args()
    runner = Probe_split(args) 
    runner.run()

if __name__ == '__main__':
    main()
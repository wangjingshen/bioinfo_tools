import os
import glob
import argparse
import pandas as pd
import scanpy as sc

class Hbv_virus_polyT:
    def __init__(self, args):
        self.fj_path = args.fj_path
        self.probe_anno = args.probe_anno
        self.outdir = args.outdir

        if not os.path.exists(self.outdir):
            os.system(f"mkdir -p {self.outdir}")

    def analysis(self):
        # read virus matrix -- 
        virus_mtx = sc.read_10x_mtx(glob.glob(f"{self.fj_path}/08.count_capture_virus_mtx/*_virus_matrix/")[0])
        virus_mtx_df = pd.DataFrame(virus_mtx.X.todense())
        virus_mtx_df.columns = virus_mtx.var_names
        virus_mtx_df.index = virus_mtx.obs_names

        # read 08 count detail --
        count_detail = pd.read_csv(glob.glob(f"{self.fj_path}/08.count_capture_virus_mtx/*_count_detail.txt")[0], sep="\t", header=0)
        count_detail = count_detail.loc[count_detail.Barcode.isin(virus_mtx_df.index),]  # filter HBV- barcode
        count_detail['barcode_umi'] = count_detail['Barcode'] + "_" + count_detail['UMI'] 
        # read probe anno --
        probe_anno = pd.read_csv(self.probe_anno, sep="\t", header = 0)
        probe_anno['barcode_umi'] = probe_anno['barcode'] + "_" + probe_anno['UMI']    
        probe_anno_virus = probe_anno.loc[ probe_anno.barcode_umi.isin(count_detail.barcode_umi),]

        def sum_probe(x):
            return(sum(x != "no_probe"))     # sum !no_probe == sum probe; 5 probe name, such as S,X
        def sum_polyT(x):
            return(sum(x == "polyT"))
        probe_anno_virus_agg = probe_anno_virus.groupby(['barcode','UMI','barcode_umi']).agg({
            'probe_name': sum_probe,
            'is_polyT': sum_polyT
        }).reset_index()

        probe_anno_virus_agg['polyT_umi'] = 0
        probe_anno_virus_agg.loc[(probe_anno_virus_agg.probe_name==0) & (probe_anno_virus_agg.is_polyT>0), 'polyT_umi'] = 1   # all reads no probe and at least one reads polyT

        probe_anno_virus_agg_file =  f"{self.outdir}/probe_anno_virus.tsv"
        probe_anno_virus_agg.to_csv(probe_anno_virus_agg_file, sep="\t", index=False)

        barcode_polyT_umi = probe_anno_virus_agg.groupby(['barcode']).agg({
            'polyT_umi': 'sum'
        })
        virus_mtx_polyT = pd.merge(virus_mtx_df, barcode_polyT_umi, left_index=True, right_index=True).reset_index()
        virus_mtx_polyT = virus_mtx_polyT.rename(columns={'index':'Barcode'}) # rename
        virus_mtx_polyT_file = f"{self.outdir}/virus_mtx_polyT.tsv"
        virus_mtx_polyT.to_csv(virus_mtx_polyT_file, sep="\t", index=False)

        # add gene to probe anno --
        probe_anno_virus_gene = pd.merge(probe_anno_virus_agg,count_detail,on="barcode_umi")
        probe_anno_virus_gene_polyT = probe_anno_virus_gene.groupby(['Barcode','geneID']).agg({'polyT_umi' : 'sum'}).reset_index().pivot(index="Barcode", columns="geneID", values='polyT_umi').fillna(0)
        probe_anno_virus_gene_polyT.columns = 'polyT_umi_from_' + probe_anno_virus_gene_polyT.columns
        
        virus_mtx_polyT_gene = pd.merge(virus_mtx_polyT, probe_anno_virus_gene_polyT, on="Barcode")
        virus_mtx_polyT_gene_file =  f"{self.outdir}/virus_mtx_polyT_gene.tsv"
        virus_mtx_polyT_gene.to_csv(virus_mtx_polyT_gene_file, sep="\t", index=False)

    def run(self):
        self.analysis()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--fj_path', help='fj celescope path', required=True)
    parsers.add_argument('--probe_anno', help='probe anno', required=True)
    parsers.add_argument('--outdir', help='dir to save result file', required=True)

    args = parsers.parse_args()
    runner = Hbv_virus_polyT(args) 
    runner.run()

if __name__ == '__main__':
    main()
#-*-coding:utf-8-*-

import pandas as pd
import argparse

# args ----
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--counts', type=str, help = 'counts file', default = None,required=True)
parser.add_argument('--count_detail', type=str, help = 'count_tail file', default = None,required=True)
parser.add_argument('--outdir',help='output dir file', default = None,required=True)
args = parser.parse_args()

# read file --
count_detail = pd.read_table(args.count_detail,sep="\t",header=0)
counts = pd.read_table(args.counts,sep="\t",header = 0)

# set mark for count_detail--
cell_bc = counts.loc[counts['mark']=="CB","Barcode"]
count_detail["mark"]="UB"
count_detail.loc[count_detail["Barcode"].isin(cell_bc),"mark"]="CB"

# calc fraction_reads_in_cells --
def calc(df):
    return (df.loc[df["mark"]=="CB","count"].sum())/(df["count"].sum())

count_fraction = count_detail.groupby('geneID').apply(calc)
count_fraction = count_fraction.to_frame().reset_index()
count_fraction.columns = ["geneID","fraction_reads_in_cells"]

# output --
count_fraction.to_csv(args.outdir + "/gene_fraction_reads_in_cell.tsv",index=0,sep="\t")
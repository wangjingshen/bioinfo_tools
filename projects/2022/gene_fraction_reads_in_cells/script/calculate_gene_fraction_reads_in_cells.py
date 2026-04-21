#-*-coding:utf-8-*-
import os
import pandas as pd
import argparse

class Gene_fraction():
    def __init__(self,args):
        self.counts = args.counts
        self.count_detail = args.count_detail
        self.outdir = args.outdir
    
    def calculate_fraction(self):
        counts_file = pd.read_table(self.counts, sep="\t", header = 0)
        count_detail_file = pd.read_table(self.count_detail, sep="\t", header=0)
        cell_bc = counts_file.loc[counts_file['mark'] == 'CB', 'Barcode'] 
        count_detail_file['mark'] = 'UB'     # set mark for count_detail
        count_detail_file.loc[count_detail_file['Barcode'].isin(cell_bc), 'mark'] = 'CB'    # set CB
 
        calculate_cb_fraction = lambda df : round(float(df.loc[df['mark'] == 'CB','count'].sum())/(df['count'].sum()),4)  # calculate fraction
        count_fraction = count_detail_file.groupby('geneID').apply(calculate_cb_fraction).to_frame().reset_index()
        count_fraction.columns = ['geneID', 'fraction_reads_in_cells']
        count_fraction.to_csv(self.outdir + "/gene_fraction_reads_in_cells.tsv", index = 0, sep = "\t")  #output
    
    def run(self):
        self.calculate_fraction()

def main():
    parser = argparse.ArgumentParser(description = 'manual to this script')
    parser.add_argument('--counts', type = str, help = '05.count/xxx_counts.txt, a output of CeleScope', default = None, required = True)
    parser.add_argument('--count_detail', type = str, help = '05.count/xxx_count_detail.txt, a output of CeleScope', default = None, required = True)
    parser.add_argument('--outdir', help = 'output dir', default = None, required = True)
    args = parser.parse_args()
    gene_fraction = Gene_fraction(args) 
    gene_fraction.run()

if __name__ == "__main__":
    main()
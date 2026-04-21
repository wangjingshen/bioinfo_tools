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
faulthandler.enable()  # debug
from Bio.Align.Applications import MuscleCommandline
from io import StringIO
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from xopen import xopen
from itertools import groupby
import sys
sys.path.append('/SGRNJ06/randd/USER/wangjingshen/script/tools')

import utils as utils



class Virus_gene_consensus:
    def __init__(self, args):
        self.fj_path = args.fj_path
        self.max_n_reads = int(args.max_n_reads)
        self.small_threshold = args.small_threshold
        self.large_threshold = args.large_threshold
        self.min_consensus_read = int(args.min_consensus_read)
        self.outdir = args.outdir
    
    def mkdir(self):
        if not os.path.exists(f"{self.outdir}/1.bam/"):
            os.system(f"mkdir -p {self.outdir}/1.bam/")
        if not os.path.exists(f"{self.outdir}/2.fasta/"):
            os.system(f"mkdir -p {self.outdir}/2.fasta/")
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


    def bam_to_fa(self):
        '''
        bam to fa for consensus
        '''
        def bam_to_fa_run(gene):
            #cmd = f'''samtools sort -n {self.outdir}/1.bam/{gene}.bam | samtools view -O SAM | awk -F '\t' '{print ">"$1"\n"$10}' > {self.outdir}/2.fasta/{gene}.fa'''
            cmd = (f'samtools sort -n {self.outdir}/1.bam/{gene}.bam | samtools view -O SAM | '
                     r''' awk -F '\t' '{print ">"$1"\n"$10}' > '''
                     f'{self.outdir}/2.fasta/{gene}.fa'
                     )
            subprocess.check_call(cmd, shell = True)
        bam_to_fa_run("S")
        bam_to_fa_run("X")
        bam_to_fa_run("pgRNA")
        bam_to_fa_run("rcDNA")
        bam_to_fa_run("cccDNA")
    
    @utils.add_log
    def run_consensus(self):       
        def keyfunc(read):
            attr = read.name.split('_')
            return (attr[0], attr[1])

        def consensus(gene):
            n_umi = 0
            out_h = xopen(f"{self.outdir}/3.consensus/{gene}_consensus.fa", 'w')

            if(os.path.getsize(f"{self.outdir}/2.fasta/{gene}.fa")==0):
                Virus_gene_consensus.run_consensus.logger.warning(f"{self.outdir} - {gene}.fa is null, skiping")
                out_h.close()
                return
            if(os.path.getsize(f"{self.outdir}/2.fasta/{gene}.fa")!=0):
                with pysam.FastxFile(f"{self.outdir}/2.fasta/{gene}.fa") as fh:
                    for (barcode, umi), g in groupby(fh, key=keyfunc):
                        read_list_s = []
                        read_list_l = []
                        for read in g:
                            read_list_l.append([read.sequence])
                            read_list_s.append(">"+read.name+"\n"+read.sequence+"\n")

                        #print(len(read_list_s))   # for check
                        if(len(read_list_s)<= self.max_n_reads):
                            '''
                            For small read list, use muslce for multiple sequence alignment, then use dumb_consensus for consensus.
                            '''
                            read_list = "".join(read_list_s)
                            records = SeqIO.parse(StringIO(read_list), "fasta")
                            handle = StringIO()
                            SeqIO.write(records, handle, "fasta")
                            muscle_cline = MuscleCommandline()
                            stdout, stderr = muscle_cline(stdin= handle.getvalue())
                            align = AlignIO.read(StringIO(stdout), "fasta")
                            summary_align = AlignInfo.SummaryInfo(align)
                            consensus_seq =  str(summary_align.dumb_consensus(threshold= self.small_threshold,ambiguous='N'))
                        
                        if(len(read_list_l) > self.max_n_reads):
                            '''
                            For large read list, use celescope consensus to reduce computing burden (such as, 200k reads of 1 UMI).
                            '''
                            consensus_seq = dumb_consensus(read_list_l, threshold = self.large_threshold, min_consensus_read = self.min_consensus_read, ambiguous="N")
                
                        n_umi += 1
                        prefix = "_".join([barcode, umi])
                        read_name = f'{prefix}_{n_umi}'
                        out_h.write(f'>{read_name}\n{consensus_seq}\n')
                        if n_umi % 1000 == 0:
                            Virus_gene_consensus.run_consensus.logger.info(f'{self.outdir} - {gene}: {n_umi} UMI done.')
                    out_h.close()
                    Virus_gene_consensus.run_consensus.logger.info(f'{self.outdir} - {gene} done.')
        
        consensus("S")
        consensus("pgRNA")
        consensus("rcDNA")
        consensus("cccDNA")
        consensus("X")

    
    @utils.add_log
    def consensus_summary(self):
        def get_consensus_summary(gene):
            out_h = xopen(f"{self.outdir}/3.consensus/{gene}_consensus_summary.fa", 'w')
            
            if(os.path.getsize(f"{self.outdir}/3.consensus/{gene}_consensus.fa")==0):
                Virus_gene_consensus.consensus_summary.logger.warning(f"{self.outdir} - {gene}_consensus.fa is null, skiping")
                out_h.close()
                return
            if(os.path.getsize(f"{self.outdir}/3.consensus/{gene}_consensus.fa")!=0):
                record_dict = defaultdict(list)
                with pysam.FastxFile(f"{self.outdir}/3.consensus/{gene}_consensus.fa") as f:
                    for i in f:
                        record_dict['name'].append(i.name)
                        record_dict['sequence'].append(i.sequence)
                record_df = pd.DataFrame(record_dict)
                record_df_sum = record_df['sequence'].value_counts().to_frame().sort_values(by='sequence',axis=0,ascending=False).reset_index().rename(columns = {'index':'seq','sequence':'count'})
                record_df_sum = (
                    record_df_sum.assign(seq_id = lambda df:"variant-"+ (df.index + 1).astype('str') + " : "+ df['count'].astype('str'))
                )
                for i in range(len(record_df_sum)):
                    out_h.write(">" + record_df_sum.iloc[i,2]+ "\n" + record_df_sum.iloc[i,0]+ "\n")

                Virus_gene_consensus.consensus_summary.logger.info(f'{self.outdir} - {gene} consensus summary done.')
                        
        get_consensus_summary("S")
        get_consensus_summary("X")
        get_consensus_summary("pgRNA")
        get_consensus_summary("rcDNA")
        get_consensus_summary("cccDNA")
    
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
        self.bam_to_fa()
        self.run_consensus()
        self.consensus_summary()
        ##self.rm_and_gzip()   # for bak use
    

def dumb_consensus(read_list, threshold=0.5, min_consensus_read=1, ambiguous='N'):  #, default_qual='F'):
    '''
    This code comes from the celescope(for fastq). Quality related code is removed here(for fasta).
    ---
    This is similar to biopython dumb_consensus.
    It will just go through the sequence residue by residue and count up the number of each type
    of residue (ie. A or G or T or C for DNA) in all sequences in the
    alignment. If 
    1. the percentage of the most common residue type > threshold;
    2. most common residue reads >= min_consensus_read;
    then we will add that residue type,
    otherwise an ambiguous character will be added.
    elements of read_list: [entry.sequence,entry.quality]
    '''

    con_len = get_read_length(read_list, threshold=threshold)
    consensus_seq = ""
    #consensus_qual = ""
    ambiguous_base_n = 0
    for n in range(con_len):
        atom_dict = defaultdict(int)
        #quality_dict = defaultdict(int)
        num_atoms = 0
        for read in read_list:
            # make sure we haven't run past the end of any sequences
            # if they are of different lengths
            sequence = read[0]
            #quality = read[1]
            if n < len(sequence):
                atom = sequence[n]
                atom_dict[atom] += 1
                num_atoms = num_atoms + 1

                #base_qual = quality[n]
                #quality_dict[base_qual] += 1

        consensus_atom = ambiguous
        for atom in atom_dict:
            if atom_dict[atom] > num_atoms * threshold and atom_dict[atom] >= min_consensus_read:
                consensus_atom = atom
                break
        if consensus_atom == ambiguous:
            ambiguous_base_n += 1
        consensus_seq += consensus_atom

        #max_freq_qual = 0
        #consensus_base_qual = default_qual
        #for base_qual in quality_dict:
        #    if quality_dict[base_qual] > max_freq_qual:
        #        max_freq_qual = quality_dict[base_qual]
        #        consensus_base_qual = base_qual

        #consensus_qual += consensus_base_qual
    return consensus_seq #, consensus_qual, ambiguous_base_n, con_len

def get_read_length(read_list, threshold=0.5):
    '''
    compute read_length from read_list. 
    length = max length with read fraction >= threshold
    elements of read_list: [entry.sequence,entry.quality]
    '''

    n_read = len(read_list)
    length_dict = defaultdict(int)
    for read in read_list:
        length = len(read[0])
        length_dict[length] += 1
    for length in length_dict:
        length_dict[length] = length_dict[length] / n_read

    fraction = 0
    for length in sorted(length_dict.keys(), reverse=True):
        fraction += length_dict[length]
        if fraction >= threshold:
            return length



def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--fj_path', help='fj path', required=True)
    parsers.add_argument('--max_n_reads', help='max n reads to run biopython consensus', required=True)
    parsers.add_argument('--small_threshold', help='valid base threshold of small UMI', type=float, default=0.5)
    parsers.add_argument('--large_threshold', help='valid base threshold of large UMI', type=float, default=0.5)
    parsers.add_argument('--min_consensus_read', help='Minimum number of reads to support a base', default=1)

    parsers.add_argument('--outdir', help='outdir', required=True)


    args = parsers.parse_args()
    runner = Virus_gene_consensus(args) 
    runner.run()

if __name__ == '__main__':
    main()
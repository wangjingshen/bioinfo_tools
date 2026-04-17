import os
import pysam
import pandas as pd
import subprocess
import argparse
from collections import defaultdict


class Linker_stat:
    def __init__(self, args):
        self.R1 = args.R1
        self.mismatch = int(args.mismatch)
        self.outdir = args.outdir
        self.name = args.name

        if not os.path.exists(f"{self.outdir}"):
            os.system(f"mkdir -p {self.outdir}")

    
    def run(self):
        '''
        /SGRNJ/Public/Software/conda_env/celescope1.9.0/lib/python3.9/site-packages/celescope/data/chemistry/scopeV3.0.1/linker_4types
        linker seq1: ATCCACGTGCTTGAGATCAGCATGCGGCTACGC
        linker seq2: TCGGTGACAGCCATATCGTAGTCAGAAGCTGAC
        linker seq3: CGAACATGTAGGTCTCGACTACGTATTAGCATC
        linker seq4: GATTGTCACTAACGCGATGCTGACTCCTAGTCC
        '''
        record_dict = defaultdict(list)
        with pysam.FastxFile(f'{self.R1}') as f:
            for i in f:
                record_dict["linker1"].append(i.sequence[9:25])
                record_dict["linker2"].append(i.sequence[34:50])
                record_dict["linker3"].append(i.sequence[59:60])
        record_df = pd.DataFrame(record_dict)
        record_df["linker"] = record_df.linker1 + record_df.linker2 + record_df.linker3
        record_df['hamming_seq1'] = record_df.linker.apply(hamming_correct, string2 = 'ATCCACGTGCTTGAGATCAGCATGCGGCTACGC', threshold = self.mismatch + 1)
        record_df['hamming_seq2'] = record_df.linker.apply(hamming_correct, string2 = 'TCGGTGACAGCCATATCGTAGTCAGAAGCTGAC', threshold = self.mismatch + 1)
        record_df['hamming_seq3'] = record_df.linker.apply(hamming_correct, string2 = 'CGAACATGTAGGTCTCGACTACGTATTAGCATC', threshold = self.mismatch + 1)
        record_df['hamming_seq4'] = record_df.linker.apply(hamming_correct, string2 = 'GATTGTCACTAACGCGATGCTGACTCCTAGTCC', threshold = self.mismatch + 1)

        res_df =pd.DataFrame({'ATCCACGTGCTTGAGATCAGCATGCGGCTACGC': record_df['hamming_seq1'].sum(),
                              'TCGGTGACAGCCATATCGTAGTCAGAAGCTGAC': record_df['hamming_seq2'].sum(),
                              'CGAACATGTAGGTCTCGACTACGTATTAGCATC': record_df['hamming_seq3'].sum(),
                              'GATTGTCACTAACGCGATGCTGACTCCTAGTCC': record_df['hamming_seq4'].sum()},index=[0]).T.reset_index()
        res_df.columns = ['linker', 'counts']
        res_df['percent'] = res_df.counts/record_df.shape[0]

        res_file = f'{self.outdir}/{self.name}_linker_stat.tsv'
        res_df.to_csv(res_file, sep = "\t", index = 0)


def hamming_correct(string1, string2, threshold):
    #threshold = len(string1) / 10 + 1
    if hamming_distance(string1, string2) < threshold:
        return True
    return False

def hamming_distance(string1, string2):
    distance = 0
    length = len(string1)
    length2 = len(string2)
    if (length != length2):
        raise Exception(f"string1({length}) and string2({length2}) do not have same length")
    for i in range(length):
        if string1[i] != string2[i]:
            distance += 1
    return distance

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--R1', help='fj path', required=True)
    parsers.add_argument('--mismatch', help='mismatch', required=True)
    parsers.add_argument('--outdir', help='outdir', required=True)
    parsers.add_argument('--name', help='name', required=True)
    args = parsers.parse_args()
    runner = Linker_stat(args) 
    runner.run()

if __name__ == '__main__':
    main()
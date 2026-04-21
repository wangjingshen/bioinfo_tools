#!/usr/bin/env python

import subprocess
import json
import argparse
import logging

logger = logging.getLogger(__name__)


from step import Step, s_common
#from celescope.tools import utils

POLY_A = '-a polyA=A{18} '

LOG_METRICS_TITLE = (
    'Total reads processed',
    'Reads with adapters',
    'Reads that were too short',
    'Reads written (passing filters)',
    'Total basepairs processed',
    'Quality-trimmed',
    'Total written (filtered)',
)


def get_cutadapt_cmd(args, input_file, output_file):
    cmd = (
        'cutadapt '
        f'{POLY_A} '
        f'-j {args.thread} '
        f'-m {args.minimum_length} '
        f'--nextseq-trim={args.nextseq_trim} '
        f'--overlap {args.overlap} '
        f'{args.cutadapt_param} '
        f'--json {args.sample}_cutadapt.json '
        f'-o {output_file} '
        f'{input_file} '
        )
    return cmd

class Cutadapt(Step):
    """
    ## Features
    - Trim poly A tails and user-provided adapters in R2 reads with cutadapt.

    ## Output
    - `cutadapt.log` Cutadapt output log file.
    - `{sample}_clean_2.fq.gz` R2 reads file without adapters.
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        # out files
        self.out_fq2 = f'{self.sample}_clean_2.fq'
        self.json_log = f'{self.sample}_cutadapt.json'


    def add_cutadapt_metrics(self):
        with open(self.json_log) as f:
            log_dict = json.load(f)

        total_reads = log_dict['read_counts']['input']
        reads_with_adapters = log_dict['read_counts']['read1_with_adapter']
        reads_too_short = log_dict['read_counts']['filtered']['too_short']
        total_base_pairs = log_dict['basepair_counts']['input']
        quality_trimmed = log_dict['basepair_counts']['quality_trimmed']

        self.add_metric(
            name='Reads with Adapters',
            value=reads_with_adapters,
            total=total_reads,
            help_info='reads with poly A tails and user-provided adapters(if any) are trimmed'
        )
        self.add_metric(
            name='Reads too Short',
            value=reads_too_short,
            total=total_reads,
            help_info=f'reads with read length less than {self.args.minimum_length} bp after trimming'
        )
        self.add_metric(
            name='Base Pairs Quality-Trimmed',
            value=quality_trimmed,
            total=total_base_pairs,
            help_info='bases pairs removed from the end of the read whose quality is smaller than the given threshold'
        )

    def run(self):
        input_file = self.args.fq
        output_file = self.out_fq2
        cmd = get_cutadapt_cmd(self.args, input_file, output_file)
        #logger.info(cmd)  #self.run.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
        self.add_cutadapt_metrics()


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--fq', help='input', required=True)
    parser.add_argument('--minimum_length',
                        help='Discard processed reads that are shorter than LENGTH.',
                        default=20)    
    parser.add_argument('--nextseq_trim',
                        help="""Quality trimming of reads using two-color chemistry (NextSeq). 
                                Some Illumina instruments use a two-color chemistry to encode the four bases. 
                                This includes the NextSeq and the NovaSeq. 
                                In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. 
                                However, dark cycles also occur when sequencing “falls off” the end of the fragment.
                                The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.""",
                        default=20,)    
    parser.add_argument('--overlap',
                        help="""Since Cutadapt allows partial matches between the read and the adapter sequence,
                                short matches can occur by chance, leading to erroneously trimmed bases. 
                                For example, roughly 0.25 of all reads end with a base that is identical to the first base of the adapter. 
                                To reduce the number of falsely trimmed bases, the alignment algorithm requires that 
                                at least {overlap} bases match between adapter and read. """, 
                        default=10)    
    parser.add_argument('--cutadapt_param', 
                        help='Other cutadapt parameters. For example, --cutadapt_param "-a p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" ', 
                        default="")


    # add version
    parser.add_argument('--version', action='version', version='1.0')
    parser = s_common(parser)  #
    args = parser.parse_args()

    with Cutadapt(args, display_title="Trimming") as runner:
        runner.run()
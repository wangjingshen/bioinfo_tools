#!/usr/bin/env python

import argparse
import json
import itertools
import math
import os
import glob
import re
import sys
from collections import Counter, defaultdict
from itertools import combinations, product
import logging

import pysam
import pyfastx
from xopen import xopen

import utils
from step import Step, s_common

logger = logging.getLogger(__name__)

MIN_T = 10


def get_seq_str(seq, sub_pattern):
    """
    join seq slices.

    Args:
        seq: usually R1 read
        sub_pattern: [slice(0,8),slice(16,24)]

    Returns:
        joined intervals seq

    Raises:
        IndexError: if seq length is not enough.

    >>> sub_pattern_dict = [slice(0,2)]
    >>> seq = "A" * 2 + "T" * 2
    >>> get_seq_str(seq, sub_pattern_dict)
    'AA'
    """
    seq_len = len(seq)
    expect_len = sub_pattern[-1].stop
    if seq_len < expect_len:
        raise IndexError(f"read length({seq_len} bp) less than expected length({expect_len} bp) in read: {seq}")
    return "".join([seq[x] for x in sub_pattern])


def findall_mismatch(seq, n_mismatch=1, bases='ACGTN'):
    """
    choose locations where there's going to be a mismatch using combinations
    and then construct all satisfying lists using product

    Return:
    all mismatch <= n_mismatch set.

    >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
    >>> seq_set = findall_mismatch("ACG")
    >>> seq_set == answer
    True
    """
    seq_set = set()
    seq_len = len(seq)
    if n_mismatch > seq_len:
        n_mismatch = seq_len
    for locs in itertools.combinations(range(seq_len), n_mismatch):
        seq_locs = [[base] for base in seq]
        for loc in locs:
            seq_locs[loc] = list(bases)
        for poss in itertools.product(*seq_locs):
            seq_set.add(''.join(poss))
    return seq_set


def get_mismatch_dict(seq_list, n_mismatch=1):
    """
    Return:
    mismatch dict. Key: mismatch seq, value: seq in seq_list

    >>> seq_list = ["AACGTGAT", "AAACATCG"]
    >>> mismatch_dict = get_mismatch_dict(seq_list)
    >>> mismatch_dict["AACGTGAA"] == "AACGTGAT"
    True
    """
    mismatch_dict = {}
    for seq in seq_list:
        seq = seq.strip()
        if seq == '':
            continue
        for mismatch_seq in findall_mismatch(seq, n_mismatch):
            mismatch_dict[mismatch_seq] = seq
    return mismatch_dict


def read_one_col(fn):
    """read one column file into list"""
    with open(fn) as f:
        return [x.strip() for x in f]


def parse_pattern(pattern, allowed="CLUNT"):
    """
    >>> pattern_dict = parse_pattern("C8L16C8L16C8L1U12T18")
    >>> pattern_dict['C']
    [slice(0, 8, None), slice(24, 32, None), slice(48, 56, None)]
    >>> pattern_dict['L']
    [slice(8, 24, None), slice(32, 48, None), slice(56, 57, None)]
    """
    pattern_dict = {}
    p = re.compile(r'([A-Z])(\d+)')
    tmp = p.findall(pattern)
    if not tmp:
        sys.exit(f'Invalid pattern: {pattern}')
    start = 0
    for x, length in tmp:
        if x not in allowed:
            sys.exit(f'Invalid pattern: {pattern}')
        if x not in pattern_dict:
            pattern_dict[x] = []
        end = start + int(length)
        pattern_dict[x].append(slice(start,end))
        start = end
    return pattern_dict


def get_protocol_dict(assets_dir):
    """
    Return:
    protocol_dict. Key: protocol name, value: protocol dict

    >>> protocol_dict = get_protocol_dict("./assets/")
    >>> protocol_dict["GEXSCOPE-MicroBead"]["pattern_dict"]
    {'C': [slice(0, 12, None)], 'U': [slice(12, 20, None)]}
    """
    json_file = os.path.join(assets_dir, "protocols.json")
    protocol_dict = json.load(open(json_file))
    whitelist_dir = os.path.join(assets_dir, "whitelist")
    # add folder prefix
    for protocol in protocol_dict:
        cur = protocol_dict[protocol]
        bc = cur.get("bc", [])
        linker = cur.get("linker", [])
        if bc:
            cur["bc"] = [os.path.join(whitelist_dir, protocol, x) for x in bc]
        if linker:
            cur["linker"] = [os.path.join(whitelist_dir, protocol, x) for x in linker]
        cur["pattern_dict"] = parse_pattern(cur["pattern"])
    return protocol_dict


def get_seq_list(seq, pattern_dict, abbr):
    """
    >>> pattern_dict = parse_pattern("C2L3C2")
    >>> seq = "AAGGGTT"
    >>> get_seq_list(seq, pattern_dict, "C")
    ['AA', 'TT']
    """    
    return [seq[item.start: item.stop] for item in pattern_dict[abbr]]    # slice.start; stop


def get_abbr_len(pattern_dict, abbr):
    """
    >>> pattern_dict = parse_pattern("C8L16C8L16C8L1U12T18")
    >>> get_abbr_len(pattern_dict, 'C')
    24
    >>> get_abbr_len(pattern_dict, 'L')
    33
    """
    length = 0
    for item in pattern_dict[abbr]:
        length += item.stop - item.start

    return length


def get_seq_str_no_exception(seq, sub_pattern_dict):
    """get subseq with intervals in arr and concatenate"""
    return ''.join([seq[item.start: item.stop] for item in sub_pattern_dict])


class Protocol:
    """
    Auto detect singleron protocols from R1-read
    GEXSCOPE-MicroBead
    GEXSCOPE-V1
    GEXSCOPE-V2
    sCircle-V1
    """

    def __init__(self, fq1_list, sample, assets_dir='assets/', max_read=10000):
        '''
        Args:
            assets_dir: Expects file 'protocols.json' and 'whitelist/{protocol}' folder under assets_dir
        '''
        self.fq1_list = fq1_list
        self.max_read = max_read
        self.sample = sample
        self.protocol_dict = get_protocol_dict(assets_dir)
        self.scircle_R1_LEN = self.protocol_dict["GEXSCOPE-sCircle"]["pattern_dict"]["U"][-1].stop


    def get_protocol(self):
        """check protocol in the fq1_list"""
        fq_protocol = {}
        for fastq1 in self.fq1_list:
            protocol = self.get_fq_protocol(fastq1)
            fq_protocol[fastq1] = protocol
        if len(set(fq_protocol.values())) != 1:
            sys.exit(f'Error: multiple protocols are not allowed for one sample: {self.sample}! \n' + str(fq_protocol))
        protocol = list(fq_protocol.values())[0]
        return protocol


    def check_linker(self, seq, protocol, i):
        """check if seq matches the linker i of protocol"""
        linker_fn = self.protocol_dict[protocol]["linker"]
        linker = read_one_col(linker_fn[i-1])
        linker_mismatch_dict = get_mismatch_dict(linker)

        pattern = self.protocol_dict[protocol]["pattern"]
        pattern_dict = parse_pattern(pattern)
        seq_linker = seq[pattern_dict["L"][i-1]]

        return seq_linker in linker_mismatch_dict


    def seq_protocol(self, seq):
        """
        Returns: protocol or None

        >>> runner = Protocol("fake_fq1_string", "fake_sample")
        >>> seq = "TCGACTGTC" + "ATCCACGTGCTTGAGA" + "TTCTAGGAT" + "TCAGCATGCGGCTACG" + "TGCACGAGA" + "C" + "CATATCAATGGG" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-V2'
        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA" + "C" + "TCCGAAGCCCAT" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-V1'
        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA"  + "TCCGAAGCC" + "CTGTCT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-sCircle'
        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA"  + "TCCGAAGCC"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-sCircle'
        >>> seq = "ATCGATCGATCG" + "ATCGATCG" + "C" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-MicroBead'
        """
        # check "GEXSCOPE-sCircle" first as it is similar to "GEXSCOPE-V1"
        protocol = "GEXSCOPE-sCircle"
        if self.check_linker(seq, protocol, 1) and (len(seq) == self.scircle_R1_LEN or self.check_linker(seq, protocol, 3)):
            return protocol

        for protocol in ["GEXSCOPE-V1", "GEXSCOPE-V2"]:
            if self.check_linker(seq, protocol, 1):
                return protocol

        # check if it is MicroBead
        if seq[16:20] != "TTTT" and seq[22:26] == "TTTT":
            return "GEXSCOPE-MicroBead"


    def get_fq_protocol(self, fq1):
        results = defaultdict(int)

        fq = pyfastx.Fastx(fq1)
        n = 0
        for name,seq,qual in fq:
            n += 1
            protocol = self.seq_protocol(seq)
            if protocol:
                results[protocol] += 1
            if n == self.max_read:
                break
        sorted_counts = sorted(results.items(), key=lambda x: x[1], reverse=True)
        logger.info(sorted_counts)

        protocol, read_counts = sorted_counts[0]
        percent = float(read_counts) / n
        if percent < 0.5:
            logger.warning("Valid protocol read counts percent < 0.5")
        if percent < 0.1:
            logger.error("Valid protocol read counts percent < 0.1")
            raise Exception(
                'Auto protocol detection failed! '
            )
        logger.info(f'{fq1}: {protocol}')

        return protocol


class Barcode(Step):
    """
    ## Features

    - Demultiplex barcodes.
    - Filter invalid R1 reads, which includes:
        - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.  
        - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.  
        - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
        - Low quality reads: low sequencing quality in barcode and UMI regions.

    ## Output

    - `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
    the read name is `{barcode}_{UMI}_{read ID}`.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")
        self.fq_number = len(self.fq1_list)
        if self.fq_number != len(self.fq2_list):
            raise Exception('fastq1 and fastq2 do not have same file number!')

        self.assets_dir = args.assets_dir
        self.assay = args.assay
        if args.protocol == 'auto':
            protocol = Protocol(self.fq1_list, args.sample, assets_dir = self.assets_dir).get_protocol()
        else:
            protocol = args.protocol
        protocol_dict = get_protocol_dict(args.assets_dir)
        
        self.protocol_dict = protocol_dict
        self.protocol = protocol

        if protocol == 'new':
            self.pattern = args.pattern
            self.whitelist_file = args.whitelist
            self.linker_file = args.linker
        else:
            self.pattern = protocol_dict[protocol]["pattern"]
            self.whitelist_file = protocol_dict[protocol]["bc"]
            self.linker_file = protocol_dict[protocol]["linker"]

        self.barcode_corrected_num = 0
        self.linker_corrected_num = 0
        self.total_num = 0
        self.clean_num = 0
        self.no_polyT_num = 0
        self.lowQual_num = 0
        self.no_linker_num = 0
        self.no_barcode_num = 0
        self.barcode_qual_Counter = Counter()
        self.umi_qual_Counter = Counter()
        
        self.lowNum = args.lowNum
        self.lowQual = args.lowQual
        self.filterNoPolyT = args.filterNoPolyT
        self.allowNoLinker = args.allowNoLinker
        self.nopolyT = args.nopolyT  # true == output nopolyT reads
        self.noLinker = args.noLinker
        self.output_R1 = args.output_R1
        self.bool_flv = False

        # flv_trust4, flv_CR
        if self.assay in ('flv_CR', 'flv_trust4'):
            self.bool_flv = True
            self.barcode_read_Counter = Counter()
            if self.assay == 'flv_trust4':
                if args.match_dir == 'None':
                    raise FileNotFoundError('Match_dir required when running flv_trust4')
                self.match_barcodes = set(utils.get_barcode_from_match_dir(args.match_dir)[0]) # barcode set of flv_rna.
                self.match_num = 0 # record read number match with flv_rna.
                self.match_cbs = set() # record barcode number match with flv_rna.

        # out file
        self.out_fq2 = f'{self.sample}_2.fq'
        self.out_fq1 = f'{self.sample}_1.fq'
        if self.nopolyT:
            self.nopolyT_1 = f'{self.sample}_noPolyT_1.fq'
            self.nopolyT_2 = f'{self.sample}_noPolyT_2.fq'
        if self.noLinker:
            self.noLinker_1 = f'{self.sample}_noLinker_1.fq'
            self.noLinker_2 = f'{self.sample}_noLinker_2.fq'

        self.open_files()

    @staticmethod
    def ord2chr(q, offset=33):
        return chr(int(q) + offset)

    @staticmethod
    def qual_int(char, offset=33):
        return ord(char) - offset

    @staticmethod
    def low_qual(quals, minQ, num):
        # print(ord('/')-33)           14
        return True if len([q for q in quals if Barcode.qual_int(q) < minQ]) > num else False

    @staticmethod
    def check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list):
        '''
        Return bool_valid, bool_corrected, corrected_seq

        >>> seq_list = ['ATA', 'AAT', 'ATA']
        >>> correct_set_list = [{'AAA'},{'AAA'},{'AAA'}]
        >>> mismatch_dict_list = [Barcode.get_mismatch_dict(['AAA'])] * 3

        >>> Barcode.check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
        (True, True, 'AAA_AAA_AAA')

        >>> seq_list = ['AAA', 'AAA', 'AAA']
        >>> Barcode.check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
        (True, False, 'AAA_AAA_AAA')
        '''
        bool_valid = True
        bool_corrected = False
        corrected_seq_list = []
        for index, seq in enumerate(seq_list):
            if seq not in correct_set_list[index]:
                if seq not in mismatch_dict_list[index]:
                    bool_valid = False
                    return bool_valid, bool_corrected, ""
                else:
                    bool_corrected = True
                    corrected_seq_list.append(mismatch_dict_list[index][seq])
            else:
                corrected_seq_list.append(seq)

        return bool_valid, bool_corrected, '_'.join(corrected_seq_list)

    @staticmethod
    def parse_whitelist_file(files: list, n_pattern: int, n_mismatch: int):
        """
        files: file paths
        n_pattern: number of sections in pattern
        n_mismatch: allowed number of mismatch bases
        Returns:
            white_set_list
            mismatch_list
        """
        n_files = len(files)
        if n_files == 1 and n_pattern > 1:
            files = [files[0]] * n_pattern
        elif n_files != n_pattern:
            sys.exit(f'number of whitelist files({n_files} files:{files}) != n_pattern({n_pattern})')
        
        white_set_list, mismatch_list = [], []
        for f in files:
            barcodes, _ = utils.read_one_col(f)
            white_set_list.append(set(barcodes))
            barcode_mismatch_dict = get_mismatch_dict(barcodes, n_mismatch)
            mismatch_list.append(barcode_mismatch_dict)

        return white_set_list, mismatch_list

    @staticmethod
    def check_polyT(seq, pattern_dict, min_polyT_count=MIN_T):
        """
        Return:
            True if polyT is found
        """
        seq_polyT = get_seq_str(seq, pattern_dict['T'])
        n_polyT_found = seq_polyT.count('T')
        if n_polyT_found >= min_polyT_count:
            return True
        return False

    def open_files(self):
        if self.output_R1 or self.bool_flv:
            self.fh_fq1 = xopen(self.out_fq1, 'w')
        if not self.args.stdout:
            self.fh_fq2 = xopen(self.out_fq2, 'w')

        if self.nopolyT:
            self.fh_nopolyT_fq1 = xopen(self.nopolyT_1, 'w')
            self.fh_nopolyT_fq2 = xopen(self.nopolyT_2, 'w')

        if self.noLinker:
            self.fh_nolinker_fq1 = xopen(self.noLinker_1, 'w')
            self.fh_nolinker_fq2 = xopen(self.noLinker_2, 'w')

    def close_files(self):
        if self.output_R1 or self.bool_flv:
            self.fh_fq1.close()
        if not self.args.stdout:
            self.fh_fq2.close()

        if self.nopolyT:
            self.fh_nopolyT_fq1.close()
            self.fh_nopolyT_fq2.close()
        
        if self.noLinker:
            self.fh_nolinker_fq1.close()
            self.fh_nolinker_fq2.close()

    def add_step_metrics(self):

        self.add_metric(
            name='Raw Reads',
            value=self.total_num,
            help_info='total reads from FASTQ files'
        )
        self.add_metric(
            name='Valid Reads',
            value=self.clean_num,
            total=self.total_num,
            help_info='reads pass filtering(filtered: reads without poly T, reads without linker, reads without correct barcode or low quality reads)'
        )

        BarcodesQ30 = sum([self.barcode_qual_Counter[k] for k in self.barcode_qual_Counter if k >= self.ord2chr(
            30)]) / float(sum(self.barcode_qual_Counter.values())) * 100
        BarcodesQ30 = round(BarcodesQ30, 2)
        BarcodesQ30_display = f'{BarcodesQ30}%'
        self.add_metric(
            name='Q30 of Barcodes',
            value=BarcodesQ30,
            display=BarcodesQ30_display,
            help_info='percent of barcode base pairs with quality scores over Q30',
        )

        UMIsQ30 = sum([self.umi_qual_Counter[k] for k in self.umi_qual_Counter if k >= self.ord2chr(
            30)]) / float(sum(self.umi_qual_Counter.values())) * 100
        UMIsQ30 = round(UMIsQ30, 2)
        UMIsQ30_display = f'{UMIsQ30}%'
        self.add_metric(
            name='Q30 of UMIs',
            value=UMIsQ30,
            display=UMIsQ30_display,
            help_info='percent of UMI base pairs with quality scores over Q30',
        )

        self.add_metric(
            name='No PolyT Reads',
            value=self.no_polyT_num,
            total=self.total_num,
            show=False
        )

        self.add_metric(
            name='Low Quality Reads',
            value=self.lowQual_num,
            total=self.total_num,
            show=False,
        )

        self.add_metric(
            name='No Linker Reads',
            value=self.no_linker_num,
            total=self.total_num,
            show=False,
        )

        self.add_metric(
            name='No Barcode Reads',
            value=self.no_barcode_num,
            total=self.total_num,
            show=False,
        )

        self.add_metric(
            name='Corrected Linker Reads',
            value=self.linker_corrected_num,
            total=self.total_num,
            show=False,
        )

        self.add_metric(
            name='Corrected Barcode Reads',
            value=self.barcode_corrected_num,
            total=self.total_num,
            show=False,
        )

        if self.clean_num == 0:
            raise Exception('no valid reads found! please check the --protocol parameter.' + '`--protocol auto` can auto-detect GEXSCOPE-MicroBead,GEXSCOPE-V1,GEXSCOPE-V2,sCircle-V1. `--protocol new` is used for user defined combinations that you need to provide `--pattern`, `--whitelist` and `--linker` at the same time.')
        
        if self.assay == 'flv_trust4':
            self.add_metric(
                name='Valid Matched Reads',
                value=self.match_num,
                total=self.total_num,
                help_info='reads match with flv_rna cell barcodes'
            )

            self.add_metric(
                name='Matched Barcodes',
                value=len(self.match_cbs),
                help_info='barcodes match with flv_rna'
            )


    def run(self):
        """
        Extract barcode and UMI from R1. Filter reads with 
            - invalid polyT
            - low quality in barcode and UMI
            - invalid inlinker
            - invalid barcode
            
        for every sample
            get chemistry
            get linker_mismatch_dict and barcode_mismatch_dict
            for every read in read1
                filter
                write valid R2 read to file
        """

        for i in range(self.fq_number):

            protocol = self.protocol[i]
            lowNum = int(self.lowNum)
            lowQual = int(self.lowQual)

            pattern_dict = parse_pattern(self.pattern)
            whitelist_file = self.whitelist_file
            linker_file = self.linker_file

            bool_T = True if 'T' in pattern_dict else False
            bool_L = True if 'L' in pattern_dict else False
            bool_whitelist = (whitelist_file is not None) and whitelist_file != "None"
            C_len = sum([item.stop - item.start for item in pattern_dict['C']])

            if bool_whitelist:
                barcode_set_list, barcode_mismatch_list = Barcode.parse_whitelist_file(whitelist_file,
                                                n_pattern=len(pattern_dict['C']), n_mismatch=1)
            if bool_L:
                if(self.protocol == "GEXSCOPE-sCircle"):
                    linker_set_list, linker_mismatch_list = Barcode.parse_whitelist_file(linker_file,
                                                   n_pattern=len(pattern_dict['L']), n_mismatch=2)
                else:
                    linker_set_list, linker_mismatch_list = Barcode.parse_whitelist_file(linker_file,
                                                    n_pattern=len(pattern_dict['L']) - 1, n_mismatch=2)  # L16L16L1

            with pysam.FastxFile(self.fq1_list[i], persist=False) as fq1, \
                    pysam.FastxFile(self.fq2_list[i], persist=False) as fq2:
                for entry1, entry2 in zip(fq1, fq2):
                    header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                    self.total_num += 1

                    # polyT filter
                    if bool_T and self.filterNoPolyT:
                        if not Barcode.check_polyT(seq1, pattern_dict):
                            self.no_polyT_num += 1
                            if self.nopolyT:
                                self.fh_nopolyT_fq1.write('@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                                self.fh_nopolyT_fq2.write('@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                            continue

                    # lowQual filter
                    C_U_quals_ascii = get_seq_str(qual1, pattern_dict['C'] + pattern_dict['U'])
                    if lowQual > 0 and Barcode.low_qual(C_U_quals_ascii, lowQual, lowNum):
                        self.lowQual_num += 1
                        continue

                    # linker filter
                    if bool_L and (not self.allowNoLinker):
                        if self.protocol == "GEXSCOPE-sCircle":
                            seq_list = get_seq_list(seq1, pattern_dict, 'L')
                        if self.protocol == "GEXSCOPE-V1" or self.protocol == "GEXSCOPE-V2":
                            seq_list = get_seq_list(seq1, pattern_dict, 'L')[0:2]    # L16L16L1, del L1

                        bool_valid, bool_corrected, _ = Barcode.check_seq_mismatch(seq_list, linker_set_list, linker_mismatch_list)
                        if not bool_valid:
                            self.no_linker_num += 1
                            if self.noLinker:
                                self.fh_nolinker_fq1.write(f'@{header1}\n{seq1}\n{seq_str}\n{qual1}\n')
                                self.fh_nolinker_fq2.write('@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                            continue
                        elif bool_corrected:
                            self.linker_corrected_num += 1

                    # barcode filter
                    seq_list = get_seq_list(seq1, pattern_dict, 'C')
                    if self.bool_flv:
                        seq_list = [utils.reverse_complement(seq) for seq in seq_list[::-1]]
                    if bool_whitelist:
                        bool_valid, bool_corrected, corrected_seq = Barcode.check_seq_mismatch(
                            seq_list, barcode_set_list, barcode_mismatch_list)

                        if not bool_valid:
                            self.no_barcode_num += 1
                            continue
                        elif bool_corrected:
                            self.barcode_corrected_num += 1
                        cb = corrected_seq
                    else:
                        cb = "_".join(seq_list)

                    self.clean_num += 1
                    self.barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
                    self.umi_qual_Counter.update(C_U_quals_ascii[C_len:])

                    umi = get_seq_str(seq1, pattern_dict['U'])
                    if not umi:
                        continue

                    if self.bool_flv:
                        qual1 = 'F' * len(cb + umi)
                        self.barcode_read_Counter.update(cb)
                        if self.assay == 'flv_trust4' and cb in self.match_barcodes:
                            self.match_num += 1
                            self.match_cbs.add(cb)
                            if self.barcode_read_Counter[cb] <= 80000:
                                self.fh_fq2.write(f'@{cb}:{umi}:{self.total_num}\n{seq2}\n+\n{qual2}\n')
                                self.fh_fq1.write(f'@{cb}:{umi}:{self.total_num}\n{cb}{umi}\n+\n{qual1}\n')
                    else:
                        if self.args.stdout:
                            print(f'@{cb}:{umi}:{self.total_num}\n{seq2}\n+\n{qual2}')
                        else:
                            self.fh_fq2.write(f'@{cb}:{umi}:{self.total_num}\n{seq2}\n+\n{qual2}\n')
                        if self.output_R1:
                            self.fh_fq1.write(f'@{cb}:{umi}:{self.total_num}\n{seq1}\n+\n{qual1}\n')                   
            
                logger.info(self.fq1_list[i] + ' finished.')

        self.close_files()
        self.add_step_metrics()


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--fq1', help='R1 fastq file. Multiple files are separated by comma.', required=True)
    parser.add_argument('--fq2', help='R2 fastq file. Multiple files are separated by comma.', required=True)
    #parser.add_argument('--sample', required=True)
    parser.add_argument('--assets_dir', required=True)
    parser.add_argument('--assay', required=True)
    parser.add_argument('--protocol', default='auto')
    parser.add_argument('--pattern')
    parser.add_argument('--whitelist', help='Cell barcode whitelist file path, one cell barcode per line.')
    parser.add_argument('--linker', help='Linker whitelist file path, one linker per line.')
    parser.add_argument('--lowNum',help='The maximum allowed lowQual bases in cell barcode and UMI.',
                        type=int,default=2)
    parser.add_argument('--lowQual',
                        help='Default 0. Bases in cell barcode and UMI whose phred value are lower than lowQual will be regarded as low-quality bases.',
                        type=int,default=0)
    parser.add_argument('--filterNoPolyT',help="Filter reads without PolyT.",action='store_true')
    parser.add_argument('--allowNoLinker',help="Allow valid reads without correct linker.",action='store_true')
    parser.add_argument('--nopolyT',help='Outputs R1 reads without polyT.',action='store_true')
    parser.add_argument('--noLinker',help='Outputs R1 reads without correct linker.',action='store_true')
    parser.add_argument('--output_R1',help="Output valid R1 reads.",action='store_true')
    parser.add_argument('--match_dir', help='Matched scRNA-seq directory, required for flv_trust4')
    parser.add_argument('--stdout',help="Write output to standard output",action='store_true')

    # add version
    parser.add_argument('--version', action='version', version='1.0')
    parser = s_common(parser)  #
    args = parser.parse_args()

    with Barcode(args, display_title='Demultiplexing') as runner:
        runner.run()
    #runner = Barcode(args)
    #runner.run()

import gzip
import os
import pandas as pd
import pysam
from collections import Counter, defaultdict
import glob

import logging
logger = logging.getLogger(__name__)

## init
OUTS_DIR = 'outs'
RAW_MATRIX_DIR_SUFFIX = 'raw'
FILTERED_MATRIX_DIR_SUFFIX = 'filtered'
BARCODE_FILE_NAME = 'barcodes.tsv.gz'

def openfile(file_name, mode='rt', **kwargs):
    """open gzip or plain file"""
    if file_name.endswith('.gz'):
        file_obj = gzip.open(file_name, mode=mode, **kwargs)
    else:
        file_obj = open(file_name, mode=mode, **kwargs)
    return file_obj


def read_one_col(file):
    """
    Read file with one column. Strip each line.
    Returns col_list, line number
    """
    df = pd.read_csv(file, header=None)
    col1 = list(df.iloc[:, 0])
    col1 = [item.strip() for item in col1]
    num = len(col1)
    return col1, num


def read_fasta(fasta_file, equal=False):
    """
    Args:
        equal: if True, seq in fasta must have equal length
    Returns:
        {seq_id: seq} dict
    """
    fa_dict = {}
    length = None
    with pysam.FastxFile(fasta_file) as infile:
        for index, record in enumerate(infile):
            seq = record.sequence
            if index == 0:
                length = len(seq)
            if equal:
                if length != len(seq):
                    raise Exception(f"{fasta_file} have different seq length")
            fa_dict[record.name] = seq
    return fa_dict, length


def reverse_complement(seq):
    """Reverse complementary sequence

    :param original seq
    :return Reverse complementary sequence
    """
    return str(Seq(seq).reverse_complement())


def genDict(dim=3, valType=int):
    if dim == 1:
        return defaultdict(valType)
    else:
        return defaultdict(lambda: genDict(dim - 1, valType=valType))


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


def hamming_correct(string1, string2):
    threshold = len(string1) / 10 + 1
    if hamming_distance(string1, string2) < threshold:
        return True
    return False


def check_arg_not_none(args, arg_name):
    """
    check if args.arg_name is not None
    Args:
        args: argparser args
        arg_name: argparser arg name
    Return:
        bool
    """
    arg_value = getattr(args, arg_name, None)
    if arg_value and arg_value.strip() != 'None':
        return True
    else:
        return False


def glob_file(pattern_list: list):
    """
    glob file among pattern list
    Returns:
        PosixPath object
    Raises:
        FileNotFoundError: if no file found
        MultipleFileFound: if more than one file is found
    """
    if not isinstance(pattern_list, list):
        raise TypeError('pattern_list must be a list')

    match_list = []
    for pattern in pattern_list:
        files = glob.glob(pattern)
        if files:
            for f in files:
                match_list.append(f)
    
    if len(match_list) == 0:
        raise FileNotFoundError(f'No file found for {pattern_list}')
    
    if len(match_list) > 1:
        raise MultipleFileFoundError(
            f'More than one file found for pattern: {pattern_list}\n'
            f'File found: {match_list}'
        )
    
    return match_list[0]


def get_matrix_dir_from_match_dir(match_dir):
    """
    Returns:
        matrix_dir: PosixPath object
    """
    matrix_dir = f"{match_dir}/{OUTS_DIR}/{FILTERED_MATRIX_DIR_SUFFIX}"
    if not os.path.exists(matrix_dir):
        raise FileNotFoundError(f'{matrix_dir} not found')
    
    return matrix_dir


def parse_match_dir(match_dir):
    '''
    return dict
    keys: 'match_barcode', 'n_match_barcode', 'matrix_dir', 'tsne_coord'

    OUTS_DIR = 'outs'
    '''
    match_dict = {}

    pattern_dict = {
        'tsne_coord': [f'{match_dir}/outs/tsne_coord.tsv'],   # {match_dir}/{OUTS_DIR}
        'markers': [f'{match_dir}/outs/markers.tsv'],
        'h5ad': [f'{match_dir}/outs/rna.h5ad'],
    }

    for file_key in pattern_dict:
        file_pattern= pattern_dict[file_key]
        try:
            match_file = glob_file(file_pattern)
        except FileNotFoundError:
            logger.warning(f"No {file_key} found in {match_dir}")
        else:
            match_dict[file_key] = match_file

    match_dict['matrix_dir'] = get_matrix_dir_from_match_dir(match_dir)
    match_barcode, n_match_barcode = get_barcode_from_match_dir(match_dir)
    match_dict['match_barcode'] = match_barcode
    match_dict['n_match_barcode'] = n_match_barcode

    return match_dict


def get_matrix_file_path(matrix_dir, file_name):
    """
    compatible with non-gzip file
    """
    non_gzip = file_name.strip('.gz')
    file_path_list = [f'{matrix_dir}/{file_name}', f'{matrix_dir}/{non_gzip}']
    for file_path in file_path_list:
        if os.path.exists(file_path):
            return file_path


def get_barcode_from_matrix_dir(matrix_dir):
    """
    Returns:
        match_barcode: list
        no_match_barcode: int
    """
  
    match_barcode_file = get_matrix_file_path(matrix_dir, BARCODE_FILE_NAME)
    match_barcode, n_match_barcode = read_one_col(match_barcode_file)

    return match_barcode, n_match_barcode


def get_barcode_from_match_dir(match_dir):
    '''
    multi version compatible
    Returns:
        match_barcode: list
        no_match_barcode: int
    '''
    matrix_dir = get_matrix_dir_from_match_dir(match_dir)
    return get_barcode_from_matrix_dir(matrix_dir)
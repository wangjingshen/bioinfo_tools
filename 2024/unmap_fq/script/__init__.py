import os
from collections import OrderedDict

__VERSION__ = "2.0.7"
__version__ = __VERSION__

ASSAY_LIST = [
    "rna",
    'vdj',
    'tag',
    'dynaseq',
    'snp',
    'capture_virus',
    'fusion',
    'citeseq',
    'flv_trust4',
    'sweetseq',
    'utils',
    'bulk_vdj',
    'bulk_rna',
]

ROOT_PATH = os.path.dirname(__file__)

RELEASED_ASSAYS = ['rna', 'vdj', 'tag', 'dynaseq', 'snp', 'capture_virus', 'fusion', 'citeseq', 'flv_trust4', 'sweetseq','bulk_vdj','bulk_rna']

# argument help
HELP_DICT = {
    'match_dir': 'Match celescope scRNA-Seq directory.',
    'gene_list': 'Required. Gene list file, one gene symbol per line. Only results of these genes are reported. Conflict with `--panel`',
    'genomeDir': 'Genome directory after running `celescope {assay} mkref`.',
    'thread': 'Thread to use.',
    'debug': 'If this argument is used, celescope may output addtional file for debugging.',
    'fasta': 'Required. Genome fasta file. Use absolute path or relative path to `genomeDir`.',
    'outdir': 'Output directory.',
    'matrix_dir': 'Match celescope scRNA-Seq matrix directory.',
    'panel': 'The prefix of bed file in `celescope/data/snp/panel/`, such as `lung_1`. Conflict with `--gene_list`',
    'virus_genomeDir': 'Required. Virus genome directory after running `celescope capture_virus mkref`.',
    'threshold_method': 'One of [otsu, auto, hard, none].',
    'tsne_file': 'match_dir t-SNE coord file. Do not required when `--match_dir` is provided.',
    'df_marker_file': 'match_dir df_marker_file. Not required when `--match_dir` is provided.',
    'cell_calling_method': 'Default `EmptyDrops_CR`. Choose from [`auto`, `EmptyDrops_CR`]',
    'additional_param': 'Additional parameters for the called software. Need to be enclosed in quotation marks. For example, `--{software}_param "--param1 value1 --param2 value2"`.',
    'genomeSAindexNbases': '''For small genomes, the parameter --genomeSAindexNbases must to be scaled down, with a typical 
value of min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal 
to 9, for 100 kiloBase genome, this is equal to 7.''',
    'chemistry': '`--chemistry auto` can auto-detect scopeV2 mRNA, scopeV3 mRNA, full length VDJ mRNA(flv_rna) and full length VDJ(flv). You need to explicitly use `--chemistry scopeV1` for legacy chemistry scopeV1. `--chemistry customized` is used for user defined combinations that you need to provide `--pattern`, `--whitelist` and `--linker` at the same time.',   
}

# report metrics help
HELP_INFO_DICT = {
    'matched_barcode_number': {
        'display': 'Number of Matched Cells',
        'info': 'cell barcode number of matched scRNA-Seq sample',
    }
}


# from celescope tools __init__.py
# barcode
PATTERN_DICT = {
    'auto': None,
    'scopeV1': 'C12U8T18',
    'scopeV2.0.0': 'C8L16C8L16C8U8T18',
    'scopeV2.0.1': 'C8L16C8L16C8L1U8T18',
    'scopeV2.1.0': 'C8L16C8L16C8U12T18',
    'scopeV2.1.1': 'C8L16C8L16C8L1U12T18',
    'scopeV2.2.1': 'C8L16C8L16C8L1U12T18',
    'scopeV3.0.1': 'C9L16C9L16C9L1U12T18',
    # flv_rna is actually U9L16, but the last 10 bases can not be sequenced accurately.
    'flv_rna': 'C8L16C8L16C8U9L6',
    'flv' : 'U9C8L16C8L16C8',
    'bulk_vdj': 'L18C6U16',
    'bulk_rna': 'C9U12',
    'customized': None,
}


# count
OUTS_DIR = 'outs'
RAW_MATRIX_DIR_SUFFIX = 'raw'
FILTERED_MATRIX_DIR_SUFFIX = 'filtered'
MATRIX_FILE_NAME = 'matrix.mtx.gz'
FEATURE_FILE_NAME = 'features.tsv.gz'
BARCODE_FILE_NAME = 'barcodes.tsv.gz'
STAR_BAM_SUFFIX = 'Aligned.out.bam'
TAG_BAM_SUFFIX = 'aligned_posSorted_addTag.bam'
STARSOLO_BAM_SUFFIX = 'Aligned.sortedByCoord.out.bam'
COUNTS_FILE_NAME = 'counts.tsv'

# mkref
GENOME_CONFIG = 'celescope_genome.config'

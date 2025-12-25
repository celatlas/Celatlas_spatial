import os

__VERSION__ = "1.7.0"
__version__ = __VERSION__

ASSAY_LIST = [
    "rna",
]
RELEASED_ASSAYS = ["rna"]

# environment variables
genome = {
    "human": "./reference/Homo_sapiens/",
    "mouse": "./reference/Mus_musculus/",
    "rat": "./reference/Rattus_norvegicus/",
}
ROOT_PATH = os.path.dirname(__file__)

# mapping
# BC_WIDTH: Barcode length in base pairs (chemistry-dependent)
# - BBV2: 12bp (C4+C4+C4) → BIT_WIDTH=24
# - BBV3: 15bp (C5+C5+C5) → BIT_WIDTH=30
# - BBV0: 18bp (C6+C6+C6) → BIT_WIDTH=36
# Default to maximum for backward compatibility
BC_WIDTH = 18  # Maximum supported barcode length
BIT_WIDTH = 2 * BC_WIDTH  # 36 bits (2 bits per base)
BIT_WIDTH_2 = 4 * BC_WIDTH  # 72 bits (4 bits per base for alternative encoding)

# Chemistry-specific barcode widths
CHEMISTRY_BC_WIDTH = {
    'BBV0': 18,      # C6+C6+C6
    'BBV1': 12,      # C4+C4+C4
    'BBV2': 12,      # C4+C4+C4
    'BBV2.4': 12,    # C4+C4+C4
    'BBV3': 15,      # C5+C5+C5
    'BBV3.1': 15,    # C5+C5+C5 (same as BBV3)
    'strnaV2.2': 12, # C4+C4+C4
}

DECODE_MAPPING = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
ENCODE_MAPPING = {v: k for k, v in DECODE_MAPPING.items()}

DECODE_MAPPING_2 = {'0000': 'N', '0001': 'A', '0010': 'G', '0100': 'C', '1000': 'T'}
ENCODE_MAPPING_2 = {v: k for k, v in DECODE_MAPPING_2.items()}

# argument help
HELP_DICT = {
    'match_dir': 'Match celatlas scRNA-Seq directory.',
    'gene_list': 'Required. Gene list file, one gene symbol per line. Only results of these genes are reported. '
                 'Conflict with `--panel`',
    'genomeDir': 'Required. Genome directory after running `celatlas {assay} mkref`.',
    'thread': 'Thread to use.',
    'debug': 'If this argument is used, celatlas may output addtional file for debugging.',
    'fasta': 'Required. Genome fasta file. Use absolute path or relative path to `genomeDir`.',
    'outdir': 'Output directory.',
    'matrix_dir': 'Match celatlas scRNA-Seq matrix directory.',
    'panel': 'The prefix of bed file in `celatlas/data/snp/panel/`, such as `lung_1`. Conflict with `--gene_list`',
    'virus_genomeDir': 'Required. Virus genome directory after running `celatlas capture_virus mkref`.',
    'threshold_method': 'One of [otsu, auto, hard, none].',
    'tsne_file': 'match_dir t-SNE coord file. Do not required when `--match_dir` is provided.',
    'df_marker_file': 'match_dir df_marker_file. Not required when `--match_dir` is provided.',
    'cell_calling_method': 'Default `EmptyDrops_CR`. Choose from [`auto`, `EmptyDrops_CR`]',
    'additional_param': 'Additional parameters for the called software. Need to be enclosed in quotation marks. '
                        'For example, `--{software}_param "--param1 value1 --param2 value2"`.',
    'genomeSAindexNbases': '''For small genomes, the parameter --genomeSAindexNbases must to be scaled down, 
                        with a typical value of min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, 
                        this is equal to 9, for 100 kiloBase genome, this is equal to 7.''',
    'chemistry': '`--chemistry auto` can auto-detect chemistry from R1 reads. '
                 'Chemistries: BBV0 (C6-space, scRNA only), BBV1 (C4-fixed, universal), BBV2 (C4-space, universal, production), BBV3 (C5-space, universal, future). '
                 'Legacy chemistries (strnaV1.0, strnaV2.x, BBV2.2, BBV2.4, BBV3.1) are deprecated but still supported. '
                 '`--chemistry customized` is used for user defined combinations that you need to provide `--pattern`, `--whitelist` and `--linker` at the same time. '
                 'Use `--mode scrna` for single-cell RNA-seq or `--mode strna` for spatial transcriptomics. '
                 'Note: BBV2/BBV3 support both scrna and strna modes.',
}

import os
import re
import cv2
import sys
import glob
import gzip
import time
import uuid
import dnaio
import pysam
import difflib
import logging
import functools
import warnings
import resource
import tifffile
import importlib
import subprocess
import numpy as np
import pandas as pd

from tqdm import tqdm
from Bio.Seq import Seq
from collections import Counter, defaultdict
from xopen import xopen
from functools import wraps
from datetime import timedelta
from multiprocessing import Pool
from celatlas_spatial.tools.avltree import AVLTree, OptimizedAVLTree
from celatlas_spatial.tools.bktree import BKTree
from celatlas_spatial.tools.trietree import Trie
from celatlas_spatial.tools.__init__ import BARCODE_FILE_NAME, FILTERED_MATRIX_DIR_SUFFIX
from celatlas_spatial.__init__ import ROOT_PATH, DECODE_MAPPING, DECODE_MAPPING_2, ENCODE_MAPPING_2, ENCODE_MAPPING, \
    BIT_WIDTH, BC_WIDTH


def add_log(func):
    """
    logging start and done.
    """
    logFormatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    module = func.__module__
    name = func.__name__
    logger_name = f'{module}.{name}'
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler(sys.stdout)
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    @wraps(func)
    def wrapper(*args, **kwargs):

        logger.info('start...')
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info('done. time used: %s', used)
        return result

    wrapper.logger = logger
    return wrapper

def add_mem(func):
    """
    logging mem.
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    module = func.__module__
    name = func.__name__
    logger_name = f'{module}.{name}'
    logger = logging.getLogger(logger_name)

    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.info(using("before"))
        result = func(*args, **kwargs)
        logger.info(using("after"))
        return result

    wrapper.logger = logger
    return wrapper

def deprecated(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        warnings.warn(f"Call to deprecated function {func.__name__}.", category=DeprecationWarning, stacklevel=2)
        return func(*args, **kwargs)
    return wrapper


def dict_gmt_txt(output_file, dict_data):
    """Write dictionary to gmt txt format

    Parameters:
    -----------
    dict_data : dict
        Dictionary where keys are pathways and values are gene lists
    output_file : str
        Path to output file
    """
    with open(output_file, 'w') as f:
        for pathway, genes in dict_data.items():
            if isinstance(genes, list) and genes:
                genes_str = '\t'.join(str(gene) for gene in genes)
                f.write(f"{pathway}\t\t{genes_str}\t\n")

def using(point=""):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return '''%s: usertime=%s systime=%s mem=%s mb
        ''' % (point, usage[0], usage[1],
               usage[2]/1024.0)

def check_mkdir(dir_name):
    """if dir_name is not exist, make one"""
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")

def find_assay_init(assay):
    init_module = importlib.import_module(f"celatlas_spatial.{assay}.__init__")
    return init_module

def find_step_module(assay, step):
    file_path_dict = {
        'assay': f'{ROOT_PATH}/{assay}/{step}.py',
        'tools': f'{ROOT_PATH}/tools/spatial_{step}.py',
    }

    init_module = find_assay_init(assay)
    if os.path.exists(file_path_dict['assay']):
        step_module = importlib.import_module(f"celatlas_spatial.{assay}.{step}")
    elif hasattr(init_module, 'IMPORT_DICT') and step in init_module.IMPORT_DICT:
        module_path = init_module.IMPORT_DICT[step]
        step_module = importlib.import_module(f"{module_path}.{step}")
    elif os.path.exists(file_path_dict['tools']):
        step_module = importlib.import_module(f"celatlas_spatial.tools.spatial_{step}")
    else:
        raise ModuleNotFoundError(f"No module found for {assay}.{step}")

    return step_module

def sort_bam(input_bam, output_bam, threads=1, by='coord'):
    cmd = (
        f'samtools sort {input_bam} '
        f'-o {output_bam} '
        f'--threads {threads} '
        '2>&1 '
    )
    if by == "name":
        cmd += " -n"
    subprocess.check_call(cmd, shell=True)

def index_bam(input_bam):
    cmd = f"samtools index {input_bam} 2>&1 "
    subprocess.check_call(cmd, shell=True)

def format_number(number: int) -> str:
    return format(number, ",")

def genDict(dim=3, valType=int):
    if dim == 1:
        return defaultdict(valType)
    else:
        return defaultdict(lambda: genDict(dim - 1, valType=valType))

def generic_open(file_name, *args, **kwargs):
    if file_name.endswith('.gz'):
        file_obj = gzip.open(file_name, *args, **kwargs)
    else:
        file_obj = open(file_name, *args, **kwargs)
    return file_obj

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

def int_to_uuid(int_id):
    return uuid.UUID(int=int_id, version=4)

def cat(input, output):

    import glob
    import os
    from xopen import xopen
    import tempfile
    import shutil
    
    if isinstance(input, str):
        if '*' in input:
            input_files = sorted(glob.glob(input))
        else:
            input_files = [input]
    else:
        input_files = input
    
    if not input_files:
        raise FileNotFoundError(f"No files found matching: {input}")
    
    valid_files = []
    for input_file in input_files:
        if os.path.exists(input_file) and os.path.getsize(input_file) > 0:
            valid_files.append(input_file)
        elif os.path.exists(input_file):
            print(f"Warning: Empty input file skipped: {input_file}")
        else:
            print(f"Warning: Input file does not exist: {input_file}")
    
    if not valid_files:
        print(f"Warning: No valid files to merge, creating empty output: {output}")

        with xopen(output, 'wt') as f:
            pass
        return
    
    temp_output = None
    try:

        output_dir = os.path.dirname(output)
        output_basename = os.path.basename(output)
        if output_basename.endswith('.gz'):
            temp_suffix = '.tmp.gz'  
        else:
            temp_suffix = '.tmp'
        temp_fd, temp_output = tempfile.mkstemp(suffix=temp_suffix, dir=output_dir)
        os.close(temp_fd)  
        
        print(f"Merging {len(valid_files)} files to temporary file: {temp_output}")
        
        with xopen(temp_output, 'wt') as outf:  
            for i, input_file in enumerate(valid_files):
                print(f"Processing file {i+1}/{len(valid_files)}: {input_file}")
                try:
                    with xopen(input_file, 'rt') as inf: 

                        line_count = 0
                        for line in inf:
                            outf.write(line)
                            line_count += 1
                        print(f"  Processed {line_count} lines from {input_file}")
                except Exception as e:
                    print(f"Error processing input file {input_file}: {e}")
                    raise
        
        print(f"Moving temporary file to final output: {output}")
        shutil.move(temp_output, output)
        temp_output = None  
        
        print(f"Successfully merged {len(valid_files)} files to {output}")
        
    except Exception as e:
        print(f"Error during file merge: {e}")
        import traceback
        traceback.print_exc()

        if temp_output and os.path.exists(temp_output):
            try:
                os.remove(temp_output)
                print(f"Cleaned up temporary file: {temp_output}")
            except:
                pass
        raise
    finally:

        if temp_output and os.path.exists(temp_output):
            try:
                os.remove(temp_output)
            except:
                pass

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

def compare_strings(seq, seq_list):
    result = []
    for s in seq_list:
        diff = list(difflib.ndiff(seq, s))
        diff_chars = [c[-1] for c in diff if c.startswith('- ')]
        result.append(''.join(diff_chars))
    return result

def decoder(binary_str):
    return ''.join(DECODE_MAPPING.get(binary_str[i:i+2], 'N') for i in range(0, len(binary_str), 2))

def decoder_2(binary_str):
    return ''.join(DECODE_MAPPING_2.get(binary_str[i:i+4], 'N') for i in range(0, len(binary_str), 4))

def encoder(bc, bc_width=None):
    """
    Encode barcode string to integer

    Args:
        bc: barcode string
        bc_width: expected barcode width (optional, for validation)

    Returns:
        np.uint64: encoded integer
    """
    if bc_width is not None and len(bc) != bc_width:
        raise ValueError(f"Barcode length must be {bc_width}bp. Got {len(bc)}bp.")

    # Dynamic encoding without length constraint if bc_width not specified
    binary_str = ''.join(ENCODE_MAPPING.get(char, '11') for char in bc)
    decimal = np.uint64(int(binary_str, 2))
    return decimal

def encoder_2(bc, bc_width=None):
    """
    Encode barcode string to integer (alternative encoding)

    Args:
        bc: barcode string
        bc_width: expected barcode width (optional, for validation)

    Returns:
        np.uint64: encoded integer
    """
    if bc_width is not None and len(bc) != bc_width:
        raise ValueError(f"Barcode length must be {bc_width}bp. Got {len(bc)}bp.")

    # Dynamic encoding without length constraint if bc_width not specified
    binary_str = ''.join(ENCODE_MAPPING_2.get(char, '0000') for char in bc)
    decimal = np.uint64(int(binary_str, 2))
    return decimal

def matrix_decoder(mtx, bit_width=None):
    """
    Decode H5 matrix to barcodes

    Args:
        mtx: H5 matrix containing encoded barcodes
        bit_width: bit width for decoding (optional, uses BIT_WIDTH if not specified)

    Returns:
        bcs: list of barcode strings
        bcs_tree: Trie tree for fast lookup
        bcs_df: DataFrame with barcode and coordinates
        n_count: number of barcodes containing 'N'
    """
    if bit_width is None:
        from celatlas_spatial import BIT_WIDTH as default_bit_width
        bit_width = default_bit_width

    mtx = mtx[:].squeeze()
    non_zero_indices = np.nonzero(mtx)
    bcs_deci = mtx[non_zero_indices].astype(np.int64)
    binary_repr_func = np.vectorize(np.binary_repr)
    bcs_bin = binary_repr_func(bcs_deci, width=bit_width)
    bcs_series = pd.Series(bcs_bin)
    bcs = bcs_series.apply(decoder)

    n_count = bcs.str.contains('N').sum()
    bcs = bcs.tolist()
    indices = np.array(list(zip(*non_zero_indices)))[:, ::-1]
    bcs_df = pd.DataFrame({'barcode': bcs, 'coord': list(map(tuple, indices))})

    # bk tree build
    # bcs_tree = BKTree(items=bcs_deci)
    # bcs_tree = OptimizedAVLTree(items=bcs_deci)
    bcs_tree = Trie(items=bcs_deci)
    # bc dict build
    # bcs_bpmap = {x: x for x in bcs_deci}
    return bcs, bcs_tree, bcs_df, n_count

def calc_hamming_distance(string1, string2):
    distance = 0
    length = len(string1)
    length2 = len(string2)
    if length != length2:
        raise Exception(f"string1({length}) and string2({length2}) do not have same length")
    for i in range(length):
        if string1[i] != string2[i]:
            distance += 1
    return distance

def read_fastq_files(fastq_file1, fastq_file2, entry_queue, worker_id, num_workers):
    """
    Read fastq files and put entry to entry_queue
    """
    n_idx = 0
    with dnaio.open(fastq_file1, fastq_file2) as f:
        for line in f:
            if n_idx % num_workers == worker_id:
                fq1 = line[0]
                fq2 = line[1]
                entry_data = (
                    (fq1.name, fq1.sequence, fq1.qualities),
                    (fq2.name, fq2.sequence, fq2.qualities)
                )
                entry_queue.put(entry_data)
            n_idx += 1

def reverse_complement(seq):
    """
    Reverse complementary sequence

    :param original seq
    :return Reverse complementary sequence
    """
    return str(Seq(seq).reverse_complement())

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

def get_chunk_paired_zip(fastq_files, file_path, thread, resume=False):
    """
    Get chunk paired zip (resume is supported)
    Additionally, to prevent the program from terminating due to unexpected circumstances, which may result
    in all chunks of fastq files being constructed but subsequent processes failing to execute properly and
    causing the program to exit, re-executing this module would consume a considerable amount of time.
    Therefore, we have introduced support for a resume operation. After verifying that all chunks have been
    properly divided, the program can proceed directly to execute subsequent processing, thus avoiding redundant
    operations.
    """
    def get_total_reads(fastq_file):
        cmd = f"zcat {fastq_file} | wc -l"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return int(result.stdout.strip()) // 4

    chunk_id = 0
    if resume:
        print('resume mode')
        for idx in os.listdir(file_path):
            chunk_path = os.path.join(file_path, idx)
            if os.path.isdir(chunk_path):
                chunk_id += 1
        if chunk_id != 0:
            return chunk_id

    # total_reads = get_total_reads(fastq_files[0])
    total_reads = os.path.getsize(fastq_files[0]) // 55
    chunk_reads = np.ceil(total_reads / thread)
    print(f'chunk size {chunk_reads:,} reads')

    with dnaio.open(fastq_files[0], fastq_files[1]) as f:
        chunk_id, chunk_num = 0, 0
        chunk_path = os.path.join(file_path, f'{chunk_id}')
        os.makedirs(chunk_path, exist_ok=True)
        chunk1_path = os.path.join(chunk_path, f'{chunk_id}_1.fq.gz')
        chunk2_path = os.path.join(chunk_path, f'{chunk_id}_2.fq.gz')

        for fq1, fq2 in tqdm(f):
            header1, seq1, qual1 = fq1.name, fq1.sequence, fq1.qualities
            header2, seq2, qual2 = fq2.name, fq2.sequence, fq2.qualities
            if chunk_num == 0:
                f1 = xopen(chunk1_path, 'w')
                f2 = xopen(chunk2_path, 'w')

            f1.write(f'@{header1}\n{seq1}\n+\n{qual1}\n')
            f2.write(f'@{header2}\n{seq2}\n+\n{qual2}\n')
            chunk_num += 1

            if chunk_num == chunk_reads and chunk_id < thread - 1:
                f1.close()
                f2.close()

                chunk_id += 1
                chunk_num = 0
                chunk_path = os.path.join(file_path, f'{chunk_id}')
                os.makedirs(chunk_path, exist_ok=True)
                chunk1_path = os.path.join(chunk_path, f'{chunk_id}_1.fq.gz')
                chunk2_path = os.path.join(chunk_path, f'{chunk_id}_2.fq.gz')

        f1.close()
        f2.close()
    thread = chunk_id + 1
    return thread

def get_matrix_file_path(matrix_dir, file_name):
    """
    compatible with gzip file
    """
    file_path_list = [f'{matrix_dir}/{file_name}', f'{matrix_dir}/{file_name}.gz']
    for file_path in file_path_list:
        if os.path.exists(file_path):
            return file_path

@add_log
def get_barcode_from_matrix_dir(matrix_dir):
    """
    Returns:
        match_barcode: list
        no_match_barcode: int
    """

    match_barcode_file = get_matrix_file_path(matrix_dir, BARCODE_FILE_NAME)
    match_barcode, n_match_barcode = read_one_col(match_barcode_file)

    return match_barcode, n_match_barcode

@add_log
def get_matrix_dir_from_match_dir(match_dir):
    """
    Returns:
        matrix_dir: PosixPath object
    """
    matrix_dir_pattern_list = []
    for matrix_dir_suffix in FILTERED_MATRIX_DIR_SUFFIX:
        matrix_dir_pattern_list.append(f"{match_dir}/*count/*{matrix_dir_suffix}")

    matrix_dir = glob_file(matrix_dir_pattern_list)
    get_matrix_dir_from_match_dir.logger.info(f"Matrix_dir :{matrix_dir}")

    return matrix_dir

@add_log
def get_barcode_from_match_dir(match_dir):
    """
    multi version compatible
    Returns:
        match_barcode: list
        no_match_barcode: int
    """
    matrix_dir = get_matrix_dir_from_match_dir(match_dir)
    return get_barcode_from_matrix_dir(matrix_dir)

def read_CID(CID_file):
    """
    return df_index, df_valid
    """
    df_index = pd.read_csv(CID_file, sep='\t', index_col=0).reset_index()
    df_valid = df_index[df_index['valid'] == True]
    return df_index, df_valid

def get_assay_text(assay):
    """
    Deprecated
    add spatial transcriptomics prefix
    deprecated
    """
    return 'Spatial-transcriptomics' + assay

@add_log
def parse_match_dir(match_dir):
    """
    return dict
    keys: 'match_barcode', 'n_match_barcode', 'matrix_dir', 'tsne_coord'
    """
    match_dict = {}

    pattern_dict = {
        'tsne_coord': [f'{match_dir}/*analysis*/*tsne_coord.tsv'],
        'markers': [f'{match_dir}/*analysis*/*markers.tsv'],
        'h5ad': [f'{match_dir}/*analysis*/*.h5ad'],
    }

    for file_key in pattern_dict:
        file_pattern= pattern_dict[file_key]
        try:
            match_file = glob_file(file_pattern)
        except FileNotFoundError:
            parse_match_dir.logger.warning(f"No {file_key} found in {match_dir}")
        else:
            match_dict[file_key] = match_file

    match_dict['matrix_dir'] = get_matrix_dir_from_match_dir(match_dir)
    match_barcode, n_match_barcode = get_barcode_from_match_dir(match_dir)
    match_dict['match_barcode'] = match_barcode
    match_dict['n_match_barcode'] = n_match_barcode

    return match_dict

def possible_mislen(bl, mismatch):
    """
    Calculate possible mismatch barcode length
    """
    if mismatch == 1:
        return bl * 3
    elif mismatch == 2:
        return bl * 3 + bl * (bl - 1) * 3 * 3 // 2
    elif mismatch == 3:
        return bl * 3 + bl * (bl - 1) * 3 * 3 // 2 + bl * (bl - 1) * (bl - 2) * 3 * 3 * 3 // 6
    else:
        return 0

def find_directory_by_name(directory, target_name):
    """
    search directory by name, return the first match.If no match, return None
    input:
        directory: str, directory to search
        target_name: str, target name to search
    return:
        match_dict: dict(including FilterBarcodes and TissueBbox)
    """
    match_dict = {}

    for root, dirs, files in os.walk(directory):
        target_pattern = os.path.join(root, f'*{target_name}*')
        matches = glob.glob(target_pattern)

        for match in matches:
            if os.path.isfile(match):
                if 'FilterBarcodes' in match:
                    match_dict['FilterBarcodes'] = match
                elif 'tissue_bbox' in match:
                    match_dict['TissueBbox'] = match

    return match_dict if match_dict else None

def hole_fill(binary_image):
    """hole fill."""
    hole = binary_image.copy()  # copy the binary image
    hole = cv2.copyMakeBorder(hole, 1, 1, 1, 1, cv2.BORDER_CONSTANT, value=[0])  # add a border
    hole2 = hole.copy()
    cv2.floodFill(hole, None, (0, 0), 255)  # floodFill function to fill the hole
    hole = cv2.bitwise_not(hole)
    binary_hole = cv2.bitwise_or(hole2, hole)[1:-1, 1:-1]
    return binary_hole

def multi_process(n, method, args):
    p = Pool(n)
    results = p.map(method, args)
    p.close()
    p.join()
    return results

def read_tiff_with_metadata(tiff_path):
    with tifffile.TiffFile(tiff_path) as tif:
        # read image data
        image_data = tif.asarray()
        # read metadata
        metadata = {}

        # read tags
        for tag in tif.pages[0].tags.values():
            metadata[tag.name] = tag.value

        # read ImageJ metadata (if exists)
        if hasattr(tif, 'imagej_metadata'):
            metadata['imagej_metadata'] = tif.imagej_metadata
        # read OME metadata (if exists)
        if hasattr(tif, 'ome_metadata'):
            metadata['ome_metadata'] = tif.ome_metadata

    type = image_data.dtype.name
    if type != 'uint8':
        image_data = convert_8bit(image_data)

    return image_data, metadata

def convert_8bit(image):
    image = image - image.min()
    image = image / image.max() * 255
    return np.uint8(image)


class MultipleFileFoundError(Exception):
    pass


class Gtf_dict(dict):
    """
    key: gene_id
    value: gene_name
    If the key does not exist, return key. This is to avoid the error:
        The gtf file contains one exon lines with a gene_id, but do not contain a gene line with the same gene_id. FeatureCounts
        work correctly under this condition, but the gene_id will not appear in the Gtf_dict.
    """

    def __init__(self, gtf_file):
        super().__init__()
        self.gtf_file = gtf_file
        self.load_gtf()


    @add_log
    def load_gtf(self):
        """
        get gene_id:gene_name from gtf file
            - one gene_name with multiple gene_id: "_{count}" will be added to gene_name.
            - one gene_id with multiple gene_name: error.
            - duplicated (gene_name, gene_id): ignore duplicated records and print a warning.
            - no gene_name: gene_id will be used as gene_name.

        Returns:
            {gene_id: gene_name} dict
        """
        gene_id_pattern = re.compile(r'gene_id "(\S+)";')
        gene_name_pattern = re.compile(r'gene_name "(\S+)"')
        id_name = {}
        c = Counter()
        with generic_open(self.gtf_file, mode='rt') as f:
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('#'):
                    continue
                tabs = line.split('\t')
                gtf_type, attributes = tabs[2], tabs[-1]
                if gtf_type == 'gene':
                    try:
                        gene_id = gene_id_pattern.findall(attributes)[-1]
                    except IndexError:
                        print(line)
                    gene_names = gene_name_pattern.findall(attributes)
                    if not gene_names:
                        gene_name = gene_id
                    else:
                        gene_name = gene_names[-1]
                    c[gene_name] += 1
                    if c[gene_name] > 1:
                        if gene_id in id_name:
                            assert id_name[gene_id] == gene_name, (
                                'one gene_id with multiple gene_name '
                                f'gene_id: {gene_id}, '
                                f'gene_name this line: {gene_name}'
                                f'gene_name previous line: {id_name[gene_id]}'
                            )
                            self.load_gtf.logger.warning(
                                'duplicated (gene_id, gene_name)'
                                f'gene_id: {gene_id}, '
                                f'gene_name {gene_name}'
                            )
                            c[gene_name] -= 1
                        else:
                            gene_name = f'{gene_name}_{c[gene_name]}'
                    id_name[gene_id] = gene_name
        self.update(id_name)

    def __getitem__(self, key):
        """if key not exist, return key"""
        return dict.get(self, key, key)


class Samtools:
    def __init__(self, in_bam, out_bam, threads=1, debug=False):
        self.in_bam = in_bam
        self.out_bam = out_bam
        self.threads = threads
        self.temp_sam_file = f"{self.out_bam}_sam.temp"
        self.debug = debug

    @add_log
    def samtools_sort(self, in_file, out_file, by='coord'):
        cmd = f"samtools sort {in_file} -o {out_file} --threads {self.threads}"
        if by == "name":
            cmd += " -n"
        self.samtools_sort.logger.debug(cmd)
        subprocess.check_call(cmd, shell=True)

    @add_log
    def samtools_index(self, in_file):
        cmd = f"samtools index {in_file}"
        self.samtools_index.logger.debug(cmd)
        subprocess.check_call(cmd, shell=True)

    def sort_bam(self, by='coord'):
        """sort in_bam"""
        self.samtools_sort(self.in_bam, self.out_bam, by=by)

    def index_bam(self):
        """index out_bam"""
        self.samtools_index(self.out_bam)

    @add_log
    def add_tag(self, gtf_file):
        """
        - CB cell barcode
        - UB UMI
        - GN gene name
        - GX gene id
        """
        gtf_dict = Gtf_dict(gtf_file)

        with pysam.AlignmentFile(self.in_bam, "rb") as original_bam:
            header = original_bam.header
            with pysam.AlignmentFile(self.temp_sam_file, "w", header=header) as temp_sam:
                for read in original_bam:
                    attr = read.query_name.split('_')
                    barcode = attr[0]
                    umi = attr[1]
                    read.set_tag(tag='CB', value=barcode, value_type='Z')
                    read.set_tag(tag='UB', value=umi, value_type='Z')
                    # assign to some gene
                    if read.has_tag('XT'):
                        gene_id = read.get_tag('XT')
                        # if multi-mapping reads are included in original bam,
                        # there are multiple gene_ids
                        if ',' in gene_id:
                            gene_name = [gtf_dict[i] for i in gene_id.split(',')]
                            gene_name = ','.join(gene_name)
                        else:
                            gene_name = gtf_dict[gene_id]
                        read.set_tag(tag='GN', value=gene_name, value_type='Z')
                        read.set_tag(tag='GX', value=gene_id, value_type='Z')
                    temp_sam.write(read)

    @add_log
    def add_RG(self, barcodes):
        """
        barcodes list
        """

        with pysam.AlignmentFile(self.in_bam, "rb") as original_bam:
            header = original_bam.header.to_dict()
            header['RG'] = []
            for index, barcode in enumerate(barcodes):
                header['RG'].append({
                    'ID': barcode,
                    'SM': index + 1,
                })

            with pysam.AlignmentFile(self.temp_sam_file, "w", header=header) as temp_sam:
                for read in original_bam:
                    read.set_tag(tag='RG', value=read.get_tag('CB'), value_type='Z')
                    temp_sam.write(read)

    def temp_sam2bam(self, by=None):
        self.samtools_sort(self.temp_sam_file, self.out_bam, by=by)
        self.rm_temp_sam()

    def rm_temp_sam(self):
        cmd = f"rm {self.temp_sam_file}"
        subprocess.check_call(cmd, shell=True)

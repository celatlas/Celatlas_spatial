import os
import re
import sys
import glob
import json
import pickle
import hashlib
import time

import h5py
import dnaio
import pandas as pd
import pysam
import shutil
import argparse
import numpy as np

from tqdm import tqdm
from xopen import xopen
from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter, defaultdict
from itertools import combinations, product
from celatlas_spatial.celatlas import ArgFormatter
from celatlas_spatial.tools import utils
from celatlas_spatial.tools.step import Step, s_common
from celatlas_spatial.tools.__init__ import PATTERN_DICT
from celatlas_spatial.__init__ import HELP_DICT, ROOT_PATH, BIT_WIDTH, BIT_WIDTH_2
from celatlas_spatial.tools.fastq_utils import FastqFileParser, FastqFileInfo

MIN_T = 10
IS_PARALLEL = True

class ProgressTracker:
    """Progress tracker for resumable processing"""
    
    def __init__(self, outdir, sample):
        self.outdir = outdir
        self.sample = sample
        self.progress_file = os.path.join(outdir, f'{sample}_progress.json')
        self.stats_file = os.path.join(outdir, f'{sample}_stats.pkl')
        self.completed_files = set()
        self.file_stats = {}
        
    def get_file_checksum(self, filepath):
        """Get file checksum based on size and modification time"""
        if not os.path.exists(filepath):
            return None
        stat = os.stat(filepath)
        return f"{stat.st_size}_{int(stat.st_mtime)}"
    
    def save_progress(self, completed_files, file_stats):
        """Save processing progress"""
        progress_data = {
            'completed_files': list(completed_files),
            'file_checksums': {},
            'timestamp': time.time()
        }
        
        # Save checksum for each file
        all_files = set()
        for fq1_list, fq2_list in completed_files:
            all_files.update(fq1_list.split(','))
            all_files.update(fq2_list.split(','))
        
        for filepath in all_files:
            progress_data['file_checksums'][filepath] = self.get_file_checksum(filepath)
        
        with open(self.progress_file, 'w') as f:
            json.dump(progress_data, f, indent=2)
        
        # Save statistics
        with open(self.stats_file, 'wb') as f:
            pickle.dump(file_stats, f)
    
    def load_progress(self):
        """Load processing progress"""
        if not os.path.exists(self.progress_file) or not os.path.exists(self.stats_file):
            return False
        
        try:
            with open(self.progress_file, 'r') as f:
                progress_data = json.load(f)
            
            # Verify if files have changed
            for filepath, old_checksum in progress_data['file_checksums'].items():
                current_checksum = self.get_file_checksum(filepath)
                if current_checksum != old_checksum:
                    print(f"File {filepath} has changed, cannot resume")
                    return False
            
            self.completed_files = set(tuple(pair) for pair in progress_data['completed_files'])
            return True
            
        except Exception as e:
            print(f"Failed to load progress: {e}")
            return False
    
    def load_stats(self):
        """Load statistics"""
        if not os.path.exists(self.stats_file):
            return {}
        
        try:
            with open(self.stats_file, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            print(f"Failed to load stats: {e}")
            return {}
    
    def cleanup(self):
        """Clean up progress files"""
        for filepath in [self.progress_file, self.stats_file]:
            if os.path.exists(filepath):
                os.remove(filepath)

class Chemistry:
    """
    Auto-detect chemistry from R1-read

    Chemistry Versions:
    - BBV0: C6-space (scRNA only, testing version)
    - BBV1: C4-fixed (Universal: scrna + strna, development version)
    - BBV2: C4-space (Universal: scrna + strna, production version)
    - BBV3: C5-space (Universal: scrna + strna, future extended version)
    - Legacy: strnaV1.0, strnaV2.x (deprecated)
    """

    def __init__(self, fq1):
        """
        Initialize chemistry detection with pattern dictionaries

        Chemistry Patterns:
        - 'BBV0': C6-space-L15-C6-L15-C6-U10 (scRNA only, with displacement)
        - 'BBV1': C4-L15-C4-L15-C4-U10 (Universal, fixed linker)
        - 'BBV2': C4-space-L15-C4-L15-C4-U10 (Universal, with displacement, production)
        - 'BBV3': C5-space-L15-C5-L15-C5-U10 (Universal, with displacement, future)

        Deprecated (kept for compatibility):
        - 'strnaV1.0': C6C6C6U6
        - 'strnaV2.0': C6L15C6L15C6U6T18
        - 'strnaV2.1': C6L15C6L15C6U8T18
        - 'strnaV2.2': C4L15C4L15C4U10T18
        """
        self.fq1 = fq1
        self.fq1_list = fq1.split(',')
        self.n_read = 10000

        # ========== Legacy strnaV patterns (deprecated) ==========
        self.pattern_dict_v1, * \
            _, self.linker_v1_set_list, self.linker_v1_mismatch_list = Barcode.parse_chemistry('strnaV1.0')
        self.pattern_dict_v2, * \
            _, self.linker_1_v2_set_list, self.linker_1_v2_mismatch_list = Barcode.parse_chemistry('strnaV2.0')
        self.pattern_dict_v2, * \
            _, self.linker_4_v2_set_list, self.linker_4_v2_mismatch_list = Barcode.parse_chemistry('strnaV2.1')
        self.pattern_dict_v3, * \
            _, self.linker_v3_set_list, self.linker_v3_mismatch_list = Barcode.parse_chemistry('strnaV2.2')

        # ========== FLV patterns (for VDJ) ==========
        try:
            self.pattern_dict_flv, * \
                _, self.linker_flv_set_list, self.linker_flv_mismatch_list = Barcode.parse_chemistry('flv')
            self.pattern_dict_flv_rna, * \
                _, self.linker_flv_rna_set_list, self.linker_flv_rna_mismatch_list = Barcode.parse_chemistry('flv_rna')
        except:
            # FLV patterns may not be defined
            self.pattern_dict_flv = None
            self.pattern_dict_flv_rna = None

    @utils.add_log
    def check_chemistry(self):
        """check chemistry in the fq1_list"""
        chemistry_list = []
        for fastq1 in self.fq1_list:
            print(fastq1)
            chemistry = self.get_chemistry(fastq1)
            chemistry_list.append(chemistry)
        if len(set(chemistry_list)) != 1:
            Chemistry.check_chemistry.logger.warning('multiple chemistry found!' + str(chemistry_list))
        return chemistry_list

    def seq_chemistry(self, seq):
        """
        Detect chemistry type from a single sequence

        Detection Priority:
        1. FLV patterns (if defined)
        2. New BBV patterns (BBV0, BBV1, BBV2, BBV3)
        3. Legacy strnaV patterns (for backward compatibility)

        Returns: chemistry string or None
        """

        # Check flv_rna first (otherwise may be considered as strnaV2.1.1 and strnaV2.2.1)
        if self.pattern_dict_flv_rna:
            try:
                linker_flv_rna = Barcode.get_seq_str(seq, self.pattern_dict_flv_rna["L"])
                bool_valid, _, _, = Barcode.check_seq_mismatch(
                    [linker_flv_rna], self.linker_flv_rna_set_list, self.linker_flv_rna_mismatch_list)
                if bool_valid:
                    return "flv_rna"
            except:
                pass

        # ========== Check new BBV patterns ==========
        # These patterns should be checked before legacy strnaV patterns
        # to ensure new chemistries are properly detected

        # Try to detect BBV patterns by checking segment lengths and linker positions
        # BBV patterns have characteristic structure: C-L-C-L-C-U
        try:
            # Check if this looks like a BBV pattern (has linker sequences at expected positions)
            # BBV0: C6-L15-C6-L15-C6-U10 (total ~60-70bp)
            # BBV1: C4-L15-C4-L15-C4-U10 (total ~60-70bp)
            # BBV2: C4-space-L15-C4-L15-C4-U10 (variable linker)
            # BBV3: C5-space-L15-C5-L15-C5-U10 (variable linker)

            if len(seq) >= 60:
                # Try BBV2/BBV2.4 pattern (C4-based with displacement)
                # Check both BBV2 and BBV2.4 since they have same linker sequences
                try:
                    pattern_dict_bbv2, *_, linker_bbv2_set_list, linker_bbv2_mismatch_list = Barcode.parse_chemistry('BBV2.4')
                    linker_bbv2 = Barcode.get_seq_str(seq, pattern_dict_bbv2['pattern']['L15']["L"])
                    bool_valid, _, _ = Barcode.check_seq_mismatch(
                        [linker_bbv2], linker_bbv2_set_list, linker_bbv2_mismatch_list)
                    if bool_valid:
                        # Return BBV2.4 as it's the most commonly used variant
                        return "BBV2.4"
                except:
                    pass

                # Try BBV3/BBV3.1 pattern (C5-based with displacement)
                try:
                    pattern_dict_bbv3, *_, linker_bbv3_set_list, linker_bbv3_mismatch_list = Barcode.parse_chemistry('BBV3.1')
                    linker_bbv3 = Barcode.get_seq_str(seq, pattern_dict_bbv3['pattern']['L15']["L"])
                    bool_valid, _, _ = Barcode.check_seq_mismatch(
                        [linker_bbv3], linker_bbv3_set_list, linker_bbv3_mismatch_list)
                    if bool_valid:
                        return "BBV3.1"
                except:
                    pass

                # Try BBV0 pattern (C6-based with displacement, scRNA only)
                try:
                    pattern_dict_bbv0, *_, linker_bbv0_set_list, linker_bbv0_mismatch_list = Barcode.parse_chemistry('BBV0')
                    linker_bbv0 = Barcode.get_seq_str(seq, pattern_dict_bbv0['pattern']['L15']["L"])
                    bool_valid, _, _ = Barcode.check_seq_mismatch(
                        [linker_bbv0], linker_bbv0_set_list, linker_bbv0_mismatch_list)
                    if bool_valid:
                        return "BBV0"
                except:
                    pass

                # Try BBV2.2 pattern (C4-based, fixed linker, different from BBV2/BBV2.4)
                try:
                    pattern_dict_bbv22, *_, linker_bbv22_set_list, linker_bbv22_mismatch_list = Barcode.parse_chemistry('BBV2.2')
                    linker_bbv22 = Barcode.get_seq_str(seq, pattern_dict_bbv22['pattern']['L15']["L"])
                    bool_valid, _, _ = Barcode.check_seq_mismatch(
                        [linker_bbv22], linker_bbv22_set_list, linker_bbv22_mismatch_list)
                    if bool_valid:
                        return "BBV2.2"
                except:
                    pass
        except:
            pass

        # ========== Check legacy strnaV patterns ==========

        # Check legacy strnaV2.0/2.1 patterns
        linker_v2 = Barcode.get_seq_str(seq, self.pattern_dict_v2["L"])
        bool_valid, _, _, = Barcode.check_seq_mismatch(
            [linker_v2], self.linker_1_v2_set_list, self.linker_1_v2_mismatch_list)
        if bool_valid:
            if seq[65:69] == "TTTT":
                return "strnaV2.0.1"
            else:
                return "strnaV2.1.1"

        # Check legacy strnaV3 pattern (C5 based)
        linker_v3 = Barcode.get_seq_str(seq, self.pattern_dict_v3["L"])
        bool_valid, _, _, = Barcode.check_seq_mismatch(
            [linker_v3], self.linker_v3_set_list, self.linker_v3_mismatch_list)
        if bool_valid:
            return "strnaV3.0.1"

        # Check legacy strnaV2.2 pattern
        linker_v2 = Barcode.get_seq_str(seq, self.pattern_dict_v2["L"])
        bool_valid, _, _, = Barcode.check_seq_mismatch(
            [linker_v2], self.linker_4_v2_set_list, self.linker_4_v2_mismatch_list)
        if bool_valid:
            return "strnaV2.2.1"

        # Check flv pattern
        if self.pattern_dict_flv:
            try:
                linker_flv = Barcode.get_seq_str(seq, self.pattern_dict_flv["L"])
                bool_valid, _, _, = Barcode.check_seq_mismatch(
                    [linker_flv], self.linker_flv_set_list, self.linker_flv_mismatch_list)
                if bool_valid:
                    return "flv"
            except:
                pass

        return None

    @utils.add_log
    def get_chemistry(self, fq1):
        results = defaultdict(int)

        with pysam.FastxFile(fq1) as fh:
            for _ in range(self.n_read):
                entry = fh.__next__()
                seq = entry.sequence
                chemistry = self.seq_chemistry(seq)
                if chemistry:
                    results[chemistry] += 1
        # if it is 0, then no other linker types
        if results["strnaV2.2.1"] != 0:
            results["strnaV2.2.1"] += results["strnaV2.1.1"]
        sorted_counts = sorted(results.items(), key=lambda x: x[1], reverse=True)
        self.get_chemistry.logger.info(sorted_counts)

        chemistry, read_counts = sorted_counts[0][0], sorted_counts[0][1]
        percent = float(read_counts) / self.n_read
        if percent < 0.5:
            self.get_chemistry.logger.warning("Valid chemistry read counts percent < 0.5")
        if percent < 0.1:
            self.get_chemistry.logger.error("Valid chemistry read counts percent < 0.1")
            raise Exception(
                'Auto chemistry detection failed! ' + HELP_DICT['chemistry']
            )
        Chemistry.get_chemistry.logger.info(f'chemistry: {chemistry}')

        return chemistry


class Barcode(Step):
    """
    Barcode Demultiplexing and Filtering for scRNA and stRNA modes

    ## Features

    - Demultiplex barcodes with support for both scRNA and stRNA modes
    - Support for displacement matching (space) in BBV0, BBV2, BBV3
    - Filter invalid R1 reads, including:
        - Reads without linker: mismatch > 2 with whitelist linkers
        - Reads without correct barcode: mismatch > 1 with whitelist barcodes
        - Reads without polyT: T bases in polyT region < 10
        - Low quality reads: low sequencing quality in barcode and UMI regions

    ## Chemistry Versions

    - **BBV0**: C6-space (scRNA only, testing version)
    - **BBV1**: C4-fixed (Universal, development version, less important)
    - **BBV2**: C4-space (Universal: scrna + strna, production version)
    - **BBV3**: C5-space (Universal: scrna + strna, future extended version)

    ## Mode-Specific Processing

    ### scRNA Mode (--mode scrna):
    - **Segment-wise validation**: Each C4/C5/C6 segment validated independently
    - **Example (BBV2)**: C4-C4-C4, each C4 segment matches against bclist
    - **All segments must pass**: Only reads with all three segments valid are kept
    - **Supported chemistries**: BBV0, BBV1, BBV2, BBV3
    - **Whitelist**: Uses fragment file (e.g., 208 4bp segments in /data/chemistry/BBV2/bclist)

    ### stRNA Mode (--mode strna):
    - **Complete barcode matching**: Full barcode matched against spatial whitelist
    - **Supported chemistries**: BBV1, BBV2, BBV3
    - **Whitelist**: Uses H5 file with spatial coordinates (barcodeToPos.h5)
    - **Spatial mapping**: Barcode linked to physical position on chip

    ## Multi-Lane/Multi-Sequencing Support

    - **Standard Format**: `{Sample}_S{N}_L{NNN}_R{N}_001.fastq.gz`
      * L001, L002, L003... represent multiple sequencing runs (lanes)
      * Use case: Supplementary sequencing (加测) for increased coverage
      * Processing: Each lane processed independently, then merged

    - **Legacy Format**: `{Sample}_fold{N}_R{N}.fq.gz` (backward compatible)
      * fold1, fold2... mapped to L001, L002...

    - **Input Methods**:
      * Comma-separated: `--fq1 sample_L001_R1.fq.gz,sample_L002_R1.fq.gz`
      * Glob pattern: `--fq1 'sample_L*_R1.fq.gz'`

    - **Output**: Combined results from all lanes for integrated analysis

    ## Output

    - `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads
    - Read name format: `{barcode}_{UMI}_{read ID}`
    - Optional: `{sample}_1.fq(.gz)` if --output_R1 is specified
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # Use new FastqFileParser to auto-detect and parse files
        # Supports both standard format and legacy fold format
        self.fq1_list, self.fq2_list = FastqFileParser.auto_detect_files(
            args.fq1, args.fq2)

        self.fq_number = len(self.fq1_list)
        if self.fq_number != len(self.fq2_list):
            raise Exception('fastq1 and fastq2 do not have same file number!')

        # Parse file information for better naming
        self.fq1_infos = FastqFileParser.parse_file_list(self.fq1_list)
        self.fq2_infos = FastqFileParser.parse_file_list(self.fq2_list)

        # Print file detection summary
        FastqFileParser.print_file_summary(self.fq1_list, self.fq2_list)

        if args.chemistry == 'auto':
            # Use comma-separated for Chemistry detection (legacy format)
            fq1_str = ','.join(self.fq1_list)
            ch = Chemistry(fq1_str)
            self.chemistry_list = ch.check_chemistry()
        else:
            self.chemistry_list = [args.chemistry] * self.fq_number
        self.mismatch = 1
        self.barcode_corrected_num = 0
        self.linker_corrected_num = 0
        self.exactly_matched_num = 0
        self.mismatched_num = 0
        self.total_num = 0
        self.clean_num = 0
        self.no_polyT_num = 0
        self.lowQual_num = 0
        self.no_linker_num = 0
        self.no_barcode_num = 0
        self.barcodes_corrds = dict()
        self.barcode_val_exact_Counter = Counter()
        self.barcode_val_mis_Counter = Counter()
        self.barcode_inval_Counter = Counter()
        self.barcode_qual_Counter = Counter()
        self.umi_qual_Counter = Counter()
        self.pattern = args.pattern
        self.linker = args.linker
        self.whitelist = args.whitelist
        self.lowNum = args.lowNum
        self.lowQual = args.lowQual
        self.filterNoPolyT = args.filterNoPolyT
        self.allowNoLinker = args.allowNoLinker
        self.nopolyT = args.nopolyT  # true == output nopolyT reads
        self.noLinker = args.noLinker
        self.output_R1 = args.output_R1
        self.bool_flv = False
        self.resume = args.resume
        self.mode = args.mode  # scRNA or stRNA
        self.tree = 'trie'  # trie, avl or bk are available
        
        # Parallel processing of multiple files
        self.max_parallel_files = getattr(args, 'max_parallel_files', 3)
        
        # Resuming from breakpoint
        self.progress_tracker = None
        if self.resume:
            self.progress_tracker = ProgressTracker(self.outdir, self.sample)

        # flv_trust4, flv_CR
        if self.assay in ('flv_CR', 'flv_trust4'):
            self.bool_flv = True
            self.barcode_read_Counter = Counter()
            if self.assay == 'flv_trust4':
                if args.match_dir == 'None':
                    raise FileNotFoundError('Match_dir required when running flv_trust4')
                self.match_barcodes = set(
                    utils.get_barcode_from_match_dir(args.match_dir)[0])  # barcode set of flv_rna.
                self.match_num = 0  # record read number match with flv_rna.
                self.match_cbs = set()  # record barcode number match with flv_rna.

        # out file
        if args.gzip:
            self.suffix = ".gz"
        else:
            self.suffix = ""
        self.tmp_dir = f'{self.outdir}/tmp'
        self.out_fq2 = f'{self.out_prefix}_2.fq{self.suffix}'
        self.out_fq1 = f'{self.out_prefix}_1.fq{self.suffix}'
        if self.nopolyT:
            self.nopolyT_1 = f'{self.out_prefix}_noPolyT_1.fq'
            self.nopolyT_2 = f'{self.out_prefix}_noPolyT_2.fq'
        if self.noLinker:
            self.noLinker_1 = f'{self.out_prefix}_noLinker_1.fq'
            self.noLinker_2 = f'{self.out_prefix}_noLinker_2.fq'

        self.open_files()

    @staticmethod
    def get_seq_str_no_exception(seq, sub_pattern_dict):
        """get subseq with intervals in arr and concatenate"""
        return ''.join([seq[item[0]: item[1]] for item in sub_pattern_dict])

    @staticmethod
    def get_seq_str(seq, sub_pattern_dict):
        """
        Get subseq with intervals in arr and concatenate

        Args:
            seq: str
            sub_pattern_dict: [[0, 8], [24, 32], [48, 56]]

        Returns:
            str

        Raise:
            IndexError: if sequence length is not enough
        """
        seq_len = len(seq)
        ans = []
        for item in sub_pattern_dict:
            start, end = item[0], item[1]
            if end > seq_len:
                raise IndexError(f"sequence length is not enough in R1 read: {seq}")
            else:
                ans.append(seq[start:end])
        return ''.join(ans)

    @staticmethod
    def get_seq_list(seq, pattern_dict, abbr):
        return [seq[item[0]: item[1]] for item in pattern_dict[abbr]]

    @staticmethod
    def get_displace_seq_list(seq, pattern_dict, abbr, clen, lk_range):
        def get_fmp(sq, lk_ptn):
            positions = [sq.find(p) for p in lk_ptn]
            valid_positions = [pos for pos in positions if pos != -1]
            if valid_positions:
                fmp = min(valid_positions)
                valid_count = len(valid_positions)
                return fmp, valid_count
            else:
                return -1, 0

        dis_pattern = pattern_dict['pattern']
        linker = pattern_dict['linker']
        sub_seq_len = int(np.ceil(clen / 3))

        # Search for second linker (p2) position
        # For BBV2.4: lk_range=[15,18], search from position 23 to 41
        # to detect L15/L16/L17/L18 patterns
        search_seq = seq[sub_seq_len*2+lk_range[0]:sub_seq_len*2+lk_range[0]+lk_range[1]]
        first_match_position, _ = get_fmp(search_seq, linker['p2'])
        if first_match_position != -1:
            idx = f'L{lk_range[0] + first_match_position}'
        else:
            idx = 'L15'
            # return -1, None

        # Fix: Directly use the detected pattern without any offset
        # dis_pattern[idx] already contains the correct barcode positions for L15/L16/L17/L18
        # No need to apply offset - that was the bug causing reads to be filtered

        # Safety check: if detected pattern doesn't exist, fall back to L15
        if idx not in dis_pattern:
            idx = 'L15'

        return idx, [seq[item[0]: item[1]] for item in dis_pattern[idx][abbr]]

    @staticmethod
    def parse_pattern(pattern):
        pattern_dict = defaultdict(list)
        p = re.compile(r'([CLUNT])(\d+)')
        tmp = p.findall(pattern)
        if not tmp:
            Barcode.parse_pattern.logger.error(f'Invalid pattern: {pattern}')
            sys.exit()
        start = 0
        for item in tmp:
            end = start + int(item[1])
            pattern_dict[item[0]].append([start, end])
            start = end
        return pattern_dict

    @staticmethod
    @utils.add_log
    def parse_displacement_pattern(barcode_dict):
        pattern_dict = defaultdict()
        pattern_dict['pattern'] = defaultdict()
        pattern_dict['linker'] = defaultdict()
        for key, value in barcode_dict['pattern'].items():
            pattern_dict['pattern'][key] = Barcode.parse_pattern(value)
        for key, value in barcode_dict['linker'].items():
            pattern = set()
            for p in value:
                pattern |= Barcode.findall_mismatch(p, n_mismatch=1)
            pattern_dict['linker'][key] = pattern
        return pattern_dict

    @staticmethod
    def get_abbr_len(pattern_dict, abbr):
        length = 0
        for item in pattern_dict[abbr]:
            length += item[1] - item[0]

        return length

    @staticmethod
    def get_scope_bc(chemistry):
        """Return (linker file path, whitelist file path)"""
        try:
            linker_f = glob.glob(f'{ROOT_PATH}/data/chemistry/{chemistry}/linker*')[0]
            whitelist_f = f'{ROOT_PATH}/data/chemistry/{chemistry}/bclist'
        except IndexError:
            return None, None
        return linker_f, whitelist_f

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
    def findall_mismatch(seq, n_mismatch=1, bases='ACGTN'):
        """
        choose locations where there's going to be a mismatch using combinations
        and then construct all satisfying lists using product

        Return:
        all mismatch <= n_mismatch set.
        """
        seq_set = set()
        seq_len = len(seq)
        if n_mismatch > seq_len:
            n_mismatch = seq_len
        for locs in combinations(range(seq_len), n_mismatch):
            seq_locs = [[base] for base in seq]
            for loc in locs:
                seq_locs[loc] = list(bases)
            for poss in product(*seq_locs):
                seq_set.add(''.join(poss))
        return seq_set

    @staticmethod
    def get_mis_overlap(tree, mis_mask_lens, mis_mask, barcode_int, mismatch):
        mis_mask_index = 0
        mis_count = 0
        overlap_iter = None

        # mis_match_list = []
        for mis in range(mismatch):
            while mis_mask_index < mis_mask_lens[mis]:
                mis_barcode_int = barcode_int ^ mis_mask[mis_mask_index]
                mis_mask_index += 1
                if (0, mis_barcode_int) in tree.find(mis_barcode_int, 0):
                    overlap_iter = mis_barcode_int
                    mis_count += 1
                    # mis_match_list.append(mis_barcode_int)
                    if mis_count > 1:
                        return None

        return overlap_iter if mis_count == 1 else None

    @staticmethod
    def get_n_overlap(tree, barcode_string, mismatch):
        n_idx = barcode_string.find('N')
        if n_idx != -1 and mismatch > 0:
            mis_count = 0
            barcode_int = utils.encoder(barcode_string)
            overlap_iter = None

            if (0, barcode_int) in tree.find(barcode_int, 0):
                mis_count += 1
                overlap_iter = barcode_int

            for j in range(1, 4):
                mis_barcode_int = barcode_int ^ np.uint64(j << (len(barcode_string) - 1 - n_idx) * 2)
                if (0, mis_barcode_int) in tree.find(mis_barcode_int, 0):
                    mis_count += 1
                    if mis_count > 1:
                        return None
                    overlap_iter = mis_barcode_int

            if mis_count == 1:
                return overlap_iter
        return None

    @staticmethod
    def get_mismask(barcode_len, mismatch, verbose=False):
        mismask_len = utils.possible_mislen(barcode_len, mismatch)
        mismask_lens = [utils.possible_mislen(barcode_len, i + 1) for i in range(mismatch)]
        mis_mask = np.zeros(mismask_len, dtype=np.uint64)

        index = 0
        if mismatch > 0:
            for i in range(barcode_len):
                for j in range(1, 4):
                    mis_mask_int = j << i * 2
                    mis_mask[index] = mis_mask_int
                    index += 1
            if verbose:
                print(f"1 mismatch mask barcode number: {index}")

        if mismatch == 2:
            mis_mask_set = set()
            for i in range(barcode_len):
                for j in range(1, 4):
                    mis_mask_int1 = j << i * 2
                    for k in range(barcode_len):
                        if k == i:
                            continue
                        for j2 in range(1, 4):
                            mis_mask_int2 = j2 << k * 2
                            mis_mask_int = mis_mask_int1 | mis_mask_int2
                            mis_mask_set.add(mis_mask_int)
            for mis_mask_int in mis_mask_set:
                mis_mask[index] = mis_mask_int
                index += 1
            if verbose:
                print(f"2 mismatch mask barcode number: {index}")

        if mismatch == 3:
            mis_mask_set = set()
            for i in range(barcode_len):
                for j in range(1, 4):
                    mis_mask_int1 = j << i * 2
                    for k in range(barcode_len):
                        if k == i:
                            continue
                        for j2 in range(1, 4):
                            mis_mask_int2 = j2 << k * 2
                            for h in range(barcode_len):
                                if h == k or h == i:
                                    continue
                                for j3 in range(1, 4):
                                    mis_mask_int3 = j3 << h * 2
                                    mis_mask_int = mis_mask_int1 | mis_mask_int2 | mis_mask_int3
                                    mis_mask_set.add(mis_mask_int)
            for mis_mask_int in mis_mask_set:
                mis_mask[index] = mis_mask_int
                index += 1
            if verbose:
                print(f"3 mismatch mask barcode number: {index}")

        if verbose:
            print(f"total mismatch mask barcode number: {index}")
        return mis_mask, mismask_lens

    @staticmethod
    @utils.add_log
    def get_mismatch_dict(seq_list, n_mismatch=1):
        """
        Return:
        mismatch dict. Key: mismatch seq, value: seq in seq_list
        """
        mismatch_dict = {}

        for seq in seq_list:
            seq = seq.strip()
            if seq == '':
                continue
            for mismatch_seq in Barcode.findall_mismatch(seq, n_mismatch):
                mismatch_dict[mismatch_seq] = seq

        return mismatch_dict

    @staticmethod
    def get_c_len(pattern_dict):
        if 'pattern' not in pattern_dict.keys():
            return sum([item[1] - item[0] for item in pattern_dict['C']])
        else:
            c_len = []
            for key, value in pattern_dict['pattern'].items():
                c_len.append(sum([item[1] - item[0] for item in value['C']]))
            c_len = set(c_len)
            return c_len.pop()


    @staticmethod
    def check_seq_mismatch_st(seq, correct_tree, mismask_lens, mismask, mismatch, bit_width=None):
        """
        Perform mismatch correction on the given sequence.

        Args:
            seq (str): The input sequence.
            correct_tree (tree): The correction tree.
            mismask_lens (list): The mismatch mask lengths.
            mismask (ndarray): The mismatch mask.
            mismatch (int): The allowed number of mismatches.
            bit_width (int): Bit width for encoding (2 * barcode_length). If None, uses global BIT_WIDTH.

        Returns:
            tuple: A tuple containing the following values:
                - is_valid (bool): Indicates whether the sequence is valid.
                - is_corrected (bool): Indicates whether the sequence has been corrected.
                - corrected_seq (str): The corrected sequence.
        """
        # Use chemistry-specific bit_width if provided, otherwise fallback to global BIT_WIDTH
        _bit_width = bit_width if bit_width is not None else BIT_WIDTH

        def decode_seq(seq):
            return utils.decoder(np.binary_repr(seq, width=_bit_width))

        bool_valid = True
        bool_corrected = False
        corrected_seq = ''

        if 'N' in seq:
            # Correct sequences with 'N' mismatches
            corrected_seq = Barcode.get_n_overlap(correct_tree, seq, mismatch)
            if corrected_seq:
                # Correction was needed for 'N' bases
                bool_corrected = True
                corrected_seq = decode_seq(corrected_seq)
            else:
                bool_valid = False
                corrected_seq = ''
        else:
            seq_int = utils.encoder(seq)
            match_list = correct_tree.find(seq_int, 0)
            if len(match_list) != 0:
                # Exact match - no correction needed
                bool_corrected = False
                corrected_seq = decode_seq(match_list[0][1])
            elif mismatch > 0:
                # Correct sequences with mismatches
                corrected_seq = Barcode.get_mis_overlap(correct_tree, mismask_lens, mismask, seq_int, mismatch)
                if corrected_seq:
                    # Correction was needed and succeeded
                    bool_corrected = True
                    corrected_seq = utils.decoder(np.binary_repr(corrected_seq, width=_bit_width))
                else:
                    bool_valid = False
                    corrected_seq = ''
            else:
                bool_valid = False
        return bool_valid, bool_corrected, corrected_seq

    @staticmethod
    def check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list):
        """
        Return bool_valid, bool_corrected, corrected_seq
        """
        bool_valid = True
        bool_corrected = False
        corrected_seq = ''
        for index, seq in enumerate(seq_list):
            if seq not in correct_set_list[index]:
                if seq not in mismatch_dict_list[index]:
                    bool_valid = False
                    return bool_valid, bool_corrected, corrected_seq
                else:
                    bool_corrected = True
                    corrected_seq += mismatch_dict_list[index][seq]
            else:
                corrected_seq += seq
        return bool_valid, bool_corrected, corrected_seq

    @staticmethod
    @utils.add_log
    def generate_combinatorial_whitelist(fragment_file, n_segments=3):
        """
        Generate complete whitelist by combining all possible barcode fragments.

        Args:
            fragment_file: Path to file containing barcode fragments (e.g., 208 4bp fragments)
            n_segments: Number of segments to combine (default: 3 for C4L15C4L15C4)

        Returns:
            barcodes_set: Set of all valid complete barcodes
            barcodes_tree: Tree structure for fast matching
        """
        from itertools import product

        # Read fragment list
        fragments, _ = utils.read_one_col(fragment_file)
        fragments = [f.strip() for f in fragments if f.strip()]

        print(f"Generating combinatorial whitelist from {len(fragments)} fragments...")
        print(f"This will create {len(fragments)**n_segments:,} complete barcodes")

        # Generate all combinations
        all_combinations = product(fragments, repeat=n_segments)
        complete_barcodes = [''.join(combo) for combo in all_combinations]

        barcode_len = len(complete_barcodes[0])
        print(f"Generated {len(complete_barcodes):,} {barcode_len}bp barcodes")

        # Build tree structure for fast matching
        if 'trie' == 'trie':
            from celatlas_spatial.tools.trietree import Trie
            barcodes_tree = Trie()
            for barcode in complete_barcodes:
                barcode_int = utils.encoder(barcode)
                barcodes_tree.insert(barcode_int, barcode_int)

        return set(complete_barcodes), barcodes_tree

    def parse_whitelist_file(self, files: list, n_pattern: int, n_mismatch: int, mode='strna', chemistry=None):
        """
        Parse whitelist file based on mode

        Mode-specific processing:
        - scrna mode: Segment-wise validation
          * Reads fragment file (e.g., 208 C4 fragments)
          * Each segment validated independently
          * Example: C4-C4-C4, each C4 matches against bclist

        - strna mode: Complete barcode matching
          * Reads H5 whitelist with spatial coordinates
          * Complete barcode matched against spatial map

        Args:
            files: whitelist file paths
            n_pattern: number of segments in pattern (e.g., 3 for C4-C4-C4)
            n_mismatch: allowed number of mismatch bases
            mode: 'strna' or 'scrna'
            chemistry: chemistry name (e.g., 'BBV3') for determining bit width

        Returns:
            white_set_list: list of barcode sets (one per segment for scrna, one complete set for strna)
            mismatch_list: list of mismatch correction dicts (scrna) or tree structure (strna)
        """
        n_files = len(files)
        if n_files == 1 and n_pattern > 1:
            files = [files[0]] * n_pattern
        elif n_files == 2 and n_pattern > 1:
            files = [files[0]] * (n_pattern - 1) + [files[1]]
        elif n_files != n_pattern:
            sys.exit(f'Invalid whitelist file number: {n_files}')

        white_set_list, mismatch_list = [], []

        if files[0].endswith('.h5'):
            # stRNA mode: H5 whitelist with spatial coordinates
            # Determine bit_width based on chemistry
            from celatlas_spatial import CHEMISTRY_BC_WIDTH
            bit_width = None
            if chemistry and chemistry in CHEMISTRY_BC_WIDTH:
                bc_width = CHEMISTRY_BC_WIDTH[chemistry]
                bit_width = bc_width * 2
                print(f"Using chemistry-specific bit_width: {bit_width} for {chemistry} (barcode: {bc_width}bp)")

            with h5py.File(files[0], 'r') as h5_mask:
                seq_mtx = h5_mask['bpMatrix_1']
                barcodes, barcodes_tree, barcodes_corrds, _ = utils.matrix_decoder(seq_mtx, bit_width=bit_width)
                white_set_list.append(barcodes)
                self.barcodes_corrds = barcodes_corrds
                return white_set_list, barcodes_tree
        else:
            # scrna mode: Segment-wise validation
            # Read fragment file and create whitelist for each segment
            # Example: 208 C4 fragments → 3 sets of 208 fragments (one per C4 position)
            for f in files:
                barcodes, _ = utils.read_one_col(f)
                white_set_list.append(set(barcodes))
                barcode_mismatch_dict = Barcode.get_mismatch_dict(barcodes, n_mismatch)
                mismatch_list.append(barcode_mismatch_dict)
            return white_set_list, mismatch_list

    @staticmethod
    def parse_chemistry(chemistry, mode='scrna'):
        """
        Parse chemistry configuration and return pattern dict and whitelists

        Returns: pattern_dict, barcode_set_list, barcode_mismatch_list, linker_set_list, linker_mismatch_list
        """
        pattern = PATTERN_DICT[chemistry]

        # Handle both simple patterns (str) and displacement patterns (dict)
        if isinstance(pattern, dict):
            # Displacement pattern (BBV0, BBV2, BBV2.4, BBV3, BBV3.1, etc.)
            pattern_dict = Barcode.parse_displacement_pattern(pattern)
        else:
            # Simple pattern (BBV1, strnaV2.2, etc.)
            pattern_dict = Barcode.parse_pattern(pattern)

        linker_file, whitelist_file = Barcode.get_scope_bc(chemistry)

        # Get number of barcode segments (works for both simple and displacement patterns)
        if 'pattern' in pattern_dict:
            # Displacement pattern: get from any L15/L16/L17/L18 entry
            first_pattern_key = list(pattern_dict['pattern'].keys())[0]
            n_pattern = len(pattern_dict['pattern'][first_pattern_key]['C'])
        else:
            # Simple pattern
            n_pattern = len(pattern_dict['C'])

        # Create temporary instance to call instance methods
        temp_instance = Barcode.__new__(Barcode)

        barcode_set_list, barcode_mismatch_list = temp_instance.parse_whitelist_file([whitelist_file],
                                                                            n_pattern=n_pattern,
                                                                            n_mismatch=1,
                                                                            mode=mode,
                                                                            chemistry=chemistry)
        linker_set_list, linker_mismatch_list = temp_instance.parse_whitelist_file([linker_file], n_pattern=1, n_mismatch=2,
                                                                          mode='linker',
                                                                          chemistry=chemistry)

        return pattern_dict, barcode_set_list, barcode_mismatch_list, linker_set_list, linker_mismatch_list

    @staticmethod
    def check_polyT(seq, pattern_dict, min_polyT_count=MIN_T):
        """
        Return:
            True if polyT is found
        """
        seq_polyT = Barcode.get_seq_str(seq, pattern_dict['T'])
        n_polyT_found = seq_polyT.count('T')
        if n_polyT_found >= min_polyT_count:
            return True
        return False

    def open_files(self):
        if self.output_R1 or self.bool_flv:
            self.fh_fq1 = xopen(self.out_fq1, 'w')  
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
        self.fh_fq2.close()

        if self.nopolyT:
            self.fh_nopolyT_fq1.close()
            self.fh_nopolyT_fq2.close()

        if self.noLinker:
            self.fh_nolinker_fq1.close()
            self.fh_nolinker_fq2.close()

    def merge_chunk_files(self, file_index=None):
        """Deprecated: Now directly merge all chunk files at the end, no intermediate steps needed"""
        print("Warning: merge_chunk_files method is deprecated in new workflow")
        print("All chunk files will be merged at the end using merge_all_chunks method")
        pass
    
    def merge_all_chunks(self):
        """Merge all chunk files to final output - new unified merge method"""
        print("Starting to merge all chunk files...")

        # Find all chunk files (extension determined dynamically by suffix)
        # Support both standard lane-based (L001_R1) and legacy (fold1_1) naming
        r1_chunk_pattern_lane = os.path.join(self.tmp_dir, f'{self.sample}_L*_R1_*.fq{self.suffix}')
        r2_chunk_pattern_lane = os.path.join(self.tmp_dir, f'{self.sample}_L*_R2_*.fq{self.suffix}')
        r1_chunk_pattern_legacy = os.path.join(self.tmp_dir, f'{self.sample}_fold*_1_*.fq{self.suffix}')
        r2_chunk_pattern_legacy = os.path.join(self.tmp_dir, f'{self.sample}_fold*_2_*.fq{self.suffix}')

        # Try lane format first, fallback to legacy
        r1_chunks = sorted(glob.glob(r1_chunk_pattern_lane))
        r2_chunks = sorted(glob.glob(r2_chunk_pattern_lane))

        if not r1_chunks:  # Fallback to legacy format
            r1_chunks = sorted(glob.glob(r1_chunk_pattern_legacy))
            r2_chunks = sorted(glob.glob(r2_chunk_pattern_legacy))
        
        print(f"Found {len(r1_chunks)} R1 chunk files")
        print(f"Found {len(r2_chunks)} R2 chunk files")
        
        # Display found files list for debugging
        if r1_chunks:
            print("R1 chunks found:")
            for chunk in r1_chunks[:5]:  # Show only first 5
                print(f"  {chunk}")
            if len(r1_chunks) > 5:
                print(f"  ... and {len(r1_chunks) - 5} more")
        
        if r2_chunks:
            print("R2 chunks found:")
            for chunk in r2_chunks[:5]:  # Show only first 5
                print(f"  {chunk}")
            if len(r2_chunks) > 5:
                print(f"  ... and {len(r2_chunks) - 5} more")
        
        # Filter out existing and non-empty files
        valid_r1_chunks = [f for f in r1_chunks if os.path.exists(f) and os.path.getsize(f) > 0]
        valid_r2_chunks = [f for f in r2_chunks if os.path.exists(f) and os.path.getsize(f) > 0]
        
        print(f"Valid R1 chunks: {len(valid_r1_chunks)}")
        print(f"Valid R2 chunks: {len(valid_r2_chunks)}")
        
        # Merge R1 files
        if (self.output_R1 or self.bool_flv) and valid_r1_chunks:
            print(f"Merging R1 chunks to: {self.out_fq1}")
            self._merge_gzipped_files_properly(valid_r1_chunks, self.out_fq1)
            print(f"Successfully merged {len(valid_r1_chunks)} R1 chunks")
            
            # Clean up R1 chunk files
            for chunk_file in valid_r1_chunks:
                try:
                    os.remove(chunk_file)
                except Exception as e:
                    print(f"Warning: Failed to remove {chunk_file}: {e}")
        
        # Merge R2 files
        if valid_r2_chunks:
            print(f"Merging R2 chunks to: {self.out_fq2}")
            self._merge_gzipped_files_properly(valid_r2_chunks, self.out_fq2)
            print(f"Successfully merged {len(valid_r2_chunks)} R2 chunks")
            
            # Clean up R2 chunk files
            for chunk_file in valid_r2_chunks:
                try:
                    os.remove(chunk_file)
                except Exception as e:
                    print(f"Warning: Failed to remove {chunk_file}: {e}")
        
        # Finally clean up temporary directory
        if os.path.exists(self.tmp_dir):
            try:
                shutil.rmtree(self.tmp_dir)
                print(f"Cleaned up temporary directory: {self.tmp_dir}")
            except Exception as e:
                print(f"Warning: Failed to clean up temp directory: {e}")
        
        print("Chunk file merging completed!")

    def merge_final_files(self, temp_files):
        """Merge multiple temporary files to final output - compatibility method, now uses merge_all_chunks"""
        print("Warning: merge_final_files is deprecated in new workflow")
        print("Using merge_all_chunks method instead of processing temp_files")
        
        # Directly call the new merge method
        self.merge_all_chunks()
    
    def _merge_gzipped_files_properly(self, input_files, output_file):
      
        import gzip
        import tempfile
        from xopen import xopen
        
        valid_files = []
        for input_file in input_files:
            if os.path.exists(input_file) and os.path.getsize(input_file) > 0:
                valid_files.append(input_file)
            elif os.path.exists(input_file):
                print(f"Warning: Empty input file skipped: {input_file}")
            else:
                print(f"Warning: Input file does not exist: {input_file}")
        
        if not valid_files:
            print(f"Warning: No valid files to merge, creating empty output: {output_file}")
            with xopen(output_file, 'wt') as f:
                pass
            return

        
        temp_output = None
        try:
            
            output_dir = os.path.dirname(output_file)
            output_basename = os.path.basename(output_file)
            if output_basename.endswith('.gz'):
                temp_suffix = '.tmp.gz'  
            else:
                temp_suffix = '.tmp'
            temp_fd, temp_output = tempfile.mkstemp(suffix=temp_suffix, dir=output_dir)
            os.close(temp_fd)  
            
            
          
            with xopen(temp_output, 'wt') as outf:
                for i, input_file in enumerate(valid_files):

                    try:
                        with xopen(input_file, 'rt') as inf:

                            line_count = 0
                            for line in inf:
                                outf.write(line)
                                line_count += 1

                    except Exception as e:
                        print(f"Error processing input file {input_file}: {e}")
                        raise
            
            shutil.move(temp_output, output_file)
            temp_output = None  
            
            print(f"Successfully merged {len(valid_files)} files to {output_file}")
            
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
    def process_single_file_pair(self, fq1, fq2, file_index, chemistry, lane_number=None):
        """
        Process a single pair of FASTQ files

        Args:
            fq1: Path to R1 file
            fq2: Path to R2 file
            file_index: 0-based index of this file pair
            chemistry: Chemistry type
            lane_number: Lane number (for standard format naming), defaults to file_index+1
        """
        if lane_number is None:
            lane_number = file_index + 1
        
        stats = {
            'total_num': 0,
            'clean_num': 0,
            'no_polyT_num': 0,
            'lowQual_num': 0,
            'no_linker_num': 0,
            'no_barcode_num': 0,
            'barcode_corrected_num': 0,
            'linker_corrected_num': 0,
            'exactly_matched_num': 0,
            'mismatched_num': 0,
            'barcode_qual_Counter': Counter(),
            'umi_qual_Counter': Counter(),
            'barcode_val_exact_Counter': Counter(),
            'barcode_val_mis_Counter': Counter(),
            'barcode_inval_Counter': Counter()
        }
        
        pattern_dict = defaultdict()
        whitelist_file, whitelist_files = '', []
        
        lowNum = int(self.lowNum)
        lowQual = int(self.lowQual)
        n_pattern = 0
        
        if chemistry == 'strnaV1':
            lowNum = min(0, lowNum)
            lowQual = max(10, lowQual)
        

        bc_pattern = PATTERN_DICT[chemistry]
        if bc_pattern:
            if self.mode == 'strna':
                whitelist_file = self.whitelist
                whitelist_files = whitelist_file.split(',')
                if isinstance(bc_pattern, str):
                    pattern_dict = self.parse_pattern(bc_pattern)
                    n_pattern = len(pattern_dict['C'])
                elif isinstance(bc_pattern, dict):
                    pattern_dict = self.parse_displacement_pattern(bc_pattern)
                    # Get n_pattern from any pattern key (L15/L16/L17/L18)
                    first_pattern_key = list(pattern_dict['pattern'].keys())[0]
                    n_pattern = len(pattern_dict['pattern'][first_pattern_key]['C'])
            elif self.mode == 'scrna':
                linker_file, whitelist_file = self.get_scope_bc(chemistry)
                whitelist_files = [whitelist_file]
                if isinstance(bc_pattern, str):
                    pattern_dict = self.parse_pattern(bc_pattern)
                    n_pattern = len(pattern_dict['C'])
                elif isinstance(bc_pattern, dict):
                    pattern_dict = self.parse_displacement_pattern(bc_pattern)
                    # Get n_pattern from any pattern key (L15/L16/L17/L18)
                    first_pattern_key = list(pattern_dict['pattern'].keys())[0]
                    n_pattern = len(pattern_dict['pattern'][first_pattern_key]['C'])
        else:
            bc_pattern = self.pattern
            if not bc_pattern:
                raise Exception("invalid bc_pattern!")
            
            linker_file = self.linker
            whitelist_file = self.whitelist
            whitelist_files = whitelist_file.split(',')
            pattern_dict = self.parse_pattern(bc_pattern)
            n_pattern = len(pattern_dict['C'])
        
        bool_T = True if 'T' in pattern_dict else False
        bool_L = True if 'L' in pattern_dict and self.mode == 'scrna' else False
        bool_whitelist = (whitelist_file is not None) and whitelist_file != "None"
        C_len = self.get_c_len(pattern_dict)
        linker_length = [15]

        mismask, mismask_lens = self.get_mismask(C_len, self.mismatch)
        if chemistry.startswith('BB'):
            # Check if this chemistry has displacement patterns (BBV2, BBV3, BBV2.4, etc.)
            if 'pattern' in pattern_dict:
                linker_length = [int(lk.split('L')[1]) for lk in pattern_dict['pattern'].keys()]
            # else: BBV1 uses fixed L15, keep linker_length = [15]
        linker_range = [min(linker_length), max(linker_length)]
        
        if bool_whitelist:
            barcode_set_list, barcode_mismatch_list = self.parse_whitelist_file(whitelist_files,
                                                                                n_pattern=n_pattern,
                                                                                n_mismatch=1,
                                                                                mode=self.mode,
                                                                                chemistry=chemistry)
        if bool_L:
            linker_set_list, linker_mismatch_list = self.parse_whitelist_file([linker_file], n_pattern=1,
                                                                              n_mismatch=2,
                                                                              mode='linker',
                                                                              chemistry=chemistry)
        
        # Use lane number for directory naming (standard lane-based style)
        current_tmp_dir = f'{self.tmp_dir}_L{lane_number:03d}'
        os.makedirs(current_tmp_dir, exist_ok=True)
        
        if IS_PARALLEL:
            self.thread = utils.get_chunk_paired_zip((fq1, fq2), current_tmp_dir, self.thread, False)  # 不使用resume
            
            results = []
            pool = Pool(self.thread)
            for process_id in range(self.thread):
                params = {
                    'pid': process_id,
                    'fastq_files': (os.path.join(current_tmp_dir, f'{process_id}', f'{process_id}_1.fq{self.suffix}'),
                                    os.path.join(current_tmp_dir, f'{process_id}', f'{process_id}_2.fq{self.suffix}')),
                    'mode': self.mode,
                    'tmp_dir': current_tmp_dir,
                    'sample': self.sample,
                    'barcode_sets': barcode_set_list if bool_whitelist else None,
                    'barcode_tree': barcode_mismatch_list if bool_whitelist else None,
                    'mismask_lens': mismask_lens,
                    'mismask': mismask,
                    'mismatch': self.mismatch,
                    'output_R1': self.output_R1,
                    'bool_whitelist': bool_whitelist,
                    'lowQual': lowQual,
                    'lowNum': lowNum,
                    'pattern_dict': pattern_dict,
                    'barcode_inval_Counter': Counter(),
                    'barcode_val_exact_Counter': Counter(),
                    'barcode_val_mis_Counter': Counter(),
                    'barcode_qual_Counter': Counter(),
                    'umi_qual_Counter': Counter(),
                    'C_len': C_len,
                    'linker_range': linker_range,
                    'suffix': self.suffix,
                    'bit_width': C_len * 2,  # Chemistry-specific bit width for barcode encoding
                }
                result = pool.apply_async(process_fastq_chunk, args=(params,))
                results.append(result)
            
            pool.close()
            pool.join()
            

            for result in results:
                (n_total, clean, lowQual_count, no_barcode, barcode_corrected, exactly_matched, mismatched,
                 barcode_inval, barcode_val_exact, barcode_val_mis, barcode_qual, umi_qual) = result.get()
                stats['total_num'] += n_total
                stats['clean_num'] += clean
                stats['lowQual_num'] += lowQual_count
                stats['no_barcode_num'] += no_barcode
                stats['barcode_corrected_num'] += barcode_corrected
                stats['exactly_matched_num'] += exactly_matched
                stats['mismatched_num'] += mismatched
                stats['barcode_inval_Counter'].update(barcode_inval)
                stats['barcode_val_exact_Counter'].update(barcode_val_exact)
                stats['barcode_val_mis_Counter'].update(barcode_val_mis)
                stats['barcode_qual_Counter'].update(barcode_qual)
                stats['umi_qual_Counter'].update(umi_qual)
        
        chunk_files_moved = []
        

        os.makedirs(self.tmp_dir, exist_ok=True)
        
        for process_id in range(self.thread):

            id_suffix = utils.int_to_uuid(process_id)
            src_chunk1 = os.path.join(current_tmp_dir, f'{self.sample}_1_{id_suffix}.fq{self.suffix}')
            src_chunk2 = os.path.join(current_tmp_dir, f'{self.sample}_2_{id_suffix}.fq{self.suffix}')
            
            # Use standard lane naming: L001, L002, L003...
            dst_chunk1 = os.path.join(self.tmp_dir, f'{self.sample}_L{lane_number:03d}_R1_{process_id}.fq{self.suffix}')
            dst_chunk2 = os.path.join(self.tmp_dir, f'{self.sample}_L{lane_number:03d}_R2_{process_id}.fq{self.suffix}')

            if (self.output_R1 or self.bool_flv) and os.path.exists(src_chunk1):
                try:
                    shutil.move(src_chunk1, dst_chunk1)
                    chunk_files_moved.append(dst_chunk1)
 
                except Exception as e:
                    print(f"Error moving R1 chunk file {src_chunk1}: {e}")
                    raise
            elif (self.output_R1 or self.bool_flv):
                print(f"Warning: Expected R1 chunk file not found: {src_chunk1}")
            
            if os.path.exists(src_chunk2):
                try:
                    shutil.move(src_chunk2, dst_chunk2)
                    chunk_files_moved.append(dst_chunk2)

                except Exception as e:
                    print(f"Error moving R2 chunk file {src_chunk2}: {e}")
                    raise
            else:
                print(f"Warning: Expected R2 chunk file not found: {src_chunk2}")
        
        if os.path.exists(current_tmp_dir):
            shutil.rmtree(current_tmp_dir)
        
        print(f"Processed Lane {lane_number:03d} (file_index {file_index}): moved {len(chunk_files_moved)} chunk files")
        return stats

    @utils.add_log
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

        self.add_metric(
            name='Exactly Matched Reads',
            value=self.exactly_matched_num,
            total=self.total_num,
            help_info='reads exactly match in valid reads'
        )

        self.add_metric(
            name='Mismatched Reads',
            value=self.mismatched_num,
            total=self.total_num,
            help_info='reads mismatch in valid reads'
        )

        barcode_qual_total = sum(self.barcode_qual_Counter.values())
        if barcode_qual_total > 0:
            BarcodesQ30 = sum([self.barcode_qual_Counter[k] for k in self.barcode_qual_Counter if k >= self.ord2chr(
                30)]) / float(barcode_qual_total) * 100
            BarcodesQ30 = round(BarcodesQ30, 2)
        else:
            BarcodesQ30 = 0.0
        BarcodesQ30_display = f'{BarcodesQ30}%'
        self.add_metric(
            name='Q30 of Barcodes',
            value=BarcodesQ30,
            display=BarcodesQ30_display,
            help_info='percent of barcode base pairs with quality scores over Q30',
        )

        umi_qual_total = sum(self.umi_qual_Counter.values())
        if umi_qual_total > 0:
            UMIsQ30 = sum([self.umi_qual_Counter[k] for k in self.umi_qual_Counter if k >= self.ord2chr(
                30)]) / float(umi_qual_total) * 100
            UMIsQ30 = round(UMIsQ30, 2)
        else:
            UMIsQ30 = 0.0
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
            raise Exception('no valid reads found! please check the --chemistry parameter.' + HELP_DICT['chemistry'])

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

    @utils.add_log
    def run(self):
        """
        Extract barcode and UMI from R1. Filter reads with
            - invalid polyT
            - low quality in barcode and UMI
            - invalid inlinker
            - invalid barcode
        """
        
        completed_file_pairs = set()
        accumulated_stats = {}
        
        if self.progress_tracker and self.progress_tracker.load_progress():
            completed_file_pairs = self.progress_tracker.completed_files
            accumulated_stats = self.progress_tracker.load_stats()
            print(f"Resume: Found {len(completed_file_pairs)} completed file pairs")
        
        file_pairs = []
        for i in range(self.fq_number):
            file_pair = (self.fq1_list[i], self.fq2_list[i])
            if file_pair not in completed_file_pairs:
                file_pairs.append((i, self.fq1_list[i], self.fq2_list[i], self.chemistry_list[i]))
        
        if not file_pairs:
            print("All files already processed, skipping to final merge.")

            self.merge_all_chunks()

            if accumulated_stats:
                self._aggregate_final_stats(accumulated_stats)
            return
        else:

            if self.max_parallel_files > 1 and len(file_pairs) > 1:

                self._run_parallel_files(file_pairs, completed_file_pairs, accumulated_stats)
            else:

                self._run_sequential_files(file_pairs, completed_file_pairs, accumulated_stats)
        
        if accumulated_stats:
            print(f"Aggregating final statistics...")
            self._aggregate_final_stats(accumulated_stats)
        else:
            print("Warning: No accumulated statistics found!")

            self.add_step_metrics()
        
        if self.progress_tracker:
            self.progress_tracker.cleanup()
    
    def _run_parallel_files(self, file_pairs, completed_file_pairs, accumulated_stats):

        print(f"Multi-file parallel processing: {len(file_pairs)} files with max_parallel={self.max_parallel_files}")
        
        temp_files = []
        all_stats = {}

        print("Note: Using sequential processing due to serialization issues")
        for file_index, fq1, fq2, chemistry in file_pairs:
            try:
                # Get lane number from file info
                lane_number = self.fq1_infos[file_index].lane_number or (file_index + 1)
                print(f"Processing file pair {file_index} (Lane {lane_number:03d}): {fq1}")
                stats = self.process_single_file_pair(fq1, fq2, file_index, chemistry, lane_number)
                all_stats[file_index] = stats
                
                completed_file_pairs.add((fq1, fq2))
                
                if self.progress_tracker:
                    self.progress_tracker.save_progress(completed_file_pairs, all_stats)
                
                print(f"✓ Completed file pair {file_index}: {fq1}")
                
            except Exception as e:
                print(f"✗ Failed processing file pair {file_index}: {e}")
                import traceback
                traceback.print_exc()
                continue
        
        print("All file pairs processed, merging all chunks...")
        self.merge_all_chunks()
        
        for stats in all_stats.values():
            for key, value in stats.items():
                if key.endswith('_Counter'):
                    if key not in accumulated_stats:
                        accumulated_stats[key] = Counter()
                    accumulated_stats[key].update(value)
                else:
                    accumulated_stats[key] = accumulated_stats.get(key, 0) + value
    
    def _run_sequential_files(self, file_pairs, completed_file_pairs, accumulated_stats):

        print(f"Sequential processing: {len(file_pairs)} files")
        
        for file_index, fq1, fq2, chemistry in file_pairs:
            try:
                # Get lane number from file info
                lane_number = self.fq1_infos[file_index].lane_number or (file_index + 1)
                print(f"Processing file pair {file_index} (Lane {lane_number:03d}): {fq1}")
                stats = self.process_single_file_pair(fq1, fq2, file_index, chemistry, lane_number)

                for key, value in stats.items():
                    if key.endswith('_Counter'):
                        if key not in accumulated_stats:
                            accumulated_stats[key] = Counter()
                        accumulated_stats[key].update(value)
                    else:
                        accumulated_stats[key] = accumulated_stats.get(key, 0) + value
                
                completed_file_pairs.add((fq1, fq2))
                
                if self.progress_tracker:
                    self.progress_tracker.save_progress(completed_file_pairs, accumulated_stats)
                
                print(f"✓ Completed file pair {file_index}: {fq1}")
                
            except Exception as e:
                print(f"✗ Failed processing file pair {file_index}: {e}")
                continue

        print("All file pairs processed, merging all chunks...")
        self.merge_all_chunks()
    
    def _process_file_wrapper(self, file_index, fq1, fq2, chemistry):
        lane_number = self.fq1_infos[file_index].lane_number or (file_index + 1)
        return self.process_single_file_pair(fq1, fq2, file_index, chemistry, lane_number)
    
    def _aggregate_final_stats(self, accumulated_stats):

        self.total_num = accumulated_stats.get('total_num', 0)
        self.clean_num = accumulated_stats.get('clean_num', 0)
        self.no_polyT_num = accumulated_stats.get('no_polyT_num', 0)
        self.lowQual_num = accumulated_stats.get('lowQual_num', 0)
        self.no_linker_num = accumulated_stats.get('no_linker_num', 0)
        self.no_barcode_num = accumulated_stats.get('no_barcode_num', 0)
        self.barcode_corrected_num = accumulated_stats.get('barcode_corrected_num', 0)
        self.linker_corrected_num = accumulated_stats.get('linker_corrected_num', 0)
        self.exactly_matched_num = accumulated_stats.get('exactly_matched_num', 0)
        self.mismatched_num = accumulated_stats.get('mismatched_num', 0)
        
        self.barcode_qual_Counter = accumulated_stats.get('barcode_qual_Counter', Counter())
        self.umi_qual_Counter = accumulated_stats.get('umi_qual_Counter', Counter())
        self.barcode_val_exact_Counter = accumulated_stats.get('barcode_val_exact_Counter', Counter())
        self.barcode_val_mis_Counter = accumulated_stats.get('barcode_val_mis_Counter', Counter())
        self.barcode_inval_Counter = accumulated_stats.get('barcode_inval_Counter', Counter())
        
        self.add_step_metrics()

def process_fastq_chunk(params):
    pid = params['pid']
    fastq_files = params['fastq_files']
    mode = params['mode']
    tmp_dir = params['tmp_dir']
    sample = params['sample']
    barcode_sets = params['barcode_sets']
    barcode_tree = params['barcode_tree']
    mismask_lens = params['mismask_lens']
    mismask = params['mismask']
    mismatch = params['mismatch']
    output_R1 = params['output_R1']
    bool_whitelist = params['bool_whitelist']
    lowQual = params['lowQual']
    lowNum = params['lowNum']
    pattern_dict = params['pattern_dict']
    barcode_inval_Counter = params['barcode_inval_Counter']
    barcode_val_exact_Counter = params['barcode_val_exact_Counter']
    barcode_val_mis_Counter = params['barcode_val_mis_Counter']
    barcode_qual_Counter = params['barcode_qual_Counter']
    umi_qual_Counter = params['umi_qual_Counter']
    C_len = params['C_len']
    linker_range = params['linker_range']
    suffix = params['suffix']
    bit_width = params.get('bit_width', None)  # Chemistry-specific bit width

    # statistics
    total_num = 0
    clean_num = 0
    lowQual_num = 0
    no_barcode_num = 0
    barcode_corrected_num = 0
    exactly_matched_num = 0
    mismatched_num = 0
    barcode_mismatch_list = barcode_tree

    id_suffix = utils.int_to_uuid(pid)
    sub_out_fq1 = os.path.join(tmp_dir, f"{sample}_1_{id_suffix}.fq{suffix}")
    sub_out_fq2 = os.path.join(tmp_dir, f"{sample}_2_{id_suffix}.fq{suffix}")

    if output_R1:
        sfh_fq1 = xopen(sub_out_fq1, 'w') 
    sfh_fq2 = xopen(sub_out_fq2, 'w')     

    with dnaio.open(fastq_files[0], fastq_files[1]) as f:
        for fq1, fq2 in f:
            header1, seq1, qual1 = fq1.name, fq1.sequence, fq1.qualities
            header2, seq2, qual2 = fq2.name, fq2.sequence, fq2.qualities
            total_num += 1

            # confirm seq1 pattern
            if 'pattern' not in pattern_dict.keys():
                seq_list = Barcode.get_seq_list(seq1, pattern_dict, 'C')
                umi = Barcode.get_seq_str(seq1, pattern_dict['U'])
                C_U_quals_ascii = Barcode.get_seq_str(qual1, pattern_dict['C'] + pattern_dict['U'])
            else:
                idx, seq_list = Barcode.get_displace_seq_list(seq1, pattern_dict, 'C', C_len, linker_range)
                if seq_list is None:
                    no_barcode_num += 1
                    continue
                umi = Barcode.get_seq_str(seq1, pattern_dict['pattern'][idx]['U'])
                C_U_quals_ascii = Barcode.get_seq_str(qual1, pattern_dict['pattern'][idx]['C'] + pattern_dict['pattern'][idx]['U'])

            # lowQual filter
            if lowQual > 0 and Barcode.low_qual(C_U_quals_ascii, lowQual, lowNum):
                lowQual_num += 1
                continue

            # barcode filter
            if bool_whitelist:
                # Mode-specific barcode validation
                if mode == 'scrna':
                    # scRNA mode: Segment-wise validation
                    # Each C4/C5/C6 segment validated independently against bclist
                    # Example: C4-C4-C4 → validate each C4 separately
                    # All three segments must pass validation
                    bool_valid, bool_corrected, corrected_seq = Barcode.check_seq_mismatch(
                        seq_list,              # ['ACGT', 'CGTA', 'TACG']
                        barcode_sets,          # [set_208, set_208, set_208]
                        barcode_mismatch_list  # [mismatch_dict, mismatch_dict, mismatch_dict]
                    )
                elif mode == 'strna':
                    # stRNA mode: Complete barcode matching
                    # Full 12bp barcode matched against spatial H5 whitelist
                    # Barcode linked to spatial coordinates on chip
                    bool_valid, bool_corrected, corrected_seq = Barcode.check_seq_mismatch_st(
                        ''.join(seq_list),  # 'ACGTCGTATACG'
                        barcode_tree,       # Trie tree structure
                        mismask_lens,
                        mismask,
                        mismatch,
                        bit_width=bit_width  # Use chemistry-specific bit width
                    )
                else:
                    raise ValueError(f"Unknown mode: {mode}, must be 'scrna' or 'strna'")

                if not bool_valid:
                    no_barcode_num += 1
                    # barcode_inval_Counter.update([''.join(seq_list)])
                    continue
                elif bool_corrected:
                    # Barcode needed correction (mismatch)
                    barcode_corrected_num += 1
                    mismatched_num += 1
                    # barcode_val_mis_Counter.update([corrected_seq])
                else:
                    # Barcode matched exactly
                    exactly_matched_num += 1
                    # barcode_val_exact_Counter.update([corrected_seq])
                cb = corrected_seq
            else:
                cb = "".join(seq_list)

            clean_num += 1
            barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
            umi_qual_Counter.update(C_U_quals_ascii[C_len:])

            pid_num = f'{pid}_{total_num}'
            seq1 = cb + umi
            qual1 = C_U_quals_ascii
            if output_R1:
                sfh_fq1.write(f'@{cb}_{umi}_{pid_num}\n{seq1}\n+\n{qual1}\n')
            sfh_fq2.write(f'@{cb}_{umi}_{pid_num}\n{seq2}\n+\n{qual2}\n')

    if output_R1:
        sfh_fq1.close()
    sfh_fq2.close()

    return total_num, clean_num, lowQual_num, no_barcode_num, barcode_corrected_num, \
        exactly_matched_num, mismatched_num, barcode_inval_Counter,  barcode_val_exact_Counter, \
        barcode_val_mis_Counter, barcode_qual_Counter, umi_qual_Counter

@utils.add_log
def barcode(args):
    with Barcode(args, display_title='Demultiplexing') as runner:
        runner.run()

def get_opts_barcode(parser, sub_program=True):
    parser.add_argument(
        '--chemistry',
        help='Predefined (pattern, barcode whitelist, linker whitelist) combinations. ' + HELP_DICT['chemistry'],
        choices=list(PATTERN_DICT.keys()),
        default='auto'
    )
    parser.add_argument(
        '--pattern',
        help="""The pattern of R1 reads, e.g. `C8L16C8L16C8L1U12T18`. The number after the letter represents the number 
of bases.  
        - `C`: cell barcode  
        - `L`: linker(common sequences)  
        - `U`: UMI    
        - `T`: poly T""",
    )
    parser.add_argument(
        '--whitelist',
        help='Cell barcode whitelist file path, one cell barcode per line.'
    )
    parser.add_argument(
        '--linker',
        default=None,
        help='Linker whitelist file path, one linker per line.'
    )
    parser.add_argument(
        '--lowQual',
        help="""Default 0. Bases in cell barcode and UMI whose phred value are lower than 
lowQual will be regarded as low-quality bases.""",
        type=int,
        default=0
    )
    parser.add_argument(
        '--lowNum',
        help='The maximum allowed lowQual bases in cell barcode and UMI.',
        type=int,
        default=2
    )
    parser.add_argument(
        '--nopolyT',
        help='Outputs R1 reads without polyT.',
        action='store_true',
    )
    parser.add_argument(
        '--noLinker',
        help='Outputs R1 reads without correct linker.',
        action='store_true',
    )
    parser.add_argument(
        '--filterNoPolyT',
        help="Filter reads without PolyT.",
        action='store_true'
    )
    parser.add_argument(
        '--allowNoLinker',
        help="Allow valid reads without correct linker.",
        action='store_true'
    )
    parser.add_argument(
        '--gzip',
        help="Output gzipped fastq files.",
        action='store_true'
    )
    parser.add_argument(
        '--output_R1',
        help="Output valid R1 reads.",
        action='store_true'
    )
    parser.add_argument(
        '--resume',
        help="Resume demultiplexing.",
        action='store_true'
    )
    parser.add_argument(
        '--mode',
        help="Barcode extraction analysis method, spatial transcriptomics & single-cell omics.",
        choices=['scrna', 'strna'],
        default='strna'
    )
    parser.add_argument(
        '--max_parallel_files',
        help='Maximum number of file pairs to process in parallel (default: 3)',
        type=int,
        default=3
    )
    if sub_program:
        parser.add_argument('--fq1', help='''R1 fastq file(s). Supports multiple input formats:
            - Single file: sample_R1.fq.gz
            - Comma-separated: sample_L001_R1.fq.gz,sample_L002_R1.fq.gz
            - Glob pattern (quoted): 'sample_L*_R1.fq.gz'
            - Standard format: {Sample}_S{N}_L{NNN}_R1_001.fastq.gz
            - Legacy format: {Sample}_fold{N}_R1.fq.gz''', required=True)
        parser.add_argument('--fq2', help='''R2 fastq file(s). Must match --fq1 format.
            See --fq1 for supported formats.''', required=True)
        parser.add_argument('--match_dir', help='Matched scRNA-seq directory, required for flv_trust4')
        parser = s_common(parser)
    return parser


def main():
    parser = argparse.ArgumentParser(description='Celatlas Spatial', formatter_class=ArgFormatter)
    parser.add_argument('-v', '--version', action='version', version='0.1.0')
    subparsers = parser.add_subparsers(dest='subparser_assay')
    subparser_1st = subparsers.add_parser('rna')
    subparser_2nd = subparser_1st.add_subparsers()
    parser_step = subparser_2nd.add_parser('barcode', formatter_class=ArgFormatter)

    sample = 'ST110148_A1'
    chemistry = 'BBV2.4'
    args = get_opts_barcode(parser_step).parse_args([
        '--outdir', f'/mnt/strna/work_project/pipeline/celatlas_spatial/{sample}/01.barcode',
        '--sample', f'{sample}',
        '--fq1', f'/mnt/strna/work_project/rawdata/celescope/BBV2.4/{sample}_1.fq.gz',
        '--fq2', f'/mnt/strna/work_project/rawdata/celescope/BBV2.4/{sample}_2.fq.gz',
        '--chemistry', f'{chemistry}',
        '--pattern', '',
        '--whitelist', f'/mnt/strna/work_project/rawdata/stomics/mask/{sample}.barcodeToPos.h5',
        '--output_R1',
        '--gzip',
        '--lowQual', '0',
        '--lowNum', '2',
        '--thread', '128',
        '--mode', 'strna',
        '--resume',
    ])
    args.subparser_assay = 'rna'
    barcode(args)


if __name__ == "__main__":
    main()

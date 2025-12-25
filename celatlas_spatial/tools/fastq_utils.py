"""
FASTQ File Naming Utilities

Support for industry-standard FASTQ file naming conventions.

Naming Format:
    {SampleName}_S{SampleNumber}_L{LaneNumber}_R{ReadNumber}_001.fastq.gz

Examples:
    Sample_S1_L001_R1_001.fastq.gz
    Sample_S1_L001_R2_001.fastq.gz
    Sample_S1_L002_R1_001.fastq.gz
    Sample_S1_L002_R2_001.fastq.gz

Backward Compatibility:
    Also supports legacy fold naming:
    - Sample_fold1_R1.fq.gz
    - Sample_fold2_R1.fq.gz
"""

import re
import os
import glob
from collections import defaultdict
from typing import List, Tuple, Dict, Optional


class FastqFileInfo:
    """Information parsed from a FASTQ file name"""

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.filename = os.path.basename(filepath)
        self.dirname = os.path.dirname(filepath)

        self.sample_name: Optional[str] = None
        self.sample_number: Optional[int] = None
        self.lane_number: Optional[int] = None
        self.read_number: Optional[int] = None
        self.chunk_number: Optional[int] = None

        self.fold_number: Optional[int] = None

        self.is_standard_format = False
        self.is_legacy_format = False

        self._parse()

    def _parse(self):
        """Parse the filename to extract components"""

        pattern_standard = r'^(.+?)_S(\d+)_L(\d+)_R([12])_(\d+)\.(fastq|fq)(\.gz)?$'
        match = re.match(pattern_standard, self.filename)

        if match:
            self.is_standard_format = True
            self.sample_name = match.group(1)
            self.sample_number = int(match.group(2))
            self.lane_number = int(match.group(3))
            self.read_number = int(match.group(4))
            self.chunk_number = int(match.group(5))
            return


        pattern_legacy = r'^(.+?)_fold(\d+)_([12]|R[12])\.(fastq|fq)(\.gz)?$'
        match = re.match(pattern_legacy, self.filename)

        if match:
            self.is_legacy_format = True
            self.sample_name = match.group(1)
            self.fold_number = int(match.group(2))
            read_str = match.group(3)
            if read_str.startswith('R'):
                self.read_number = int(read_str[1])
            else:
                self.read_number = int(read_str)
            self.lane_number = self.fold_number
            return

        pattern_simple = r'^(.+?)_R?([12])\.(fastq|fq)(\.gz)?$'
        match = re.match(pattern_simple, self.filename)

        if match:
            self.sample_name = match.group(1)
            self.read_number = int(match.group(2))
            self.lane_number = 1  
            return

    def __repr__(self):
        if self.is_standard_format:
            return (f"FastqFileInfo(Standard: {self.sample_name}, "
                   f"S{self.sample_number}, L{self.lane_number:03d}, "
                   f"R{self.read_number}, {self.chunk_number:03d})")
        elif self.is_legacy_format:
            return (f"FastqFileInfo(Legacy: {self.sample_name}, "
                   f"fold{self.fold_number}, R{self.read_number})")
        else:
            return f"FastqFileInfo(Unknown: {self.filename})"

    def get_sort_key(self) -> Tuple:
        """Get sorting key for grouping files"""
        return (self.sample_name or '',
                self.sample_number or 0,
                self.lane_number or 0,
                self.read_number or 0,
                self.chunk_number or 0)


class FastqFileParser:
    """Parse and group FASTQ files"""

    @staticmethod
    def parse_file_list(file_paths: List[str]) -> List[FastqFileInfo]:
        """Parse a list of file paths"""
        return [FastqFileInfo(fp) for fp in file_paths]

    @staticmethod
    def group_by_lane(file_infos: List[FastqFileInfo]) -> Dict[int, List[FastqFileInfo]]:
        """Group files by lane number"""
        groups = defaultdict(list)
        for info in file_infos:
            if info.lane_number is not None:
                groups[info.lane_number].append(info)
        return dict(groups)

    @staticmethod
    def group_by_read(file_infos: List[FastqFileInfo]) -> Dict[int, List[FastqFileInfo]]:
        """Group files by read number (R1/R2)"""
        groups = defaultdict(list)
        for info in file_infos:
            if info.read_number is not None:
                groups[info.read_number].append(info)
        return dict(groups)

    @staticmethod
    def validate_paired_files(r1_files: List[str], r2_files: List[str]) -> bool:
        """Validate that R1 and R2 files are properly paired"""
        r1_infos = FastqFileParser.parse_file_list(r1_files)
        r2_infos = FastqFileParser.parse_file_list(r2_files)

        if len(r1_infos) != len(r2_infos):
            return False

        # Sort both lists by lane/fold number
        r1_infos.sort(key=lambda x: x.lane_number or 0)
        r2_infos.sort(key=lambda x: x.lane_number or 0)

        # Check pairing
        for r1, r2 in zip(r1_infos, r2_infos):
            if r1.sample_name != r2.sample_name:
                return False
            if r1.lane_number != r2.lane_number:
                return False
            if r1.read_number != 1 or r2.read_number != 2:
                return False

        return True

    @staticmethod
    def auto_detect_files(fq1_input: str, fq2_input: str) -> Tuple[List[str], List[str]]:
        """
        Auto-detect and parse FASTQ files from input

        Supports multiple input formats:
        1. Comma-separated file paths
        2. Glob patterns (with wildcards)
        3. Single file paths

        Args:
            fq1_input: R1 file path(s), comma-separated or glob pattern
            fq2_input: R2 file path(s), comma-separated or glob pattern

        Returns:
            Tuple of (R1 file list, R2 file list), sorted by lane number

        Examples:
            # Comma-separated
            fq1 = "ST_S1_L001_R1_001.fq.gz,ST_S1_L002_R1_001.fq.gz"
            fq2 = "ST_S1_L001_R2_001.fq.gz,ST_S1_L002_R2_001.fq.gz"

            # Glob pattern
            fq1 = "ST_S1_L*_R1_001.fq.gz"
            fq2 = "ST_S1_L*_R2_001.fq.gz"

            # Legacy fold format
            fq1 = "ST_fold1_R1.fq.gz,ST_fold2_R1.fq.gz"
            fq2 = "ST_fold1_R2.fq.gz,ST_fold2_R2.fq.gz"
        """

        def expand_input(input_str: str) -> List[str]:
            """Expand input string to file list"""
            # Check if contains wildcard
            if '*' in input_str or '?' in input_str:
                # Glob pattern
                files = sorted(glob.glob(input_str))
                if not files:
                    # Fallback to treating as literal
                    files = [input_str]
                return files
            else:
                # Comma-separated or single file
                return [f.strip() for f in input_str.split(',')]

        # Expand inputs
        r1_files = expand_input(fq1_input)
        r2_files = expand_input(fq2_input)

        # Parse file infos
        r1_infos = FastqFileParser.parse_file_list(r1_files)
        r2_infos = FastqFileParser.parse_file_list(r2_files)

        # Sort by lane number
        r1_infos.sort(key=lambda x: (x.sample_name or '', x.lane_number or 0))
        r2_infos.sort(key=lambda x: (x.sample_name or '', x.lane_number or 0))

        # Return sorted file paths
        r1_sorted = [info.filepath for info in r1_infos]
        r2_sorted = [info.filepath for info in r2_infos]

        return r1_sorted, r2_sorted

    @staticmethod
    def print_file_summary(r1_files: List[str], r2_files: List[str]):
        """Print a summary of detected files"""
        print("\n" + "="*60)
        print("FASTQ File Detection Summary")
        print("="*60)

        r1_infos = FastqFileParser.parse_file_list(r1_files)
        r2_infos = FastqFileParser.parse_file_list(r2_files)

        print(f"\nTotal file pairs: {len(r1_infos)}")

        if r1_infos:
            sample_name = r1_infos[0].sample_name
            print(f"Sample name: {sample_name}")

            if r1_infos[0].is_standard_format:
                print("Format: Multi-lane")
            elif r1_infos[0].is_legacy_format:
                print("Format: Multi-fold")
            else:
                print("Format: Single")

        print("\nFile pairs:")
        for i, (r1, r2) in enumerate(zip(r1_infos, r2_infos)):
            if r1.is_standard_format:
                print(f"  [{i+1}] Lane {r1.lane_number:03d}: {r1.filename} / {r2.filename}")
            elif r1.is_legacy_format:
                print(f"  [{i+1}] Fold {r1.fold_number}: {r1.filename} / {r2.filename}")
            else:
                print(f"  [{i+1}] {r1.filename} / {r2.filename}")

        print("="*60 + "\n")


def convert_fold_to_lane_naming(sample: str, fold_index: int, read: int,
                                 suffix: str = '.fq.gz') -> str:
    """
    Convert legacy fold naming to standard lane naming

    Args:
        sample: Sample name (e.g., 'Sample123')
        fold_index: Fold number (0-based index)
        read: Read number (1 or 2)
        suffix: File suffix (default: '.fq.gz')

    Returns:
        Standard lane-based filename

    Examples:
        >>> convert_fold_to_lane_naming('Sample123', 0, 1)
        'Sample123_S1_L001_R1_001.fq.gz'

        >>> convert_fold_to_lane_naming('Sample123', 1, 2)
        'Sample123_S1_L002_R2_001.fq.gz'
    """
    lane = fold_index + 1
    return f"{sample}_S1_L{lane:03d}_R{read}_001{suffix}"


# Test code
if __name__ == "__main__":
    # Test standard format
    test_files_standard = [
        "Sample_S1_L001_R1_001.fastq.gz",
        "Sample_S1_L002_R1_001.fastq.gz",
        "Sample_S1_L001_R2_001.fastq.gz",
        "Sample_S1_L002_R2_001.fastq.gz",
    ]

    print("Testing standard format:")
    for f in test_files_standard:
        info = FastqFileInfo(f)
        print(f"  {info}")

    # Test legacy format
    test_files_legacy = [
        "Sample_fold1_R1.fq.gz",
        "Sample_fold2_R1.fq.gz",
        "Sample_fold1_R2.fq.gz",
        "Sample_fold2_R2.fq.gz",
    ]

    print("\nTesting legacy format:")
    for f in test_files_legacy:
        info = FastqFileInfo(f)
        print(f"  {info}")

    # Test validation
    r1_files = [test_files_standard[0], test_files_standard[1]]
    r2_files = [test_files_standard[2], test_files_standard[3]]

    print("\nValidation test:")
    is_valid = FastqFileParser.validate_paired_files(r1_files, r2_files)
    print(f"  Paired files valid: {is_valid}")

    # Test naming conversion
    print("\nNaming conversion test:")
    for i in range(3):
        for read in [1, 2]:
            filename = convert_fold_to_lane_naming('Sample', i, read)
            print(f"  fold{i+1}_R{read} -> {filename}")

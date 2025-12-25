# ==========================================
# Barcode Pattern Definitions
# ==========================================
# Pattern Format:
#   - C: Cell barcode segment
#   - L: Linker (common sequences)
#   - U: UMI (Unique Molecular Identifier)
#   - T: Poly-T tail
#
# Chemistry Versions:
#   - BBV0: Testing version, C6-space (scRNA only)
#   - BBV1: Development version, C4-fixed (Universal, less important)
#   - BBV2: Production version, C4-space (Universal: scrna + strna)
#   - BBV3: Extended version, C5-space (Universal: scrna + strna, future)
#
# Application Modes:
#   - scrna mode (--mode scrna): Single-cell RNA sequencing
#     * Segment-wise matching: C4-C4-C4, each segment validated independently
#     * All three segments must match against bclist
#     * Example: BBV2 validates each C4 segment separately
#
#   - strna mode (--mode strna): Spatial transcriptomics
#     * Uses H5 whitelist with spatial coordinates
#     * Complete barcode matching against spatial map
#
# Displacement Matching (space):
#   - Supports dynamic linker length detection (L15, L16, L17, L18)
#   - Used in BBV0, BBV2, BBV3
#   - BBV1 uses fixed linker length (no displacement)
#
# Multi-fold Structure:
#   - fold1, fold2, ... represent multiple sequencing runs
#   - Used for supplementary sequencing (加测) scenarios
#   - Results can be merged for combined analysis
#
# Deprecated Patterns:
#   - strnaV1.0, strnaV2.0, strnaV2.1, strnaV2.2 (legacy, will be removed)
#   - BBV2.2, BBV2.4, BBV3.1 (kept for backward compatibility)
# ==========================================

PATTERN_DICT = {
    'auto': None,

    # ========== BBV0: scRNA only (Testing version) ==========
    # BBV0: C6 with displacement matching (scRNA mode only)
    # Pattern: C6-space-L15-C6-L15-C6-U10
    'BBV0': {
        'pattern': {
            'L15': 'C6L15C6L15C6U10T18',
            'L16': 'C6L16C6L15C6U10T18',
            'L17': 'C6L17C6L15C6U10T18',
            'L18': 'C6L18C6L15C6U10T18',
        },
        'linker': {
            'p1': ['CGACTCACTACACGT'],
            'p2': ['TCGCTGACACGATCG'],
        },
    },

    # ========== BBV1: Universal (Development version, less important) ==========
    # BBV1: C4 with fixed linker (scrna + strna)
    # Pattern: C4-L15-C4-L15-C4-U10
    'BBV1': 'C4L15C4L15C4U10T18',

    # ========== BBV2: Universal (Production version) ==========
    # BBV2: C4 with displacement matching (scrna + strna)
    # Pattern: C4-space-L15-C4-L15-C4-U10
    # scrna mode: Each C4 segment validated independently against bclist
    # strna mode: Complete barcode matched against H5 whitelist
    'BBV2': {
        'pattern': {
            'L15': 'C4L15C4L15C4U10T18',
            'L16': 'C4L16C4L15C4U10T18',
            'L17': 'C4L17C4L15C4U10T18',
            'L18': 'C4L18C4L15C4U10T18',
        },
        'linker': {
            'p1': ['CGACTCACTACACGT'],
            'p2': ['TCGCTGACACGATCG'],
        },
    },

    # ========== BBV3: Universal (Extended version, future) ==========
    # BBV3: C5 with displacement matching (scrna + strna)
    # Pattern: C5-space-L15-C5-L15-C5-U10
    'BBV3': {
        'pattern': {
            'L15': 'C5L15C5L15C5U10T18',
            'L16': 'C5L16C5L15C5U10T18',
            'L17': 'C5L17C5L15C5U10T18',
            'L18': 'C5L18C5L15C5U10T18',
        },
        'linker': {
            'p1': ['CGACTCACTACACGT'],
            'p2': ['TCGCTGACACGATCG'],
        },
    },

    # ========== Legacy Patterns (Deprecated) ==========
    # These patterns are kept for backward compatibility
    # Recommend migrating to BBV0-BBV3
    'strnaV1.0': 'C6C6C6U6',  # Deprecated: Use BBV0 instead
    'strnaV2.0': 'C6L15C6L15C6U6T18',  # Deprecated: Use BBV0 instead
    'strnaV2.1': 'C6L15C6L15C6U8T18',  # Deprecated: Use BBV0 instead
    'strnaV2.2': 'C4L15C4L15C4U10T18',  # Deprecated: Use BBV1 or BBV2 instead

    'BBV2.2': {  # Legacy: Same as BBV1
        'pattern': {
            'L15': 'C4L15C4L15C4U10T18',
        },
        'linker': {
            'p1': ['CGACTCACTACAGGG'],
            'p2': ['TCGGTGACACGATCG'],
        },
    },

    'BBV2.4': {  # Legacy: Similar to BBV2 with different linker
        'pattern': {
            'L15': 'C4L15C4L15C4U10T18',
            'L16': 'C4L16C4L15C4U10T18',
            'L17': 'C4L17C4L15C4U10T18',
            'L18': 'C4L18C4L15C4U10T18',
        },
        'linker': {
            'p1': ['CGACTCACTACACGT'],
            'p2': ['TCGCTGACACGATCG'],
        },
    },

    'BBV3.1': {  # Legacy: Similar to BBV3
        'pattern': {
            'L15': 'C5L15C5L15C5U10T18',
            'L16': 'C5L16C5L15C5U10T18',
            'L17': 'C5L17C5L15C5U10T18',
            'L18': 'C5L18C5L15C5U10T18',
        },
        'linker': {
            'p1': ['CGACTCACTACACGT'],
            'p2': ['TCGCTGACACGATCG'],
        },
    },

    # Custom pattern for user-defined combinations
    'customized': None,
}

# count
RAW_MATRIX_DIR_SUFFIX = ['raw_feature_bc_matrix', 'all_matrix']
FILTERED_MATRIX_DIR_SUFFIX = ['filtered_feature_bc_matrix', 'matrix_10X']
MATRIX_FILE_NAME = 'matrix.mtx'
FEATURE_FILE_NAME = 'features.tsv'
BARCODE_FILE_NAME = 'barcodes.tsv'

# mkref
GENOME_CONFIG = 'celatlas_genome.config'

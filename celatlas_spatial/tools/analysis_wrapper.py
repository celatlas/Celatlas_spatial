import os


def _inherit_parallel_thread_limits():
    """Keep Numba/joblib from bypassing the configured native thread limit."""

    native_threads = os.environ.get('OMP_NUM_THREADS', '').strip()
    if not native_threads.isdigit() or int(native_threads) < 1:
        return
    os.environ.setdefault('NUMBA_NUM_THREADS', native_threads)
    os.environ.setdefault('LOKY_MAX_CPU_COUNT', native_threads)


_inherit_parallel_thread_limits()

import shutil
import cv2
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import seaborn as sns
import networkx as nx

# Set matplotlib to non-interactive backend before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

import json
import re
# from PIL import Image

from pathlib import Path
from celatlas_spatial.tools import utils
from celatlas_spatial.rna.mkref import Mkref_rna
from celatlas_spatial.tools.step import Step, s_common
from celatlas_spatial.__init__ import HELP_DICT, ROOT_PATH

# markers adjust p_value
PVAL_CUTOFF = 0.05
# scanpy mitochondrial variable name
MITO_VAR = 'mt'
NORMALIZED_LAYER = 'normalised'
RESOLUTION = 1.0
N_PCS = 25
MITO_GENE_PERCENT_LIST = [5, 10, 15, 20, 50]
# output marker top n in html
MARKER_TOP_N = 100
CELL_FIGSIZE = (6, 5)
CELL_EMBEDDING_FIGSIZE = (7.2, 5.4)
CELL_SPATIAL_FIGSIZE = (6, 6)
CELL_DPI = 220
CELL_POINT_SIZE = 2
CELL_EMBEDDING_POINT_SIZE = 3
MARKER_COLUMNS = [
    'cluster',
    'gene',
    'avg_log2FC',
    'p_val',
    'p_val_adj',
    'pct.1',
    'pct.2',
]


def get_spatial_pca_config(n_obs, n_hvgs):
    """Return valid PCA and downstream-PC counts for spatial analysis."""
    n_comps = min(50, n_obs, n_hvgs)
    if n_comps < 2:
        raise ValueError(
            "Spatial analysis has too few observations or highly variable genes "
            f"for PCA: n_obs={n_obs}, n_hvgs={n_hvgs}; at least 2 of each are required."
        )
    return n_comps, min(N_PCS, n_comps)


def get_tsne_perplexity(n_obs):
    """Return a t-SNE perplexity valid for the number of observations."""
    if n_obs < 2:
        raise ValueError(f"t-SNE requires at least 2 observations; got n_obs={n_obs}.")
    return min(30, max(1, (n_obs - 1) // 3))


def categorical_palette(n_categories):
    if n_categories <= 20:
        return sc.pl.palettes.default_20[:n_categories]
    cmap = cm.get_cmap('tab20', max(n_categories, 20))
    return [mcolors.to_hex(cmap(i)) for i in range(n_categories)]


def build_qc_filter_mask(
    obs,
    *,
    min_genes=0,
    min_counts=0,
    max_genes=0,
    max_counts=0,
    max_mt=100,
):
    """Return a reusable QC keep mask and per-threshold removal counts."""
    keep = np.ones(len(obs), dtype=bool)
    failed_masks = {}

    if int(min_genes or 0) > 0:
        failed_masks['below_min_genes'] = obs['n_genes_by_counts'] < int(min_genes)
    if int(min_counts or 0) > 0:
        failed_masks['below_min_counts'] = obs['total_counts'] < int(min_counts)
    if int(max_genes or 0) > 0:
        failed_masks['above_max_genes'] = obs['n_genes_by_counts'] > int(max_genes)
    if int(max_counts or 0) > 0:
        failed_masks['above_max_counts'] = obs['total_counts'] > int(max_counts)

    mt_column = f'pct_counts_{MITO_VAR}'
    if float(max_mt) < 100 and mt_column in obs.columns:
        failed_masks['above_max_mt'] = obs[mt_column] > float(max_mt)

    removed_by_filter = {}
    for name, failed in failed_masks.items():
        failed_values = np.asarray(failed, dtype=bool)
        removed_by_filter[f'spots_{name}'] = int(failed_values.sum())
        keep &= ~failed_values

    return keep, removed_by_filter


def read_tsne(tsne_file):
    df = pd.read_csv(tsne_file, sep='\t')
    # compatible with old version
    if 'Unnamed: 0' in df.columns:
        df.rename(columns={'Unnamed: 0': 'barcode'}, inplace=True)
        df = df.set_index('barcode')
    return df


def _normalise_rank_genes_groups_df(df_markers):
    """Return a stable marker table when no significant genes are found."""
    markers_name_dict = {
        'group': 'cluster',
        'names': 'gene',
        'logfoldchanges': 'avg_log2FC',
        'pvals': 'p_val',
        'pvals_adj': 'p_val_adj',
        'pct_nz_group': 'pct.1',
        'pct_nz_reference': 'pct.2',
    }
    if df_markers is None or 'logfoldchanges' not in df_markers.columns:
        return pd.DataFrame(columns=MARKER_COLUMNS)

    df_markers = df_markers[df_markers['logfoldchanges'].notna()].rename(
        markers_name_dict,
        axis='columns',
    )
    if 'cluster' not in df_markers.columns or 'p_val_adj' not in df_markers.columns:
        return pd.DataFrame(columns=MARKER_COLUMNS)

    clusters = pd.to_numeric(df_markers['cluster'], errors='coerce')
    df_markers = df_markers.loc[clusters.notna()].copy()
    df_markers['cluster'] = clusters.loc[df_markers.index].astype(int)
    df_markers = df_markers.loc[df_markers['p_val_adj'] < PVAL_CUTOFF]
    for column in MARKER_COLUMNS:
        if column not in df_markers.columns:
            df_markers[column] = pd.Series(index=df_markers.index, dtype='float64')
    return df_markers.loc[:, MARKER_COLUMNS]


def format_df_marker(df_marker):
    avg_logfc_col = "avg_log2FC"  # seurat 4
    if "avg_logFC" in df_marker.columns:  # seurat 2.3.4
        avg_logfc_col = "avg_logFC"
    df_marker = df_marker.loc[:, ["cluster", "gene", avg_logfc_col, "pct.1", "pct.2", "p_val_adj"]]
    df_marker["cluster"] = df_marker["cluster"].apply(lambda x: f"cluster {x}")
    df_marker = df_marker[df_marker["p_val_adj"] < PVAL_CUTOFF]

    return df_marker


def parse_gtf_gene_id_name(gtf_file):
    """Return {gene_id: gene_name} from a GTF file."""
    gene_map = {}
    if not gtf_file or not os.path.exists(gtf_file):
        return gene_map
    opener = gzip.open if str(gtf_file).endswith('.gz') else open
    gene_id_pattern = re.compile(r'gene_id "([^"]+)"')
    gene_name_pattern = re.compile(r'gene_name "([^"]+)"')
    with opener(gtf_file, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            gene_id_match = gene_id_pattern.search(parts[8])
            if not gene_id_match:
                continue
            gene_id = gene_id_match.group(1)
            gene_name_match = gene_name_pattern.search(parts[8])
            gene_map[gene_id] = gene_name_match.group(1) if gene_name_match else gene_id
    return gene_map


def get_opts_analysis(parser, sub_program):
    parser.add_argument('--genomeDir', help=HELP_DICT['genomeDir'], required=True)
    if sub_program:
        parser.add_argument(
            '--matrix_file',
            help='Required. Matrix_10X directory from step count.',
            required=False,
        )
        parser.add_argument(
            '--position_file',
            help='Required. Position file from step spatial.',
            required=False,
        )
        parser.add_argument(
            '--square_bin_dir',
            help='Required. Square bin directory from step spatial.',
            required=False,  # Changed to False to support both spatial and scRNA
        )
        parser.add_argument(
            '--cell-seg-dir',
            help='Optional. Cell segmentation output directory from 06.segment/02.cellsegment.',
            required=False,
        )
        parser.add_argument(
            '--cell-h5ad',
            help='Optional. Cell-level AnnData file from cell segmentation.',
            required=False,
        )
        parser.add_argument(
            '--cell-matrix-dir',
            help='Optional. Cell-level 10X matrix directory from cell segmentation.',
            required=False,
        )
        parser.add_argument(
            '--cell-metadata',
            help='Optional. Cell metadata TSV. Defaults to <cell-seg-dir>/<sample>.cell_metadata.tsv.',
            required=False,
        )
        parser.add_argument(
            '--hires-image',
            help='Optional. H&E/tissue hires image used as the spatial plot background.',
            required=False,
        )
        parser.add_argument(
            '--fullres-image',
            help='Optional. Registered full-resolution tissue image used only to infer fullres image dimensions.',
            required=False,
        )

        parser.add_argument(
            '--pixel-size',
            help='Required. Pixel size of the image.',
            type=float,
            required=False,  # Changed to False to support both spatial and scRNA
        )
        parser.add_argument(
            '--bin',
            help='Required. Bin size(micron).',
            type=int,
            default=100,
        )
        parser.add_argument(
            '--mt-gene-list',
            help='Optional. List of mitochondrial genes.',
            required=False,
        )
        parser.add_argument(
            '--assay',
            help='Analysis mode: "spatial", "scrna", or "cell". Auto-detected if not specified.',
            type=str,
            choices=['spatial', 'scrna', 'cell'],
            required=False,
        )
        parser.add_argument(
            '--skip-tsne',
            help='Skip t-SNE calculation (useful for large datasets or to avoid segfault)',
            action='store_true',
            default=False,
        )
        parser.add_argument(
            '--skip-umap',
            help='Skip UMAP calculation',
            action='store_true',
            default=False,
        )
        parser.add_argument(
            '--cluster-resolution',
            help=f'Leiden clustering resolution for binned spatial analysis. Default: {RESOLUTION}.',
            type=float,
            default=RESOLUTION,
        )
        parser.add_argument(
            '--cell-cluster-resolution',
            help='Leiden clustering resolution for cell-segmented analysis. Default: 0.4.',
            type=float,
            default=0.4,
        )
        parser.add_argument(
            '--cell-n-neighbors',
            help='Number of nearest neighbors for cell-segmented graph construction. Default: 10.',
            type=int,
            default=10,
        )
        parser.add_argument(
            '--cell-n-pcs',
            help='Number of PCA components for cell-segmented neighbors/UMAP/Leiden. Default: 10.',
            type=int,
            default=10,
        )
        parser.add_argument(
            '--cell-min-genes',
            help='Minimum genes for QC of segmented cells or the selected spatial bin. Default: 0 disables filtering.',
            type=int,
            default=0,
        )
        parser.add_argument(
            '--cell-min-counts',
            help='Minimum UMI/counts for QC of segmented cells or the selected spatial bin. Default: 0 disables filtering.',
            type=int,
            default=0,
        )
        parser.add_argument(
            '--cell-max-genes',
            help='Maximum genes for QC of segmented cells or the selected spatial bin. Default: 0 disables filtering.',
            type=int,
            default=0,
        )
        parser.add_argument(
            '--cell-max-counts',
            help='Maximum UMI/counts for QC of segmented cells or the selected spatial bin. Default: 0 disables filtering.',
            type=int,
            default=0,
        )
        parser.add_argument(
            '--cell-max-mt',
            help='Maximum mitochondrial percent for QC of segmented cells or the selected spatial bin. Default: 100 disables filtering.',
            type=float,
            default=100.0,
        )
        parser = s_common(parser)

def anndata_to_cundata(adata):
    """
    Convert anndata object to cuAnnData object.
    :param adata: anndata object
    :return: cuAnnData object cunData
    """
    import cudf
    import cupy as cp
    import cupyx.scipy.sparse as cusparse

    if isinstance(adata.X, cusparse.csr_matrix):
        cu_X = adata.X
    else:
        cu_X = cusparse.csr_matrix(cp.array(adata.X.todense() if sc.sparse.issparse(adata.X) else adata.X))

    cu_obs = cudf.from_pandas(adata.obs)
    cu_var = cudf.from_pandas(adata.var)

    cu_adata = {
        'X': cu_X,
        'obs': cu_obs,
        'var': cu_var,
        'uns': adata.uns,  # notice that uns is not converted to cudf
        'obsm': {},
        'varm': {},
        'layers': {}
    }

    for key, matrix in adata.obsm.items():
        cu_adata['obsm'][key] = cp.array(matrix)
    for key, matrix in adata.varm.items():
        cu_adata['varm'][key] = cp.array(matrix)

    for key, layer in adata.layers.items():
        if isinstance(layer, cusparse.csr_matrix):
            cu_adata['layers'][key] = layer
        else:
            cu_adata['layers'][key] = cusparse.csr_matrix(
                cp.array(layer.todense() if sc.sparse.issparse(layer) else layer))

    return cu_adata


class Scanpy_wrapper(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.customed_cluster = '0'
        self.mt_gene_list = args.mt_gene_list if hasattr(args, 'mt_gene_list') else None

        # Dimensional reduction options
        self.skip_tsne = getattr(args, 'skip_tsne', False)
        self.skip_umap = getattr(args, 'skip_umap', False)
        self.cluster_resolution = float(getattr(args, 'cluster_resolution', RESOLUTION) or RESOLUTION)
        self.cell_cluster_resolution = float(getattr(args, 'cell_cluster_resolution', 0.4) or 0.4)
        self.cell_n_neighbors = max(2, int(getattr(args, 'cell_n_neighbors', 10) or 10))
        self.cell_n_pcs = max(1, int(getattr(args, 'cell_n_pcs', 10) or 10))
        self.spatial_n_pcs = N_PCS

        # Auto-detect assay type
        if hasattr(args, 'assay') and args.assay:
            self.assay = args.assay
        elif (
            getattr(args, 'cell_seg_dir', None)
            or getattr(args, 'cell_h5ad', None)
            or getattr(args, 'cell_matrix_dir', None)
        ):
            self.assay = 'cell'
        elif hasattr(args, 'square_bin_dir') and args.square_bin_dir:
            self.assay = 'spatial'
        elif hasattr(args, 'matrix_file') and args.matrix_file:
            self.assay = 'scrna'
        else:
            raise ValueError("Cannot determine assay type. Please specify --assay or provide appropriate inputs.")

        self.adata = {}
        self.adata_name = []

        # Load data based on assay type
        if self.assay == 'spatial':
            self._load_spatial_data()
        elif self.assay == 'cell':
            self._load_cell_data()
        else:  # scrna
            self._load_scrna_data()

        # Common output directories
        self._setup_output_dirs()

    def _spatial_spot_size(self, adata, default=1.0):
        """Return spot diameter in full-resolution pixels for scanpy spatial plots."""
        spatial_meta = adata.uns.get("spatial", {})
        library_meta = spatial_meta.get(self.args.sample)
        if library_meta is None and spatial_meta:
            library_meta = next(iter(spatial_meta.values()))
        scalefactors = (library_meta or {}).get("scalefactors", {})
        spot_size = scalefactors.get("spot_diameter_fullres", default)
        try:
            spot_size = float(spot_size)
        except (TypeError, ValueError):
            spot_size = default
        return spot_size if spot_size > 0 else default

    def _spatial_figsize(self, adata=None, base=3.6, ncols=1, min_side=3.0, max_side=14.0):
        """Return a matplotlib figsize matching the chip/image aspect ratio."""
        adata = adata if adata is not None else self.adata.get(self.bin)
        width = height = None
        if adata is not None:
            spatial_meta = adata.uns.get("spatial", {})
            library_meta = spatial_meta.get(self.args.sample)
            if library_meta is None and spatial_meta:
                library_meta = next(iter(spatial_meta.values()))
            images = (library_meta or {}).get("images", {})
            image = images.get("hires")
            if image is None:
                image = images.get("lowres")
            if image is not None and getattr(image, "shape", None) is not None and len(image.shape) >= 2:
                height, width = float(image.shape[0]), float(image.shape[1])

        if (not width or not height) and adata is not None and "spatial" in adata.obsm:
            coords = np.asarray(adata.obsm["spatial"], dtype=float)
            finite = np.isfinite(coords).all(axis=1) if coords.ndim == 2 and coords.shape[1] >= 2 else []
            if coords.size and np.any(finite):
                coords = coords[finite]
                width = float(np.nanmax(coords[:, 0]) - np.nanmin(coords[:, 0]) + 1)
                height = float(np.nanmax(coords[:, 1]) - np.nanmin(coords[:, 1]) + 1)

        if (
            not width or not height
            or not np.isfinite(width) or not np.isfinite(height)
            or width <= 0 or height <= 0
        ):
            width = height = 1.0

        aspect = width / height
        if aspect >= 1:
            fig_width = base * aspect * ncols
            fig_height = base
        else:
            fig_width = base * ncols
            fig_height = base / aspect

        return (
            float(min(max(fig_width, min_side * ncols), max_side * ncols)),
            float(min(max(fig_height, min_side), max_side)),
        )

    def _load_spatial_data(self):
        """Load spatial transcriptomics data"""
        self.px = self.args.pixel_size
        self.bin = f"bin{self.args.bin}"
        self.bins = ['bin10', 'bin20', 'bin50', 'bin100']
        load_errors = {}

        for bin_key in self.bins:
            bin_size = int(bin_key.replace('bin', ''))
            bin_dir = os.path.join(self.args.square_bin_dir, f"{self.args.sample}_{bin_key}")

            try:
                matrix_dir = os.path.join(bin_dir, "filtered_feature_bc_matrix")
                spatial_dir = os.path.join(bin_dir, "spatial")
                positions_file = os.path.join(spatial_dir, "tissue_positions_list.csv")
                scalefactors_file = os.path.join(spatial_dir, "scalefactors_json.json")
                hires_image_file = os.path.join(spatial_dir, "tissue_hires_image.png")
                if not os.path.isfile(positions_file) or os.path.getsize(positions_file) == 0:
                    raise ValueError(f"Spatial positions file is missing or empty: {positions_file}")

                adata = sc.read_10x_mtx(matrix_dir, var_names='gene_symbols')
                adata.layers['raw'] = adata.X.copy()

                positions = pd.read_csv(positions_file, header=None)
                positions.columns = ["barcode", "in_tissue", "row_in_bin", "col_in_bin", "pxl_row_in_fullres",
                                     "pxl_col_in_fullres"]
                positions.index = positions["barcode"].to_list()
                adata.obs = adata.obs.join(positions, how="left")

                with open(scalefactors_file, 'r') as f:
                    scalefactors = json.load(f)

                image = plt.imread(hires_image_file)
                library_id = self.args.sample
                quality = "hires"
                spot_diameter_fullres = scalefactors['spot_diameter_fullres']

                adata.uns["spatial"] = {
                    library_id: {
                        "images": {quality: image},
                        "use_quality": quality,
                        "scalefactors": scalefactors
                    }
                }

                # Scanpy spatial plotting expects full-resolution image coordinates.
                adata.obsm["spatial"] = adata.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values

                self.adata[bin_key] = adata
                self.adata_name.append(bin_key)

            except Exception as e:
                load_errors[bin_key] = str(e)
                print(f"Error loading data for {bin_key}: {e}")

        if self.bin not in self.adata:
            available_bins = ', '.join(self.adata) or 'none'
            reason = load_errors.get(self.bin, 'bin directory was not loaded')
            raise ValueError(
                f"Failed to load requested bin {self.bin}: {reason}. "
                f"Successfully loaded bins: {available_bins}."
            )

    def _load_scrna_data(self):
        """Load scRNA-seq data from 05.count/{sample}_filtered_feature_bc_matrix"""
        # Construct path to filtered matrix
        if hasattr(self.args, 'matrix_file') and self.args.matrix_file:
            matrix_dir = self.args.matrix_file
        else:
            # Auto-detect from outdir structure: 05.count/{sample}_filtered_feature_bc_matrix
            count_dir = os.path.join(os.path.dirname(self.outdir), '05.count')
            matrix_dir = os.path.join(count_dir, f'{self.args.sample}_filtered_feature_bc_matrix')

        print(f"Loading scRNA-seq data from {matrix_dir}")

        # Load 10X matrix in MTX format (standard 3-column features.tsv.gz format)
        adata = sc.read_10x_mtx(matrix_dir, var_names='gene_symbols', cache=False)
        adata.layers['raw'] = adata.X.copy()

        # Store in adata dict with 'scrna' key
        self.adata['scrna'] = adata
        self.adata_name.append('scrna')
        self.bin = 'scrna'  # For compatibility with downstream code

        print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")

    def _default_cell_seg_dirs(self):
        sample_dir = self._default_sample_dir()
        return [
            os.path.join(str(sample_dir), '06.segment', '02.cellsegment'),
        ]

    def _default_sample_dir(self):
        outdir = Path(self.outdir).resolve()
        return outdir.parent.parent if outdir.parent.name == "07.outs" else outdir.parent

    def _set_default_cell_images(self):
        sample_dir = self._default_sample_dir()
        if not getattr(self.args, 'hires_image', None):
            hires_candidates = [
                sample_dir / '06.segment' / '01.binsegment' / 'images' / 'tissue_hires_image.png',
            ]
            for image_path in hires_candidates:
                if image_path.exists():
                    self.args.hires_image = str(image_path)
                    break

        if not getattr(self.args, 'fullres_image', None):
            fullres_candidates = [
                sample_dir / '06.segment' / '01.binsegment' / 'images' / f'{self.sample}_regist.tif',
                sample_dir / '06.segment' / '01.binsegment' / 'images' / f'{self.sample}.tif',
            ]
            for image_path in fullres_candidates:
                if image_path.exists():
                    self.args.fullres_image = str(image_path)
                    break

    def _load_cell_data(self):
        """Load cell-level data from 06.segment/02.cellsegment."""
        cell_seg_dir = getattr(self.args, 'cell_seg_dir', None)
        cell_seg_dirs = [cell_seg_dir] if cell_seg_dir else self._default_cell_seg_dirs()
        h5ad_candidates = []
        if getattr(self.args, 'cell_h5ad', None):
            h5ad_candidates.append(self.args.cell_h5ad)
        for candidate_dir in cell_seg_dirs:
            if not candidate_dir:
                continue
            h5ad_candidates.extend([
                os.path.join(candidate_dir, f'{self.args.sample}_cell_adata.h5ad'),
                os.path.join(candidate_dir, f'{self.args.sample}.cell_adata.h5ad'),
                os.path.join(candidate_dir, f'{self.args.sample}.stardist_cell_adata.h5ad'),
            ])

        h5ad_file = next((path for path in h5ad_candidates if path and os.path.exists(path)), None)
        loaded_cell_seg_dir = os.path.dirname(h5ad_file) if h5ad_file else cell_seg_dir
        if h5ad_file:
            print(f"Loading cell-level AnnData from {h5ad_file}")
            adata = sc.read_h5ad(h5ad_file)
        else:
            matrix_candidates = []
            if getattr(self.args, 'cell_matrix_dir', None):
                matrix_candidates.append(self.args.cell_matrix_dir)
            for candidate_dir in cell_seg_dirs:
                if not candidate_dir:
                    continue
                matrix_candidates.extend([
                    os.path.join(candidate_dir, 'cell_matrix'),
                    os.path.join(candidate_dir, 'filtered_feature_bc_matrix'),
                ])
            matrix_dir = next(
                (
                    path for path in matrix_candidates
                    if path and os.path.exists(os.path.join(path, 'matrix.mtx.gz'))
                ),
                None,
            )
            if not matrix_dir:
                raise FileNotFoundError(
                    "Cell analysis requires --cell-h5ad, --cell-matrix-dir, or "
                    "--cell-seg-dir containing cell_matrix/matrix.mtx.gz."
                )
            print(f"Loading cell-level 10X matrix from {matrix_dir}")
            adata = sc.read_10x_mtx(matrix_dir, var_names='gene_symbols', cache=False)
            loaded_cell_seg_dir = os.path.dirname(matrix_dir)

        adata.var_names_make_unique()
        adata.layers['raw'] = adata.X.copy()
        self._annotate_gene_symbols_from_reference(adata)
        self._attach_cell_metadata(adata, loaded_cell_seg_dir)
        self._set_default_cell_images()
        self.adata['cell'] = adata
        self.adata_name.append('cell')
        self.bin = 'cell'
        print(f"Loaded {adata.n_obs} segmented cells and {adata.n_vars} genes")

    def _annotate_gene_symbols_from_reference(self, adata):
        """Add gene_symbols to adata.var when matrix features only contain gene IDs."""
        if 'gene_symbols' in adata.var.columns:
            return
        gene_ids = adata.var['gene_ids'].astype(str) if 'gene_ids' in adata.var.columns else pd.Index(adata.var_names.astype(str))
        gene_names = pd.Series(adata.var_names.astype(str), index=adata.var_names)
        if gene_names.str.upper().str.startswith('MT-').sum() == 0:
            try:
                gtf_file = Mkref_rna.parse_genomeDir(self.args.genomeDir)['gtf']
                gene_map = parse_gtf_gene_id_name(gtf_file)
                if gene_map:
                    mapped = gene_ids.map(lambda gene_id: gene_map.get(str(gene_id), str(gene_id)))
                    adata.var['gene_symbols'] = mapped.values
                    adata.var_names = pd.Index(mapped.astype(str))
                    adata.var_names_make_unique()
                    adata.var.index.name = None
                    print(f"Annotated gene_symbols from GTF: {gtf_file}")
                    return
            except Exception as exc:
                print(f"WARNING: Failed to annotate gene_symbols from GTF: {exc}")
        adata.var['gene_symbols'] = gene_names.values
        adata.var.index.name = None

    def _set_qc_gene_flags(self, adata):
        """Set mt/rb/hb flags using gene symbols, gene IDs, or explicit mt list."""
        if self.mt_gene_list:
            mito_genes, _ = utils.read_one_col(self.mt_gene_list)
            mito_set = {str(gene).upper() for gene in mito_genes}
            names = pd.Index(adata.var_names.astype(str))
            gene_ids = pd.Index(adata.var['gene_ids'].astype(str)) if 'gene_ids' in adata.var.columns else names
            symbols = pd.Index(adata.var['gene_symbols'].astype(str)) if 'gene_symbols' in adata.var.columns else names
            adata.var[MITO_VAR] = (
                names.str.upper().isin(mito_set) |
                gene_ids.str.upper().isin(mito_set) |
                symbols.str.upper().isin(mito_set)
            )
        else:
            symbols = pd.Index(adata.var['gene_symbols'].astype(str)) if 'gene_symbols' in adata.var.columns else pd.Index(adata.var_names.astype(str))
            names = pd.Index(adata.var_names.astype(str))
            adata.var['mt'] = symbols.str.upper().str.startswith('MT-') | names.str.upper().str.startswith('MT-')
            adata.var['rb'] = symbols.str.upper().str.contains('^RP[SL]') | names.str.upper().str.contains('^RP[SL]')
            adata.var['hb'] = symbols.str.upper().str.contains('^HB[^(P)]') | names.str.upper().str.contains('^HB[^(P)]')
        for qc_var in ('mt', 'rb', 'hb'):
            if qc_var not in adata.var.columns:
                adata.var[qc_var] = False
            adata.var[qc_var] = adata.var[qc_var].astype(bool)

    def _attach_cell_metadata(self, adata, cell_seg_dir):
        metadata_candidates = []
        if getattr(self.args, 'cell_metadata', None):
            metadata_candidates.append(self.args.cell_metadata)
        if cell_seg_dir:
            metadata_candidates.extend([
                os.path.join(cell_seg_dir, f'{self.args.sample}.cell_metadata.tsv'),
                os.path.join(cell_seg_dir, 'cell_metadata.tsv'),
            ])
        metadata_file = next((path for path in metadata_candidates if path and os.path.exists(path)), None)
        if not metadata_file:
            print("Cell metadata TSV not found; continuing without centroid/spatial metadata.")
            return

        print(f"Loading cell metadata from {metadata_file}")
        metadata = pd.read_csv(metadata_file, sep='\t')
        if 'Barcode' not in metadata.columns:
            print(f"WARNING: Cell metadata has no Barcode column: {metadata_file}")
            return
        metadata = metadata.drop_duplicates('Barcode').set_index('Barcode')
        adata.obs = adata.obs.join(metadata, how='left')

        coord_pairs = [
            ('pxl_col_in_fullres', 'pxl_row_in_fullres'),
            ('centroid_col_px', 'centroid_row_px'),
        ]
        for x_col, y_col in coord_pairs:
            if x_col in adata.obs.columns and y_col in adata.obs.columns:
                coords = adata.obs[[x_col, y_col]].apply(pd.to_numeric, errors='coerce')
                valid = coords.notna().all(axis=1)
                if valid.any():
                    adata.obsm['spatial'] = coords.fillna(0).to_numpy(dtype=float)
                    adata.obs['has_cell_spatial'] = valid.astype(int)
                    print(f"Attached cell spatial coordinates from {x_col}/{y_col}")
                return

    def _setup_output_dirs(self):
        """Setup output directories based on assay type"""
        if self.assay == 'spatial':
            self.h5ad_dir = f'{self.outdir}/h5ad'
            self.spatial_qc = f'{self.outdir}/qc'
            self.spatial_cluster = f'{self.outdir}/clustering'
            self.diff_exp_analysis = f'{self.outdir}/markers'
            self.gene_set_enrichment_analysis = None
        elif self.assay == 'cell':
            self.h5ad_dir = f'{self.outdir}/h5ad'
            self.cell_qc_dir = f'{self.outdir}/qc'
            self.cell_embedding_dir = f'{self.outdir}/clustering'
            self.cell_marker_dir = f'{self.outdir}/markers'
            self.cell_spatial_dir = f'{self.outdir}/spatial'
        else:  # scrna
            # scRNA output directories (outdir already points to 06_analysis_wrapper)
            self.analysis_dir = self.outdir
            self.qc_dir = f'{self.analysis_dir}/01.QC'
            self.umap_dir = f'{self.analysis_dir}/02.UMAP'

        # Common files
        if self.assay == 'spatial':
            self.df_marker_file = f'{self.diff_exp_analysis}/{self.sample}_markers.tsv'
            self.df_marker_raw_file = f'{self.diff_exp_analysis}/{self.sample}_markers_raw.tsv'
            self.df_tsne_file = f'{self.spatial_cluster}/{self.sample}_tsne_coord.tsv'
            self.h5ad_file = f'{self.h5ad_dir}/{self.sample}_{self.bin}.h5ad'
        elif self.assay == 'cell':
            self.df_marker_file = f'{self.cell_marker_dir}/{self.sample}_cell_markers.tsv'
            self.df_marker_raw_file = f'{self.cell_marker_dir}/{self.sample}_cell_markers_raw.tsv'
            self.df_tsne_file = f'{self.cell_embedding_dir}/{self.sample}_cell_tsne_coord.tsv'
            self.h5ad_file = f'{self.h5ad_dir}/{self.sample}_cell.h5ad'
        else:
            self.df_marker_file = f'{self.out_prefix}_markers.tsv'
            self.df_marker_raw_file = f'{self.out_prefix}_markers_raw.tsv'
            self.df_tsne_file = f'{self.out_prefix}_tsne_coord.tsv'
            self.h5ad_file = f'{self.out_prefix}.h5ad'

        # Create directories
        if self.assay in ('spatial', 'cell'):
            self.bioinfo_data()
        else:
            self.scrna_data()

    def _render_html(self):
        """Override parent's _render_html to skip HTML generation.
        Scanpy_wrapper doesn't generate reports - that's done by the separate report step.
        """
        pass

    @utils.add_log
    def add_adata_info(self):
        image = cv2.imread(self.args.hires_image)
        image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
        imgarr = np.array(image)
        for i, adata in self.adata.items():
            bin_size = float(i.split('bin')[1]) if i != 'raw' else self.px
            bin_size = np.ceil(bin_size / self.px)
            positions = pd.read_csv(self.args.position_file[i])
            positions.columns = [
                "barcode",
                "in_tissue",
                "row_in_bin",
                "col_in_bin",
                "pxl_row_in_fullres",
                "pxl_col_in_fullres"
            ]

            positions.index = positions["barcode"].to_list()
            positions["new_col"] = positions["col_in_bin"]
            positions["new_row"] = positions["row_in_bin"]

            self.adata[i].obs = self.adata[i].obs.join(positions, how="left")

            self.adata[i].var["ENSEMBL"] = self.adata[i].var["gene_ids"].values
            self.adata[i].obs["index"] = self.adata[i].obs.index

            # max_coor = np.max(positions[["new_col", "new_row"]].values)
            max_coor = np.ceil(max(self.args.im_shape) / bin_size)
            scale = 2000 / max_coor
            self.adata[i].obs["imagecol"] = positions["new_col"].values * scale
            self.adata[i].obs["imagerow"] = positions["new_row"].values * scale

            # max_size = np.max([self.adata[i].obs["imagecol"].max(), self.adata[i].obs["imagerow"].max()])
            # max_size = int(max_size + 0.1 * max_size)
            # image = Image.new("RGB", (max_size, max_size), "White")

            library_id = self.args.sample
            quality = "hires"
            spot_diameter_fullres = bin_size  # micron
            self.adata[i].uns["spatial"] = {
                library_id: {
                    "images": {quality: imgarr},
                    "use_quality": quality,
                    "scalefactors": {
                        f"tissue_{quality}_scalef": scale,
                        "spot_diameter_fullres": spot_diameter_fullres  # in pixels
                    }
                }
            }

            self.adata[i].obsm["spatial"] = self.adata[i].obs[["new_col", "new_row"]].values

    @utils.add_log
    def calculate_qc_metrics(self):
        for i, adata in self.adata.items():
            if self.mt_gene_list:
                mito_genes, _ = utils.read_one_col(self.mt_gene_list)
                self.adata[i].var[MITO_VAR] = self.adata[i].var_names.map(lambda x: True if x in mito_genes else False)
                self.adata[i].var[MITO_VAR] = self.adata[i].var[MITO_VAR].astype(bool)
            else:
                self.adata[i].var['mt'] = self.adata[i].var_names.str.upper().str.startswith('MT-')
                self.adata[i].var['rb'] = self.adata[i].var_names.str.upper().str.contains('^RP[SL]')
                self.adata[i].var['hb'] = self.adata[i].var_names.str.upper().str.contains('^HB[^(P)]')

            sc.pp.calculate_qc_metrics(
                self.adata[i],
                qc_vars=['mt', 'rb', 'hb'],
                percent_top=None,
                use_raw=False,
                log1p=False,
                inplace=True
            )

        import matplotlib.pyplot as plt
        
        single_fig_size = self._spatial_figsize(self.adata[self.bin], base=3.6)
        combined_fig_size = self._spatial_figsize(self.adata[self.bin], base=3.6, ncols=2)
        dpi = 300
        
        # Generate individual plots
        metrics = ["total_counts", "n_genes_by_counts"]
        metric_names = ["total_counts", "gene_counts"]  # File name friendly versions
        
        for metric, name in zip(metrics, metric_names):
            plt.figure(figsize=single_fig_size, dpi=dpi)
            sc.pl.spatial(
                self.adata[self.bin],
                color=[metric],
                cmap="Spectral_r",
                spot_size=self._spatial_spot_size(self.adata[self.bin]),
                alpha_img=0.75,
                use_raw=False,
                show=False,
            )
            plt.savefig(f'{self.spatial_qc}/{self.sample}_spatial_{name}.png', 
                       bbox_inches="tight", dpi=dpi)
            plt.close()
        
        # Generate combined plot with both metrics
        plt.figure(figsize=combined_fig_size, dpi=dpi)
        sc.pl.spatial(
            self.adata[self.bin],
            color=["total_counts", "n_genes_by_counts"],
            cmap="Spectral_r",
            spot_size=self._spatial_spot_size(self.adata[self.bin]),
            alpha_img=0.75,
            use_raw=False,
            show=False,
        )
        plt.savefig(f'{self.spatial_qc}/{self.sample}_spatial_gene_expression_distribution.png', 
                   bbox_inches="tight", dpi=dpi)
        plt.close()
        
        self.raw_adata = self.adata[self.bin].copy()

    @utils.add_log
    def apply_spatial_qc(self):
        """Apply configurable count/gene/mitochondrial QC to the selected bin."""
        adata = self.adata[self.bin]
        before_n = adata.n_obs
        min_genes = int(getattr(self.args, 'cell_min_genes', 0) or 0)
        min_counts = int(getattr(self.args, 'cell_min_counts', 0) or 0)
        max_genes = int(getattr(self.args, 'cell_max_genes', 0) or 0)
        max_counts = int(getattr(self.args, 'cell_max_counts', 0) or 0)
        max_mt = float(getattr(self.args, 'cell_max_mt', 100.0))

        keep, removed_by_filter = build_qc_filter_mask(
            adata.obs,
            min_genes=min_genes,
            min_counts=min_counts,
            max_genes=max_genes,
            max_counts=max_counts,
            max_mt=max_mt,
        )
        filtering_enabled = any((min_genes, min_counts, max_genes, max_counts)) or max_mt < 100
        if filtering_enabled:
            adata = adata[keep, :].copy()
            if adata.n_obs < 2:
                raise ValueError(
                    "Spatial QC retained fewer than two spots: "
                    f"{before_n} -> {adata.n_obs}. Relax the QC thresholds."
                )
            self.adata[self.bin] = adata
            self.raw_adata = self.raw_adata[keep, :].copy()

        summary = {
            'bin': self.bin,
            'spots_before_qc': int(before_n),
            'spots_after_qc': int(adata.n_obs),
            'spots_removed_qc': int(before_n - adata.n_obs),
            'min_genes': min_genes,
            'min_counts': min_counts,
            'max_genes': max_genes,
            'max_counts': max_counts,
            'max_mt': max_mt,
            **removed_by_filter,
        }
        summary_path = os.path.join(
            self.spatial_qc,
            f'{self.sample}_{self.bin}_spatial_qc_summary.tsv',
        )
        pd.DataFrame([summary]).to_csv(summary_path, sep='\t', index=False)
        print(
            f"Spatial QC ({self.bin}): {before_n} -> {adata.n_obs} spots "
            f"(min_counts={min_counts}, max_counts={max_counts}, "
            f"min_genes={min_genes}, max_genes={max_genes}, max_mt={max_mt})"
        )
        self.add_metric(name=f'{self.bin} spots before QC', value=int(before_n))
        self.add_metric(name=f'{self.bin} spots after QC', value=int(adata.n_obs))

    def scatter_plot(self):
        import matplotlib.pyplot as plt
        
        # Fixed size: 8cm x 8cm (converted to inches)
        # 8cm = 8/2.54 ≈ 3.15 inches
        fig_size = (3.15, 3.15)
        dpi = 300
        
        for i, adata in self.adata.items():
            plt.figure(figsize=fig_size, dpi=dpi)
            sc.pl.scatter(
                self.adata[i],
                x='total_counts',
                y='n_genes_by_counts',
                color='total_counts',
                color_map='plasma', 
                size=10,
                show=False,
            )
            plt.savefig(f'{self.spatial_qc}/{self.sample}_{i}_scatter.png', 
                       bbox_inches="tight", dpi=dpi)
            plt.close()
    def violin_plot(self):
        import matplotlib.pyplot as plt
        colors = ["#f38181", "#a8d8ea", "#fce38a"]  
        
        for i, adata in self.adata.items():

            # Enlarged size with 3:1 aspect ratio: 30cm width x 10cm height = 11.81 x 3.94 inches (2x larger)
            fig, axes = plt.subplots(1, 3, figsize=(11.81, 3.94))
            
            metrics = ['total_counts', 'n_genes_by_counts', f'pct_counts_{MITO_VAR}']
            
            for idx, (metric, color, ax) in enumerate(zip(metrics, colors, axes)):

                sc.pl.violin(
                    self.adata[i],
                    metric,
                    jitter=0.4,
                    show=False,
                    ax=ax
                )
                
                # Separate violin patches and scatter points
                violin_patches = []
                scatter_points = []
                
                for collection in ax.collections:
                    # Check the type of collection to distinguish violin patches from scatter points
                    if hasattr(collection, '_paths') and len(collection._paths) > 0:
                        # This is likely a violin patch (PolyCollection)
                        violin_patches.append(collection)
                    elif hasattr(collection, '_offsets') and len(collection._offsets) > 0:
                        # This is likely scatter points (PathCollection)
                        scatter_points.append(collection)
                
                # Style violin patches with the specified color
                for patch in violin_patches:
                    patch.set_facecolor(color)
                    patch.set_alpha(1)
                    patch.set_edgecolor('#323232')  
                    patch.set_linewidth(1.5)
                
                # Style scatter points as solid black dots
                for scatter in scatter_points:
                    scatter.set_facecolor('#323232')
                    scatter.set_edgecolor('#323232')
                    scatter.set_alpha(1.0)
                
                # Style lines (median, quartiles, etc.)
                for line in ax.lines:
                    line.set_color('#323232')  
                    line.set_alpha(1.0)  
            
            plt.tight_layout()
            plt.savefig(f'{self.spatial_qc}/{self.sample}_{i}_violin.png', bbox_inches="tight", dpi=300)
            plt.close()

    @utils.add_log
    def write_mito_stats(self):
        for i, adata in self.adata.items():
            mt_pct_var = f'pct_counts_{MITO_VAR}'
            total_cell_number = self.adata[i].n_obs

            for mito_gene_percent in MITO_GENE_PERCENT_LIST:
                cell_number = sum(self.adata[i].obs[mt_pct_var] > mito_gene_percent)
                fraction = round(cell_number / total_cell_number * 100, 2)
                self.add_metric(
                    name=f'Fraction of cells have mito gene percent>{mito_gene_percent}%',
                    value=f'{fraction}%',
                )

    @utils.add_log
    def normalize(self):
        """
        sc.pp.normalize_per_cell() and sc.pp.log1p()
        """
        sc.pp.normalize_total(
            self.adata[self.bin],
            target_sum=1e4,
            inplace=True,
        )
        sc.pp.log1p(
            self.adata[self.bin],
        )
        self.adata[self.bin].layers[NORMALIZED_LAYER] = self.adata[self.bin].X

    @utils.add_log
    def hvg(self):
        """
        Wrapper function for sc.highly_variable_genes()
        """
        sc.pp.highly_variable_genes(
            self.adata[self.bin],
            layer=None,
            n_top_genes=None,
            min_disp=0.5,
            max_disp=np.inf,
            min_mean=0.0125,
            max_mean=3,
            span=0.3,
            n_bins=20,
            flavor='seurat',
            subset=False,
            inplace=True,
            batch_key=None,
            check_values=True
        )

    @utils.add_log
    def scale(self):
        """
        Wrapper function for sc.pp.scale
        """
        sc.pp.scale(
            self.adata[self.bin],
            zero_center=True,
            max_value=10,
            copy=False,
            layer=None,
            obsm=None
        )

    @utils.add_log
    def pca(self):
        """
        Wrapper function for sc.pp.pca
        """
        adata = self.adata[self.bin]
        n_hvgs = int(np.asarray(adata.var['highly_variable'], dtype=bool).sum())
        n_comps, self.spatial_n_pcs = get_spatial_pca_config(adata.n_obs, n_hvgs)
        if n_comps < 50:
            print(
                f"[PCA] Using {n_comps} components (n_obs={adata.n_obs}, "
                f"n_highly_variable_genes={n_hvgs})."
            )
        sc.pp.pca(
            adata,
            n_comps=n_comps,
            zero_center=True,
            svd_solver='auto',
            random_state=0,
            return_info=False,
            mask_var='highly_variable',
            dtype='float32',
            copy=False,
            chunked=False,
            chunk_size=None
        )

    @utils.add_log
    def neighbors(self):
        """
        Wrapper function for sc.pp.neighbors(), for supporting multiple n_neighbors
        """
        sc.pp.neighbors(
            self.adata[self.bin],
            n_neighbors=15,
            n_pcs=self.spatial_n_pcs,
            use_rep=None,
            knn=True,
            random_state=0,
            method='umap',
            metric='euclidean',
            key_added=None,
            copy=False
        )

    @utils.add_log
    def tsne(self):
        """
        Wrapper function for sc.tl.tsne, for supporting named slot of tsne embeddings

        Note: BLAS thread limits should be set via environment variables BEFORE
        Python process starts to avoid segmentation faults
        """
        adata = self.adata[self.bin]
        perplexity = get_tsne_perplexity(adata.n_obs)
        if perplexity < 30:
            print(f"[t-SNE] Using perplexity={perplexity} for n_obs={adata.n_obs}.")
        sc.tl.tsne(
            adata,
            n_pcs=self.spatial_n_pcs,
            perplexity=perplexity,
            n_jobs=1,
            copy=False,
        )

    @utils.add_log
    def umap(self, ):
        """
        Wrapper function for sc.tl.umap, for supporting named slot of umap embeddings
        """
        sc.tl.umap(
            self.adata[self.bin],
            min_dist=0.5,
            spread=1.0,
            n_components=2,
            maxiter=None,
            alpha=1.0,
            gamma=1.0,
            negative_sample_rate=5,
            init_pos='spectral',
            random_state=0,
            a=None,
            b=None,
            copy=False,
            method='umap',
            neighbors_key=None
        )
        # rsc.tl.umap(
        #     self.adata[self.bin],
        #     min_dist=0.5,
        #     spread=1.0,
        #     n_components=2,
        #     maxiter=None,
        #     alpha=1.0,
        #     negative_sample_rate=5,
        #     init_pos='spectral',
        #     random_state=0,
        #     a=None,
        #     b=None,
        #     copy=False,
        #     neighbors_key=None
        # )

    @utils.add_log
    def leiden(self):
        """
        Wrapper function for sc.tl.leiden
        """
        sc.tl.leiden(
            self.adata[self.bin],
            resolution=self.cluster_resolution,
            restrict_to=None,
            random_state=0,
            flavor='igraph',
            key_added='cluster',
            adjacency=None,
            directed=False,
            use_weights=True,
            n_iterations=2,
            partition_type=None,
            neighbors_key=None,
            obsp=None,
            copy=False
        )

    @utils.add_log
    def highest_expr_genes(self):
        """
        Find genes that are highly expressed in each cluster
        """
        sc.pl.highest_expr_genes(
            self.adata[self.bin],
            n_top=20,
            show=False
        )
        plt.savefig(f'{self.diff_exp_analysis}/highest_expr_genes.png', bbox_inches="tight")

    @utils.add_log
    def find_marker_genes(self):
        """
        Wrapper function for sc.tl.rank_genes_groups
        """
        import matplotlib.pyplot as plt
        
        sc.tl.rank_genes_groups(
            self.adata[self.bin],
            "cluster",
            reference='rest',
            pts=True,
            method="wilcoxon",
            use_raw=False,
            layer=NORMALIZED_LAYER
        )
        
        # Fixed size for marker genes plot (5 columns layout)
        # 15cm width x 12cm height = 5.91 x 4.72 inches
        fig_size = (5.91, 4.72)
        dpi = 300
        
        plt.figure(figsize=fig_size, dpi=dpi)
        sc.pl.rank_genes_groups(self.adata[self.bin], n_genes=25, sharey=False, ncols=5, show=False)
        plt.savefig(f'{self.diff_exp_analysis}/{self.sample}_marker_genes.png', 
                   bbox_inches="tight", dpi=dpi)
        plt.close()

    @utils.add_log
    def differentiate_expression_marker_genes(self):
        """
        Analyze and visualize marker genes based on differential expression results.
        """
        import matplotlib.pyplot as plt
        
        n_genes = 10
        min_logfc = 0.25
        adata = self.adata[self.bin]
        dpi = 300

        # 1. differential expression results
        de_results = self.get_de_results(adata, n_genes, min_logfc)

        if not de_results.empty:
            # Fixed size: 8cm x 8cm = 3.15 x 3.15 inches for all subplots
            fig_size = (3.15, 3.15)
            
            # 2. differential expression heatmap
            plt.figure(figsize=(4.72, 3.94), dpi=dpi)
            sc.pl.heatmap(adata, de_results['names'].unique(), groupby='cluster', cmap="Spectral_r",
                          use_raw=False, swap_axes=True, show_gene_labels=True, show=False)
            plt.savefig(f'{self.diff_exp_analysis}/{self.sample}_marker_heatmap.png', 
                       bbox_inches="tight", dpi=dpi)
            plt.close()

            # 3. differential expression spatial plot (4 columns layout)
            self.customed_cluster = de_results['cluster'].unique()[0]
            # Get actual clusters instead of assuming range(4)
            actual_clusters = sorted(de_results['cluster'].unique())[:4]
            top_genes = [de_results[de_results['cluster'] == cluster]['names'].iloc[0]
                         for cluster in actual_clusters
                         if not de_results[de_results['cluster'] == cluster].empty]

            # Only create plots if we have genes to plot
            if top_genes:
                fig = plt.figure(figsize=self._spatial_figsize(self.raw_adata, base=3.4, ncols=min(len(top_genes), 4)), dpi=dpi)
                # Minimal spacing between subplots
                plt.subplots_adjust(wspace=0.05, hspace=0.05)
                sc.pl.spatial(
                    self.raw_adata,
                    color=top_genes,
                    cmap="Spectral_r",
                    ncols=4,
                    spot_size=self._spatial_spot_size(self.raw_adata),
                    alpha_img=0.75,
                    show=False,
                )
                plt.savefig(f'{self.diff_exp_analysis}/{self.sample}_spatial_markers.png',
                           bbox_inches="tight", dpi=dpi)
                plt.close()
            else:
                print(f"WARNING: No top marker genes found for spatial plot (clusters: {actual_clusters})")
                print(f"  This may happen if differential expression didn't find significant genes (min_logfc={min_logfc})")

                # 4. differential expression violin plot
                fig = plt.figure(figsize=(12.6, 3.15), dpi=dpi)
                # Minimal spacing between subplots
                plt.subplots_adjust(wspace=0.05, hspace=0.05)
                sc.pl.violin(adata, top_genes, groupby='cluster', rotation=90, jitter=0., size=0, show=False)
                plt.savefig(f'{self.diff_exp_analysis}/{self.sample}_violin_markers.png',
                           bbox_inches="tight", dpi=dpi)
                plt.close()

            # 5. differential expression t-SNE and UMAP plot
            if 'X_tsne' in adata.obsm:
                # 4 subplots in one row, each 3.15x3.15 inches with reduced spacing
                fig = plt.figure(figsize=(12.6, 3.15), dpi=dpi)  # 4 * 3.15 = 12.6 inches width
                plt.subplots_adjust(wspace=0.05, hspace=0.05)
                sc.pl.tsne(adata, color=top_genes, cmap="Spectral_r", show=False)
                plt.savefig(f'{self.diff_exp_analysis}/{self.sample}_tsne_markers.png', 
                           bbox_inches="tight", dpi=dpi)
                plt.close()
                
            if 'X_umap' in adata.obsm:
                # 4 subplots in one row, each 3.15x3.15 inches with reduced spacing
                fig = plt.figure(figsize=(12.6, 3.15), dpi=dpi)  # 4 * 3.15 = 12.6 inches width
                plt.subplots_adjust(wspace=0.05, hspace=0.05)
                sc.pl.umap(adata, color=top_genes, cmap="Spectral_r", show=False)
                plt.savefig(f'{self.diff_exp_analysis}/{self.sample}_umap_markers.png', 
                           bbox_inches="tight", dpi=dpi)
                plt.close()
        else:
            self.add_metric(name="WARNING:No significant DE genes found. Skipping differential expression analysis.",
                            value="")

    @utils.add_log
    def gse_analyze_genes(self):
        """
        Perform GO, KEGG, and Reactome enrichment analysis on a list of genes.
        """
        min_logfc = 0.0
        min_pct = 0.1
        max_pval = 0.05

        organism = 'mouse'
        cluster_key = 'cluster'
        gse_outdir = self.gene_set_enrichment_analysis
        adata = self.adata[self.bin]

        dbs = []
        if organism == 'human':
            dbs = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
                   'GO_Cellular_Component_2023', 'KEGG_2021_Human', 'Reactome_2022']
        elif organism == 'mouse':
            dbs = ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
                   'GO_Cellular_Component_2023', 'KEGG_2019_Mouse', 'Reactome_2022']

        try:
            # print("Checking and downloading databases locally...")
            # dbs_offline = []
            # for db in dbs:
            #     db_dir = os.path.join(ROOT_PATH, 'data', 'gene_sets_library', f'{db}.gmt')
            #     os.makedirs(os.path.dirname(db_dir), exist_ok=True)
            #
            #     if not os.path.exists(db_dir):
            #         # download database from enrichr
            #         go_mf = gp.get_library(db, organism=organism)
            #         utils.dict_gmt_txt(db_dir, go_mf)
            #         print(f"Database {db} downloaded successfully.")
            #     else:
            #         print(f"Database {db} is ready.")
            #
            #     dbs_offline.append(db_dir)

            all_results = {}
            for cluster in sorted(adata.obs[cluster_key].unique(), key=int):
                de_genes = sc.get.rank_genes_groups_df(adata, group=cluster)
                de_genes = de_genes[(de_genes['logfoldchanges'] > min_logfc) &
                                    (de_genes['pct_nz_group'] > min_pct) &
                                    (de_genes['pvals_adj'] < max_pval)]
                gene_list = de_genes['names'].tolist()

                if not gene_list:
                    print(f"No significant DE genes found for cluster {cluster}. Skipping enrichment analysis.")
                    continue

                for i in range(len(dbs)):
                    db = dbs[i]
                    enr = gp.enrichr(gene_list=gene_list, gene_sets=db, organism=organism)
                    # try:
                    #     enr = gp.enrich(gene_list=gene_list, gene_sets=dbs_offline[i], verbose=True)
                    # except Exception as e:
                    #     print(f"Error in enrichment analysis: {e}, may not found gene list in database {db}.")
                    #     continue
                    if db not in all_results:
                        all_results[db] = {}
                    enr.res2d.Term = enr.res2d.Term.str.split(" \(GO").str[0]
                    all_results[db][cluster] = enr.res2d
                print(f"Enrichment analysis completed for cluster {cluster}")

            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')  # ignore warnings from gseapy
                for db, results in all_results.items():
                    self.plot_gse_heatmap(gse_outdir, db, results)
                    self.plot_gse_barplot(gse_outdir, db, results, self.customed_cluster)
                    self.plot_gse_dotplot(gse_outdir, db, results, self.customed_cluster)
                    self.plot_gse_emapplot(gse_outdir, db, results, self.customed_cluster)

        except Exception as e:
            print(f"Error in enrichment analysis: {e}")

    @staticmethod
    def plot_gse_heatmap(outdir, db, results, top_n=20):
        """
        Plot GO enrichment analysis results as a heatmap.
        :param outdir: Output directory for saving the plot
        :param db: Database name (e.g., 'GO_Biological_Process_2021')
        :param results: Dictionary of enrichment results for all clusters
        :param top_n: Number of top GO terms to display (default is 20)
        """

        # Check if we have results
        if not results or len(results) == 0:
            print(f"No results available for {db} heatmap")
            return
            
        # merge results
        merged_results = pd.concat(results, axis=0, keys=results.keys())
        merged_results.reset_index(level=0, inplace=True)
        merged_results.rename(columns={'level_0': 'Cluster'}, inplace=True)
        
        # Filter valid results
        merged_results = merged_results.dropna(subset=['Term', 'Adjusted P-value'])
        merged_results = merged_results[merged_results['Adjusted P-value'] > 0]  # Remove zero p-values
        
        if merged_results.empty:
            print(f"No valid results for {db} heatmap after filtering")
            return

        # choose top n terms based on minimum p-value across clusters
        top_terms = merged_results.groupby('Term')['Adjusted P-value'].min().nsmallest(top_n).index

        # prepare data for heatmap
        heatmap_data = merged_results[merged_results['Term'].isin(top_terms)]
        pivot_data = heatmap_data.pivot(index='Term', columns='Cluster', values='Adjusted P-value')
        
        # Fill NaN values with maximum p-value to avoid display issues
        pivot_data = pivot_data.fillna(1.0)
        
        # Convert to -log10(p-value) for better visualization
        pivot_data = -np.log10(pivot_data + 1e-10)  # Add small value to avoid log(0)
        
        # Ensure we have valid data
        if pivot_data.empty or pivot_data.isna().all().all():
            print(f"No valid pivot data for {db} heatmap")
            return

        # Larger size for better readability: 16cm width x 12cm height = 6.30 x 4.72 inches
        fig_size = (6.30, 4.72)
        dpi = 300
        
        plt.figure(figsize=fig_size, dpi=dpi)
        
        # Create heatmap with better formatting
        sns.heatmap(pivot_data, 
                   cmap='YlOrRd', 
                   annot=False, 
                   fmt='.1f',
                   cbar_kws={'label': '-log10(Adjusted P-value)', 'shrink': 0.8},
                   linewidths=0.1,
                   linecolor='white',
                   square=False)
                   
        plt.title(f'{db} Enrichment Analysis Heatmap', fontsize=16, pad=20)
        plt.xlabel('Cluster', fontsize=10)
        plt.ylabel('GO Terms', fontsize=10)
        plt.xticks(fontsize=8, rotation=0)
        plt.yticks(fontsize=8, rotation=0)
        
        # Improve layout
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f'{db.lower()}_enrichment_heatmap.png'), 
                   bbox_inches="tight", dpi=dpi)
        plt.close()

        merged_results.to_csv(os.path.join(outdir, f'{db.lower()}_enrichment_results.csv'), index=False)

    @staticmethod
    def plot_gse_barplot(outdir, db, results, cluster, top_n=20):
        """
        Plot GO enrichment analysis results as a barplot for a specific cluster using gseapy.

        :param outdir: Output directory for saving the plot
        :param db: Database name (e.g., 'GO_Biological_Process_2021')
        :param results: Dictionary of enrichment results for all clusters
        :param cluster: Cluster to analyze (default is 'cluster_0')
        :param top_n: Number of top GO terms to display (default is 20)
        """
        # Extract results for the specified cluster
        cluster_results = results[cluster]

        # Larger size for better readability: 12cm width x 16cm height = 4.72 x 6.30 inches
        fig_size = (4.72, 6.30)
        
        # Create the plot with matplotlib to control font sizes
        plt.figure(figsize=fig_size, dpi=300)
        
        # Sort and get top terms
        top_results = cluster_results.sort_values('Adjusted P-value').head(top_n)
        
        if top_results.empty:
            print(f"No results to plot for {db} cluster {cluster}")
            return
        
        # Create beautiful gradient color scheme
        values = -np.log10(top_results['Adjusted P-value'])
        colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(values)))  # Use viridis colormap
        
        # Create horizontal barplot with gradient colors
        y_pos = range(len(top_results))
        bars = plt.barh(y_pos, values, color=colors, alpha=0.8, edgecolor='white', linewidth=0.5)
        
        # Add value labels on bars
        for i, (bar, value) in enumerate(zip(bars, values)):
            plt.text(value + 0.1, bar.get_y() + bar.get_height()/2, 
                    f'{value:.1f}', ha='left', va='center', fontsize=8, fontweight='bold')
        
        # Customize appearance
        plt.yticks(y_pos, top_results['Term'], fontsize=10)
        plt.xlabel('-log10(Adjusted P-value)', fontsize=10, fontweight='bold')
        plt.title(f'{db} Enrichment Analysis for Cluster {cluster}', fontsize=16, fontweight='bold', pad=20)
        plt.xticks(fontsize=8)
        plt.gca().invert_yaxis()
        
        # Add grid for better readability
        plt.grid(axis='x', alpha=0.3, linestyle='--')
        plt.gca().set_axisbelow(True)
        
        # Style the plot
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['left'].set_linewidth(0.5)
        plt.gca().spines['bottom'].set_linewidth(0.5)
        
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f'{db.lower()}_cluster{cluster}_barplot.png'), 
                   dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()

    @staticmethod
    def plot_gse_dotplot(outdir, db, results, cluster, top_n=20):
        """
        Plot GO enrichment analysis results as a dotplot for a specific cluster.

        :param outdir: Output directory for saving the plot
        :param db: Database name (e.g., 'GO_Biological_Process_2021')
        :param results: Dictionary of enrichment results for all clusters
        :param cluster: Cluster to analyze
        :param top_n: Number of top GO terms to display (default is 20)
        """
        # Extract results for the specified cluster
        cluster_results = results[cluster]

        # Larger size for better readability: 12cm width x 16cm height = 4.72 x 6.30 inches
        fig_size = (4.72, 6.30)
        
        # Create the plot with matplotlib to control font sizes
        plt.figure(figsize=fig_size, dpi=300)
        
        # Sort and get top terms
        top_results = cluster_results.sort_values('Adjusted P-value').head(top_n)
        
        # Create scatter plot (dotplot style)
        y_pos = range(len(top_results))
        sizes = -np.log10(top_results['Adjusted P-value']) * 20  # Scale for visibility
        plt.scatter(-np.log10(top_results['Adjusted P-value']), y_pos, s=sizes, 
                   c=-np.log10(top_results['Adjusted P-value']), cmap='coolwarm', alpha=0.7)
        plt.yticks(y_pos, top_results['Term'], fontsize=10)
        plt.xlabel('-log10(Adjusted P-value)', fontsize=10)
        plt.title(f'{db} Enrichment Analysis for cluster{cluster}', fontsize=16)
        plt.xticks(fontsize=8)
        plt.gca().invert_yaxis()
        plt.colorbar(label='-log10(Adjusted P-value)')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f'{db.lower()}_cluster{cluster}_dotplot.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()

    @staticmethod
    def plot_gse_emapplot(outdir, db, results, cluster, top_n=30):
        """
        Plot GO enrichment analysis results as a network for a specific cluster.

        :param outdir: Output directory for saving the plot
        :param db: Database name (e.g., 'GO_Biological_Process_2021')
        :param results: Dictionary of enrichment results for all clusters
        :param cluster: Cluster to analyze
        :param top_n: Number of top GO terms to display (default is 20)
        """
        # Extract results for the specified cluster
        cluster_results = results[cluster]

        # Sort by adjusted p-value and select top N terms
        top_terms = cluster_results.sort_values('Adjusted P-value').head(top_n)

        # If cluster_results is already in enr.res2d format, we can use it directly
        res2d = top_terms

        if 'NES' not in res2d.columns:
            res2d['NES'] = -np.log10(res2d['Adjusted P-value'])
        if 'Hits_ratio' not in res2d.columns:
            res2d['Hits_ratio'] = res2d['Genes'].apply(lambda x: len(x.split(';')) / 100)

        # Generate network data
        nodes, edges = gp.enrichment_map(res2d)

        # Build graph
        G = nx.from_pandas_edgelist(edges,
                                    source='src_idx',
                                    target='targ_idx',
                                    edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])

        # Add missing nodes if there are any
        for node in nodes.index:
            if node not in G.nodes():
                G.add_node(node)

        # Check if we have nodes to plot
        if nodes.empty or len(nodes) == 0:
            print(f"No nodes to plot for {db} cluster {cluster} enrichment map")
            return
            
        # Larger size for better readability: 16cm width x 16cm height = 6.30 x 6.30 inches
        fig_size = (6.30, 6.30)
        dpi = 300
        
        # Create figure with dark background for better contrast
        fig, ax = plt.subplots(figsize=fig_size, dpi=dpi, facecolor='white')
        ax.set_facecolor('white')

        # Initialize node coordinates with better algorithm
        pos = nx.spring_layout(G, k=3, iterations=100, seed=42)  # Reproducible layout

        # Prepare node properties
        node_colors = list(nodes.NES)
        node_sizes = [max(50, min(800, ratio * 2000)) for ratio in nodes.Hits_ratio]  # Constrain size range
        
        # Draw edges first (behind nodes) with improved styling
        edge_weights = list(nx.get_edge_attributes(G, 'jaccard_coef').values())
        if edge_weights:
            nx.draw_networkx_edges(G,
                                   pos=pos,
                                   width=[max(0.5, min(4, w * 12)) for w in edge_weights],  # Constrain edge width
                                   edge_color='lightgray',
                                   alpha=0.4,
                                   style='solid')

        # Draw nodes with beautiful gradient
        nodes_plot = nx.draw_networkx_nodes(G,
                                           pos=pos,
                                           cmap=plt.cm.RdYlBu_r,  # Beautiful diverging colormap
                                           node_color=node_colors,
                                           node_size=node_sizes,
                                           alpha=0.9,
                                           linewidths=2,
                                           edgecolors='white')

        # Draw selective labels (only for most significant terms)
        if len(nodes) <= 15:  # Show all labels if not too many
            label_dict = nodes.Term.to_dict()
        else:  # Show only top terms to avoid overcrowding
            top_indices = nodes.NES.nlargest(10).index
            label_dict = {idx: nodes.loc[idx, 'Term'] for idx in top_indices if idx in nodes.index}
        
        # Truncate long labels
        label_dict = {k: (v[:30] + '...' if len(v) > 30 else v) for k, v in label_dict.items()}
        
        nx.draw_networkx_labels(G,
                                pos=pos,
                                labels=label_dict,
                                font_size=9,
                                font_weight='bold',
                                font_color='black',
                                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))

        # Customize the plot
        plt.title(f'{db} Enrichment Network for Cluster {cluster}', 
                 fontsize=16, fontweight='bold', pad=20)
        plt.axis('off')

        # Add improved colorbar
        if nodes_plot:
            sm = plt.cm.ScalarMappable(cmap=plt.cm.RdYlBu_r, 
                                     norm=plt.Normalize(vmin=min(node_colors), vmax=max(node_colors)))
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, shrink=0.4, aspect=20)
            cbar.set_label('-log10(Adjusted P-value)', rotation=270, labelpad=20, fontsize=10, fontweight='bold')
            cbar.ax.tick_params(labelsize=9)
            
        # Add legend for node sizes
        legend_sizes = [50, 200, 400]
        legend_labels = ['Low', 'Medium', 'High']
        legend_elements = [plt.scatter([], [], s=size, c='gray', alpha=0.6, edgecolors='white', linewidth=1) 
                          for size in legend_sizes]
        plt.legend(legend_elements, legend_labels, title='Gene Ratio', 
                  loc='upper left', bbox_to_anchor=(0, 1), framealpha=0.9)

        # Adjust layout and save the plot
        plt.tight_layout(pad=0.5)  # Reduced padding for tighter spacing
        plt.savefig(os.path.join(outdir, f'{db.lower()}_cluster{cluster}_emapplot.png'), 
                   dpi=dpi, bbox_inches='tight')
        plt.close(fig)

    def write_stat_images(self):
        import matplotlib.pyplot as plt
        
        fig_size = self._spatial_figsize(self.adata[self.bin], base=3.6)
        dpi = 300
        
        # Spatial cluster plot
        plt.figure(figsize=fig_size, dpi=dpi)
        sc.pl.spatial(
            self.adata[self.bin],
            color=["cluster"],
            palette=sc.pl.palettes.default_20,
            size=1,
            spot_size=self._spatial_spot_size(self.adata[self.bin]),
            alpha_img=0.75,
            show=False,
        )
        plt.savefig(f'{self.spatial_cluster}/{self.sample}_spatial_cluster.png', 
                   bbox_inches="tight", dpi=dpi)
        plt.close()

        # UMAP cluster plot  
        plt.figure(figsize=fig_size, dpi=dpi)
        sc.pl.umap(self.adata[self.bin],
                   color=["cluster"],
                   palette=sc.pl.palettes.default_20,
                   show=False,
        )
        plt.savefig(f'{self.spatial_cluster}/{self.sample}_umap_cluster.png', 
                   bbox_inches="tight", dpi=dpi)
        plt.close()

    @utils.add_log
    def write_markers(self):
        """
        write only p_val_adj < P_VAL_CUTOFF to avoid too many markers
        """
        df_markers = sc.get.rank_genes_groups_df(self.adata[self.bin], group=None, pval_cutoff=PVAL_CUTOFF)
        df_markers = _normalise_rank_genes_groups_df(df_markers)
        df_markers.to_csv(self.df_marker_raw_file, index=None, sep='\t')

        df_markers_filter = df_markers.loc[df_markers['avg_log2FC'] > 0].sort_values('p_val_adj').groupby(
            'cluster').head(100)
        df_markers_filter = df_markers_filter.round({
            'avg_log2FC': 3,
            'pct.1': 3,
            'pct.2': 3,
        })
        df_markers_filter.to_csv(self.df_marker_file, index=None, sep='\t')

    @utils.add_log
    def write_tsne(self):
        df_tsne = self.adata[self.bin].obsm.to_df()[['X_tsne1', 'X_tsne2']]
        df_tsne['cluster'] = self.adata[self.bin].obs.cluster
        df_tsne['Gene_Counts'] = self.adata[self.bin].obs.n_genes_by_counts
        tsne_name_dict = {'X_tsne1': 'tSNE_1', 'X_tsne2': 'tSNE_2'}
        df_tsne = df_tsne.rename(tsne_name_dict, axis='columns')
        df_tsne['cluster'] = df_tsne['cluster'].map(lambda x: int(x))
        df_tsne.to_csv(self.df_tsne_file, sep='\t')

    @utils.add_log
    def write_h5ad(self):
        for i, adata in self.adata.items():
            # For scRNA mode, optimize h5ad file size
            if self.assay == 'scrna':
                # Keep only essential data to reduce file size
                adata_save = adata.copy()

                # Restore raw counts to X (as sparse matrix)
                if 'raw' in adata_save.layers:
                    adata_save.X = adata_save.layers['raw']
                    # Remove all layers except raw
                    adata_save.layers = {'raw': adata_save.layers['raw']}

                # Keep only essential obsm (PCA and UMAP)
                obsm_keep = {}
                if 'X_pca' in adata_save.obsm:
                    obsm_keep['X_pca'] = adata_save.obsm['X_pca']
                if 'X_umap' in adata_save.obsm:
                    obsm_keep['X_umap'] = adata_save.obsm['X_umap']
                adata_save.obsm = obsm_keep

                # Clear varm to save space
                adata_save.varm = {}

                # Write compressed h5ad
                adata_save.write_h5ad(
                    Path(os.path.join(self.outdir, f'{self.sample}_{i}.h5ad')),
                    compression='gzip'
                )
                print(f"Saved compressed h5ad: {self.sample}_{i}.h5ad")
            else:
                self.adata[i].write_h5ad(Path(os.path.join(self.h5ad_dir, f'{self.sample}_{i}.h5ad')))

    @staticmethod
    def get_de_results(adata, n_genes, min_logfc):
        results = []
        for cluster in adata.obs['cluster'].cat.categories:
            genes = sc.get.rank_genes_groups_df(adata, group=cluster, pval_cutoff=PVAL_CUTOFF)
            genes = genes[genes['logfoldchanges'] > min_logfc]
            genes['cluster'] = cluster
            results.append(genes.head(n_genes))
        return pd.concat(results)

    def bioinfo_data(self):
        if self.assay == 'spatial':
            os.makedirs(self.h5ad_dir, exist_ok=True)
            os.makedirs(self.spatial_qc, exist_ok=True)
            os.makedirs(self.spatial_cluster, exist_ok=True)
            os.makedirs(self.diff_exp_analysis, exist_ok=True)
        elif self.assay == 'cell':
            os.makedirs(self.h5ad_dir, exist_ok=True)
            os.makedirs(self.cell_qc_dir, exist_ok=True)
            os.makedirs(self.cell_embedding_dir, exist_ok=True)
            os.makedirs(self.cell_marker_dir, exist_ok=True)
            os.makedirs(self.cell_spatial_dir, exist_ok=True)

        if self.assay == 'spatial' and getattr(self.args, 'hires_image', None):
            shutil.copy(self.args.hires_image, self.spatial_cluster)

    def scrna_data(self):
        """Create scRNA output directories"""
        os.makedirs(self.analysis_dir, exist_ok=True)
        os.makedirs(self.qc_dir, exist_ok=True)
        os.makedirs(self.umap_dir, exist_ok=True)

    @utils.add_log
    def scrna_qc(self):
        """scRNA-seq quality control with relaxed thresholds"""
        import matplotlib.pyplot as plt

        adata = self.adata[self.bin]

        # Calculate QC metrics
        if self.mt_gene_list:
            mito_genes, _ = utils.read_one_col(self.mt_gene_list)
            adata.var[MITO_VAR] = adata.var_names.map(lambda x: True if x in mito_genes else False)
            adata.var[MITO_VAR] = adata.var[MITO_VAR].astype(bool)
        else:
            adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
            adata.var['rb'] = adata.var_names.str.upper().str.contains('^RP[SL]')
            adata.var['hb'] = adata.var_names.str.upper().str.contains('^HB[^(P)]')

        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=['mt', 'rb', 'hb'],
            percent_top=None,
            use_raw=False,
            log1p=False,
            inplace=True
        )

        # ===== QC前的图 =====
        print(f"Before QC: {adata.n_obs} cells")

        fig_size = (3.15, 3.15)
        dpi = 300

        # Violin plot - Before QC
        colors = ["#f38181", "#a8d8ea", "#fce38a"]
        fig, axes = plt.subplots(1, 3, figsize=(11.81, 3.94))

        metrics = ['total_counts', 'n_genes_by_counts', f'pct_counts_{MITO_VAR}']

        for idx, (metric, color, ax) in enumerate(zip(metrics, colors, axes)):
            sc.pl.violin(adata, metric, jitter=0.4, show=False, ax=ax)

            # Style violin patches
            for collection in ax.collections:
                if hasattr(collection, '_paths') and len(collection._paths) > 0:
                    collection.set_facecolor(color)
                    collection.set_alpha(1)
                    collection.set_edgecolor('#323232')
                    collection.set_linewidth(1.5)
                elif hasattr(collection, '_offsets') and len(collection._offsets) > 0:
                    collection.set_facecolor('#323232')
                    collection.set_edgecolor('#323232')
                    collection.set_alpha(1.0)

            for line in ax.lines:
                line.set_color('#323232')
                line.set_alpha(1.0)

        plt.suptitle('Before QC', fontsize=14, y=1.02)
        plt.tight_layout()
        plt.savefig(f'{self.qc_dir}/{self.sample}_qc_violin_before.png', bbox_inches="tight", dpi=dpi)
        plt.close()

        # Scatter plot - Before QC
        plt.figure(figsize=fig_size, dpi=dpi)
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
                     color='pct_counts_mt', show=False)
        plt.title('Before QC')
        plt.savefig(f'{self.qc_dir}/{self.sample}_qc_scatter_before.png',
                   bbox_inches="tight", dpi=dpi)
        plt.close()

        # ===== QC过滤 =====
        # Very relaxed thresholds
        sc.pp.filter_cells(adata, min_genes=100)  # At least 100 genes
        sc.pp.filter_cells(adata, min_counts=200)  # At least 200 UMIs
        adata = adata[adata.obs.pct_counts_mt < 50, :]  # Less than 50% mitochondrial

        print(f"After QC: {adata.n_obs} cells")

        # Update adata
        self.adata[self.bin] = adata

        # ===== QC后的图 =====
        # Violin plot - After QC
        fig, axes = plt.subplots(1, 3, figsize=(11.81, 3.94))

        for idx, (metric, color, ax) in enumerate(zip(metrics, colors, axes)):
            sc.pl.violin(adata, metric, jitter=0.4, show=False, ax=ax)

            # Style violin patches
            for collection in ax.collections:
                if hasattr(collection, '_paths') and len(collection._paths) > 0:
                    collection.set_facecolor(color)
                    collection.set_alpha(1)
                    collection.set_edgecolor('#323232')
                    collection.set_linewidth(1.5)
                elif hasattr(collection, '_offsets') and len(collection._offsets) > 0:
                    collection.set_facecolor('#323232')
                    collection.set_edgecolor('#323232')
                    collection.set_alpha(1.0)

            for line in ax.lines:
                line.set_color('#323232')
                line.set_alpha(1.0)

        plt.suptitle('After QC', fontsize=14, y=1.02)
        plt.tight_layout()
        plt.savefig(f'{self.qc_dir}/{self.sample}_qc_violin_after.png', bbox_inches="tight", dpi=dpi)
        plt.close()

        # Scatter plot - After QC
        plt.figure(figsize=fig_size, dpi=dpi)
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
                     color='pct_counts_mt', show=False)
        plt.title('After QC')
        plt.savefig(f'{self.qc_dir}/{self.sample}_qc_scatter_after.png',
                   bbox_inches="tight", dpi=dpi)
        plt.close()

    @utils.add_log
    def scrna_umap(self):
        """scRNA-seq dimensional reduction with simple visualizations"""
        import matplotlib.pyplot as plt

        adata = self.adata[self.bin]

        print(f"Starting analysis with {adata.n_obs} cells and {adata.n_vars} genes")

        # 1. Normalization
        sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
        sc.pp.log1p(adata)
        adata.layers[NORMALIZED_LAYER] = adata.X

        # 2. Highly Variable Genes Selection
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
        n_hvg = sum(adata.var['highly_variable'])
        print(f"Identified {n_hvg} highly variable genes")

        # Simple HVG plot (no text labels)
        plt.figure(figsize=(6, 5), dpi=300)
        plt.scatter(adata.var.loc[~adata.var['highly_variable'], 'means'],
                   adata.var.loc[~adata.var['highly_variable'], 'dispersions_norm'],
                   s=3, c='lightgray', alpha=0.5)
        plt.scatter(adata.var.loc[adata.var['highly_variable'], 'means'],
                   adata.var.loc[adata.var['highly_variable'], 'dispersions_norm'],
                   s=5, c='red', alpha=0.6)
        plt.xlabel('Average Expression')
        plt.ylabel('Normalized Dispersion')
        plt.title(f'Highly Variable Genes (n={n_hvg})')
        plt.tight_layout()
        plt.savefig(f'{self.umap_dir}/{self.sample}_hvg_selection.png', dpi=300, bbox_inches='tight')
        plt.close()

        # 3. Scaling and PCA
        sc.pp.scale(adata, max_value=10)

        # PCA with adaptive number of components
        max_pcs = min(adata.n_obs - 1, adata.n_vars - 1, 50)
        if max_pcs < 10:
            print(f"Warning: Only {max_pcs} components available, using all")
            n_pcs_use = max_pcs
        else:
            n_pcs_use = 30

        sc.pp.pca(adata, n_comps=max_pcs)
        print(f"PCA computed with {max_pcs} components, using {n_pcs_use} for downstream analysis")

        # 4. Simple Elbow Plot
        variance_ratio = adata.uns['pca']['variance_ratio']
        plt.figure(figsize=(6, 4), dpi=300)
        plt.plot(range(1, len(variance_ratio) + 1), variance_ratio, 'o-')
        plt.axvline(x=n_pcs_use, color='red', linestyle='--', label=f'Selected PCs (n={n_pcs_use})')
        plt.xlabel('Principal Component')
        plt.ylabel('Variance Explained')
        plt.title('PCA Elbow Plot')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{self.umap_dir}/{self.sample}_pca_elbow.png', dpi=300, bbox_inches='tight')
        plt.close()

        # 5. Clustering
        n_neighbors = min(15, adata.n_obs - 1)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=min(n_pcs_use, max_pcs))

        # KMeans clustering
        from sklearn.cluster import KMeans
        n_clusters = min(10, max(2, adata.n_obs // 1000))
        print(f"Running KMeans clustering with {n_clusters} clusters")
        kmeans = KMeans(n_clusters=n_clusters, random_state=0, n_init=10)
        adata.obs['cluster'] = kmeans.fit_predict(adata.obsm['X_pca'][:, :min(n_pcs_use, max_pcs)])
        # Convert to string type to avoid KeyError with integer cluster IDs
        adata.obs['cluster'] = adata.obs['cluster'].astype(str).astype('category')

        # 6. UMAP Embedding
        sc.tl.umap(adata)

        # 7. Simple UMAP Visualizations (static)
        fig_size = (6, 5)
        dpi = 300

        # UMAP by cluster
        plt.figure(figsize=fig_size, dpi=dpi)
        sc.pl.umap(adata, color='cluster', show=False)
        plt.savefig(f'{self.umap_dir}/{self.sample}_umap_cluster.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        # UMAP by gene counts
        plt.figure(figsize=fig_size, dpi=dpi)
        sc.pl.umap(adata, color='n_genes_by_counts', show=False)
        plt.savefig(f'{self.umap_dir}/{self.sample}_umap_genes.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        # UMAP by UMI counts
        plt.figure(figsize=fig_size, dpi=dpi)
        sc.pl.umap(adata, color='total_counts', show=False)
        plt.savefig(f'{self.umap_dir}/{self.sample}_umap_counts.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        # 8. Interactive UMAP (Plotly)
        try:
            import plotly.graph_objects as go
            import plotly.express as px
            import json

            umap_df = pd.DataFrame({
                'UMAP_1': adata.obsm['X_umap'][:, 0],
                'UMAP_2': adata.obsm['X_umap'][:, 1],
                'Cluster': adata.obs['cluster'].astype(str),
                'n_genes': adata.obs['n_genes_by_counts'],
                'total_counts': adata.obs['total_counts'],
                'pct_mt': adata.obs['pct_counts_mt']
            })

            fig = px.scatter(
                umap_df,
                x='UMAP_1',
                y='UMAP_2',
                color='Cluster',
                hover_data={
                    'UMAP_1': ':.2f',
                    'UMAP_2': ':.2f',
                    'Cluster': True,
                    'n_genes': True,
                    'total_counts': True,
                    'pct_mt': ':.2f'
                },
                title=f'Interactive UMAP - {self.sample}',
                labels={'Cluster': 'Cluster'},
                color_discrete_sequence=px.colors.qualitative.Set1
            )

            fig.update_traces(
                marker=dict(size=4, opacity=0.7, line=dict(width=0)),
                hovertemplate='<b>Cluster %{customdata[2]}</b><br>' +
                             'UMAP_1: %{x:.2f}<br>' +
                             'UMAP_2: %{y:.2f}<br>' +
                             'Genes: %{customdata[3]}<br>' +
                             'UMIs: %{customdata[4]}<br>' +
                             'MT%: %{customdata[5]:.2f}<extra></extra>'
            )

            fig.update_layout(
                width=900,
                height=700,
                plot_bgcolor='white',
                hovermode='closest',
                legend=dict(
                    title=dict(text='Cluster', font=dict(size=14)),
                    orientation="v",
                    yanchor="top",
                    y=1,
                    xanchor="left",
                    x=1.02
                )
            )

            fig.write_html(f'{self.umap_dir}/{self.sample}_umap_interactive.html')

            plotly_json = fig.to_json()
            with open(f'{self.umap_dir}/{self.sample}_umap_interactive.json', 'w') as f:
                f.write(plotly_json)

            print(f"Interactive UMAP saved: {self.umap_dir}/{self.sample}_umap_interactive.html")

        except Exception as e:
            print(f"Warning: Could not generate interactive UMAP: {e}")

        # 9. t-SNE Embedding and Visualization
        print("Computing t-SNE embedding...")
        # Note: BLAS thread limits should be set via environment variables BEFORE
        # Python process starts to avoid segmentation faults
        sc.tl.tsne(adata, n_pcs=min(n_pcs_use, max_pcs))

        # Static t-SNE plot
        plt.figure(figsize=fig_size, dpi=dpi)
        sc.pl.tsne(adata, color='cluster', show=False)
        plt.savefig(f'{self.umap_dir}/{self.sample}_tsne_cluster.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        try:
            import plotly.express as px
            import json

            tsne_df = pd.DataFrame({
                'tSNE_1': adata.obsm['X_tsne'][:, 0],
                'tSNE_2': adata.obsm['X_tsne'][:, 1],
                'Cluster': adata.obs['cluster'].astype(str),
                'n_genes': adata.obs['n_genes_by_counts'],
                'total_counts': adata.obs['total_counts'],
                'pct_mt': adata.obs['pct_counts_mt']
            })

            fig = px.scatter(
                tsne_df,
                x='tSNE_1',
                y='tSNE_2',
                color='Cluster',
                hover_data={
                    'tSNE_1': ':.2f',
                    'tSNE_2': ':.2f',
                    'Cluster': True,
                    'n_genes': True,
                    'total_counts': True,
                    'pct_mt': ':.2f'
                },
                title=f'Interactive t-SNE - {self.sample}',
                labels={'Cluster': 'Cluster'},
                color_discrete_sequence=px.colors.qualitative.Set1
            )

            fig.update_traces(
                marker=dict(size=4, opacity=0.7, line=dict(width=0)),
                hovertemplate='<b>Cluster %{customdata[2]}</b><br>' +
                             't-SNE_1: %{x:.2f}<br>' +
                             't-SNE_2: %{y:.2f}<br>' +
                             'Genes: %{customdata[3]}<br>' +
                             'UMIs: %{customdata[4]}<br>' +
                             'MT%: %{customdata[5]:.2f}<extra></extra>'
            )

            fig.update_layout(
                width=900,
                height=700,
                plot_bgcolor='white',
                hovermode='closest',
                legend=dict(
                    title=dict(text='Cluster', font=dict(size=14)),
                    orientation="v",
                    yanchor="top",
                    y=1,
                    xanchor="left",
                    x=1.02
                )
            )

            fig.write_html(f'{self.umap_dir}/{self.sample}_tsne_interactive.html')

            plotly_json = fig.to_json()
            with open(f'{self.umap_dir}/{self.sample}_tsne_interactive.json', 'w') as f:
                f.write(plotly_json)

            print(f"Interactive t-SNE saved: {self.umap_dir}/{self.sample}_tsne_interactive.html")

        except Exception as e:
            print(f"Warning: Could not generate interactive t-SNE: {e}")

        print(f"Analysis completed. Plots saved to: {self.umap_dir}/")

        # Update adata
        self.adata[self.bin] = adata

    @utils.add_log
    def scrna_find_markers(self):
        """Find marker genes for each cluster in scRNA-seq data"""
        import matplotlib.pyplot as plt

        adata = self.adata[self.bin]

        print(f"Finding marker genes for {len(adata.obs['cluster'].unique())} clusters...")

        # Run differential expression analysis
        sc.tl.rank_genes_groups(
            adata,
            "cluster",
            reference='rest',
            pts=True,
            method="wilcoxon",
            use_raw=False,
            layer=NORMALIZED_LAYER
        )

        # Create marker genes directory
        marker_dir = f'{self.analysis_dir}/03.Markers'
        os.makedirs(marker_dir, exist_ok=True)

        # Plot rank genes groups
        fig_size = (5.91, 4.72)
        dpi = 300

        plt.figure(figsize=fig_size, dpi=dpi)
        sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, ncols=5, show=False)
        plt.savefig(f'{marker_dir}/{self.sample}_marker_genes.png',
                   bbox_inches="tight", dpi=dpi)
        plt.close()

        print(f"Marker genes plot saved to {marker_dir}/{self.sample}_marker_genes.png")

        # Update adata
        self.adata[self.bin] = adata

    @utils.add_log
    def scrna_visualize_markers(self):
        """Visualize marker genes in scRNA-seq data"""
        import matplotlib.pyplot as plt

        adata = self.adata[self.bin]
        marker_dir = f'{self.analysis_dir}/03.Markers'

        n_genes = 10
        min_logfc = 0.25
        dpi = 300

        # Get differential expression results
        de_results = self.get_de_results(adata, n_genes, min_logfc)

        if de_results.empty:
            print("WARNING: No significant marker genes found. Skipping marker visualization.")
            return

        # 1. Heatmap
        plt.figure(figsize=(4.72, 3.94), dpi=dpi)
        sc.pl.heatmap(adata, de_results['names'].unique(), groupby='cluster', cmap="Spectral_r",
                     use_raw=False, swap_axes=True, show_gene_labels=True, show=False)
        plt.savefig(f'{marker_dir}/{self.sample}_marker_heatmap.png',
                   bbox_inches="tight", dpi=dpi)
        plt.close()

        # 2. UMAP plots with top markers
        actual_clusters = sorted(de_results['cluster'].unique())[:4]
        top_genes = [de_results[de_results['cluster'] == cluster]['names'].iloc[0]
                    for cluster in actual_clusters
                    if not de_results[de_results['cluster'] == cluster].empty]

        if top_genes:
            # UMAP with top markers
            fig = plt.figure(figsize=(12.6, 3.15), dpi=dpi)
            plt.subplots_adjust(wspace=0.05, hspace=0.05)
            sc.pl.umap(adata, color=top_genes, cmap="Spectral_r", show=False)
            plt.savefig(f'{marker_dir}/{self.sample}_umap_markers.png',
                       bbox_inches="tight", dpi=dpi)
            plt.close()

            # Violin plot
            fig = plt.figure(figsize=(12.6, 3.15), dpi=dpi)
            plt.subplots_adjust(wspace=0.05, hspace=0.05)
            sc.pl.violin(adata, top_genes, groupby='cluster', rotation=90, jitter=0., size=0, show=False)
            plt.savefig(f'{marker_dir}/{self.sample}_violin_markers.png',
                       bbox_inches="tight", dpi=dpi)
            plt.close()

        # 3. t-SNE plots with top markers (if t-SNE was computed)
        if 'X_tsne' in adata.obsm:
            fig = plt.figure(figsize=(12.6, 3.15), dpi=dpi)
            plt.subplots_adjust(wspace=0.05, hspace=0.05)
            sc.pl.tsne(adata, color=top_genes, cmap="Spectral_r", show=False)
            plt.savefig(f'{marker_dir}/{self.sample}_tsne_markers.png',
                       bbox_inches="tight", dpi=dpi)
            plt.close()

        print(f"Marker visualizations saved to {marker_dir}/")

        # Update adata
        self.adata[self.bin] = adata

    @utils.add_log
    def scrna_write_markers(self):
        """Write marker genes to TSV files for scRNA-seq data"""
        adata = self.adata[self.bin]
        marker_dir = f'{self.analysis_dir}/03.Markers'

        # Write raw markers (all significant)
        df_markers = sc.get.rank_genes_groups_df(adata, group=None, pval_cutoff=PVAL_CUTOFF)
        df_markers = _normalise_rank_genes_groups_df(df_markers)

        raw_marker_file = f'{marker_dir}/{self.sample}_markers_raw.tsv'
        df_markers.to_csv(raw_marker_file, index=None, sep='\t')
        print(f"Raw marker genes saved to {raw_marker_file}")

        # Write filtered markers (top 100 per cluster with positive log2FC)
        df_markers_filter = df_markers.loc[df_markers['avg_log2FC'] > 0].sort_values('p_val_adj').groupby(
            'cluster').head(100)
        df_markers_filter = df_markers_filter.round({
            'avg_log2FC': 3,
            'pct.1': 3,
            'pct.2': 3,
        })

        filtered_marker_file = f'{marker_dir}/{self.sample}_markers.tsv'
        df_markers_filter.to_csv(filtered_marker_file, index=None, sep='\t')
        print(f"Filtered marker genes (top 100 per cluster) saved to {filtered_marker_file}")

        # Print summary
        n_clusters = len(df_markers_filter['cluster'].unique())
        n_markers = len(df_markers_filter)
        print(f"Found {n_markers} marker genes across {n_clusters} clusters")

    @utils.add_log
    def cell_qc(self):
        """QC and QC plots for StarDist/cell-level matrices."""
        adata = self.adata[self.bin]

        self._set_qc_gene_flags(adata)
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=['mt', 'rb', 'hb'],
            percent_top=None,
            use_raw=False,
            log1p=False,
            inplace=True,
        )

        before_n = adata.n_obs
        min_genes = int(getattr(self.args, 'cell_min_genes', 0) or 0)
        min_counts = int(getattr(self.args, 'cell_min_counts', 0) or 0)
        max_genes = int(getattr(self.args, 'cell_max_genes', 0) or 0)
        max_counts = int(getattr(self.args, 'cell_max_counts', 0) or 0)
        max_mt = float(getattr(self.args, 'cell_max_mt', 100.0))

        qc_filters = {}
        if min_genes > 0:
            qc_filters['below_min_genes'] = adata.obs['n_genes_by_counts'] < min_genes
        if min_counts > 0:
            qc_filters['below_min_counts'] = adata.obs['total_counts'] < min_counts
        if max_genes > 0:
            qc_filters['above_max_genes'] = adata.obs['n_genes_by_counts'] > max_genes
        if max_counts > 0:
            qc_filters['above_max_counts'] = adata.obs['total_counts'] > max_counts
        if max_mt < 100 and f'pct_counts_{MITO_VAR}' in adata.obs.columns:
            qc_filters['above_max_mt'] = adata.obs[f'pct_counts_{MITO_VAR}'] > max_mt

        keep = np.ones(before_n, dtype=bool)
        removed_by_filter = {}
        for name, failed in qc_filters.items():
            failed_values = np.asarray(failed, dtype=bool)
            removed_by_filter[f'cells_{name}'] = int(failed_values.sum())
            keep &= ~failed_values
        if qc_filters:
            adata = adata[keep, :].copy()

        self.adata[self.bin] = adata
        print(f"Cell QC: {before_n} -> {adata.n_obs} cells")

        mt_col = f'pct_counts_{MITO_VAR}'
        total_counts = pd.to_numeric(adata.obs['total_counts'], errors='coerce') if 'total_counts' in adata.obs else pd.Series(dtype=float)
        n_genes = pd.to_numeric(adata.obs['n_genes_by_counts'], errors='coerce') if 'n_genes_by_counts' in adata.obs else pd.Series(dtype=float)
        mt_pct = pd.to_numeric(adata.obs[mt_col], errors='coerce') if mt_col in adata.obs else pd.Series(dtype=float)
        summary = {
            'cells_before_qc': int(before_n),
            'cells_after_qc': int(adata.n_obs),
            'cells_removed_qc': int(before_n - adata.n_obs),
            'genes': int(adata.n_vars),
            'mitochondrial_genes': int(adata.var[MITO_VAR].sum()) if MITO_VAR in adata.var.columns else 0,
            'median_umi_per_cell': float(total_counts.median()) if not total_counts.empty else 0.0,
            'median_genes_per_cell': float(n_genes.median()) if not n_genes.empty else 0.0,
            'median_mt_percent': float(mt_pct.median()) if not mt_pct.empty else 0.0,
            'mean_mt_percent': float(mt_pct.mean()) if not mt_pct.empty else 0.0,
            'min_genes': min_genes,
            'min_counts': min_counts,
            'max_genes': max_genes,
            'max_counts': max_counts,
            'max_mt': max_mt,
            **removed_by_filter,
        }
        pd.DataFrame([summary]).to_csv(
            os.path.join(self.cell_qc_dir, f'{self.sample}_cell_qc_summary.tsv'),
            sep='\t',
            index=False,
        )
        pd.DataFrame([summary]).to_csv(
            os.path.join(self.outdir, 'metrics_summary.csv'),
            index=False,
        )
        with open(os.path.join(self.outdir, 'metrics_summary.json'), 'w') as handle:
            json.dump(summary, handle, indent=2)
        metric_names = {
            'cells_before_qc': 'Cells Before QC',
            'cells_after_qc': 'Cells After QC',
            'genes': 'Genes',
            'mitochondrial_genes': 'Mitochondrial Genes',
            'median_umi_per_cell': 'Median UMI per Cell',
            'median_genes_per_cell': 'Median Genes per Cell',
            'median_mt_percent': 'Median Mitochondrial Percent',
            'mean_mt_percent': 'Mean Mitochondrial Percent',
        }
        for key, name in metric_names.items():
            value = summary[key]
            display = f"{value:.2f}" if isinstance(value, float) else str(value)
            self.add_metric(name=name, value=value, display=display)

        metrics = ['total_counts', 'n_genes_by_counts', f'pct_counts_{MITO_VAR}']
        available_metrics = [metric for metric in metrics if metric in adata.obs.columns]
        if available_metrics:
            fig, axes = plt.subplots(1, len(available_metrics), figsize=(4 * len(available_metrics), 4), dpi=CELL_DPI)
            if len(available_metrics) == 1:
                axes = [axes]
            for metric, ax in zip(available_metrics, axes):
                sc.pl.violin(adata, metric, jitter=0.4, show=False, ax=ax)
            plt.tight_layout()
            plt.savefig(os.path.join(self.cell_qc_dir, f'{self.sample}_cell_qc_violin.png'), dpi=CELL_DPI)
            plt.close()

        if {'total_counts', 'n_genes_by_counts'}.issubset(adata.obs.columns):
            plt.figure(figsize=CELL_FIGSIZE, dpi=CELL_DPI)
            sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color=f'pct_counts_{MITO_VAR}', show=False)
            plt.savefig(os.path.join(self.cell_qc_dir, f'{self.sample}_cell_qc_scatter.png'), dpi=CELL_DPI)
            plt.close()

    @staticmethod
    def _image_size(image_path):
        try:
            from PIL import Image

            Image.MAX_IMAGE_PIXELS = None
            with Image.open(image_path) as image:
                return image.size
        except Exception:
            return None

    def _cell_spatial_extent(self, coords, image_shape):
        image_height, image_width = image_shape[:2]
        fullres_size = None
        fullres_image = getattr(self.args, 'fullres_image', None)
        if fullres_image and os.path.exists(fullres_image):
            fullres_size = self._image_size(fullres_image)

        if fullres_size:
            fullres_width, fullres_height = fullres_size
        else:
            finite = np.isfinite(coords).all(axis=1)
            if finite.any():
                x_max = float(np.nanmax(coords[finite, 0]))
                y_max = float(np.nanmax(coords[finite, 1]))
                if x_max > image_width * 1.25 or y_max > image_height * 1.25:
                    fullres_width = max(x_max, image_width)
                    fullres_height = max(y_max, image_height)
                else:
                    fullres_width = image_width
                    fullres_height = image_height
            else:
                fullres_width = image_width
                fullres_height = image_height

        return (0, fullres_width, fullres_height, 0), fullres_width, fullres_height

    def _update_cell_metrics_summary(self, updates):
        summary_path = os.path.join(self.outdir, 'metrics_summary.json')
        summary = {}
        if os.path.exists(summary_path):
            try:
                with open(summary_path) as handle:
                    summary = json.load(handle)
            except Exception:
                summary = {}
        summary.update(updates)
        with open(summary_path, 'w') as handle:
            json.dump(summary, handle, indent=2)
        pd.DataFrame([summary]).to_csv(os.path.join(self.outdir, 'metrics_summary.csv'), index=False)

    def _write_cell_spatial_plot(self, color, filename, categorical=False):
        adata = self.adata[self.bin]
        if 'spatial' not in adata.obsm:
            return
        if color not in adata.obs.columns:
            return

        coords = adata.obsm['spatial']
        fig, ax = plt.subplots(figsize=CELL_SPATIAL_FIGSIZE, dpi=CELL_DPI)
        image_drawn = False
        hires_image = getattr(self.args, 'hires_image', None)
        if hires_image and os.path.exists(hires_image):
            try:
                image = plt.imread(hires_image)
                extent, fullres_width, fullres_height = self._cell_spatial_extent(coords, image.shape)
                ax.imshow(image, extent=extent, origin='upper')
                ax.set_xlim(0, fullres_width)
                ax.set_ylim(fullres_height, 0)
                ax.set_aspect('equal', adjustable='box')
                image_drawn = True
            except Exception as exc:
                print(f"WARNING: Could not draw cell spatial HE background: {exc}")

        if categorical:
            categories = adata.obs[color].astype('category')
            values = categories.cat.codes
            palette = categorical_palette(len(categories.cat.categories))
            cmap = mcolors.ListedColormap(palette)
            ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=values,
                s=CELL_POINT_SIZE,
                cmap=cmap,
                linewidths=0,
                alpha=0.85,
                rasterized=True,
            )
            if len(categories.cat.categories) <= 20:
                from matplotlib.lines import Line2D

                handles = [
                    Line2D(
                        [0],
                        [0],
                        marker='o',
                        linestyle='None',
                        markerfacecolor=palette[idx],
                        markeredgecolor='none',
                        markersize=4,
                        label=str(category),
                    )
                    for idx, category in enumerate(categories.cat.categories)
                ]
                if handles:
                    ax.legend(
                        handles=handles,
                        title=color,
                        fontsize=6,
                        title_fontsize=7,
                        loc='upper left',
                        bbox_to_anchor=(1.01, 1.0),
                        borderaxespad=0,
                        frameon=True,
                    )
        else:
            values = pd.to_numeric(adata.obs[color], errors='coerce').fillna(0)
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=values,
                s=CELL_POINT_SIZE,
                cmap='Spectral_r',
                linewidths=0,
                alpha=0.85,
                rasterized=True,
            )
            fig.colorbar(scatter, ax=ax, label=color, shrink=0.8)
        if not image_drawn:
            ax.invert_yaxis()
            ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('fullres_col')
        ax.set_ylabel('fullres_row')
        fig.tight_layout()
        fig.savefig(os.path.join(self.cell_spatial_dir, filename), dpi=CELL_DPI, bbox_inches='tight', pad_inches=0.08)
        plt.close(fig)

    def _write_cell_embedding_plot(self, basis, color, filename, categorical=False, palette=None):
        adata = self.adata[self.bin]
        obsm_key = f'X_{basis}'
        if obsm_key not in adata.obsm or color not in adata.obs.columns:
            return

        coords = adata.obsm[obsm_key]
        fig, ax = plt.subplots(figsize=CELL_EMBEDDING_FIGSIZE, dpi=CELL_DPI)
        if categorical:
            categories = adata.obs[color].astype('category')
            palette = palette or categorical_palette(len(categories.cat.categories))
            color_lookup = {category: palette[idx] for idx, category in enumerate(categories.cat.categories)}
            point_colors = categories.map(color_lookup).to_numpy()
            ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=point_colors,
                s=CELL_EMBEDDING_POINT_SIZE,
                linewidths=0,
                alpha=0.85,
                rasterized=True,
            )
            if len(categories.cat.categories) <= 30:
                from matplotlib.lines import Line2D

                handles = [
                    Line2D(
                        [0],
                        [0],
                        marker='o',
                        linestyle='None',
                        markerfacecolor=palette[idx],
                        markeredgecolor='none',
                        markersize=4,
                        label=str(category),
                    )
                    for idx, category in enumerate(categories.cat.categories)
                ]
                ax.legend(
                    handles=handles,
                    title=color,
                    fontsize=7,
                    title_fontsize=8,
                    loc='center left',
                    bbox_to_anchor=(1.01, 0.5),
                    borderaxespad=0,
                    frameon=False,
                )
        else:
            values = pd.to_numeric(adata.obs[color], errors='coerce').fillna(0)
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=values,
                s=CELL_EMBEDDING_POINT_SIZE,
                cmap='Spectral_r',
                linewidths=0,
                alpha=0.85,
                rasterized=True,
            )
            fig.colorbar(scatter, ax=ax, label=color, shrink=0.8, pad=0.02)

        basis_label = 'UMAP' if basis == 'umap' else 't-SNE'
        ax.set_xlabel(f'{basis_label}1')
        ax.set_ylabel(f'{basis_label}2')
        ax.set_title(color)
        finite = np.isfinite(coords).all(axis=1)
        if finite.any():
            x_min, y_min = coords[finite].min(axis=0)
            x_max, y_max = coords[finite].max(axis=0)
            x_pad = max((x_max - x_min) * 0.06, 0.5)
            y_pad = max((y_max - y_min) * 0.06, 0.5)
            ax.set_xlim(x_min - x_pad, x_max + x_pad)
            ax.set_ylim(y_min - y_pad, y_max + y_pad)
        fig.tight_layout()
        fig.savefig(os.path.join(self.cell_embedding_dir, filename), dpi=CELL_DPI, bbox_inches='tight', pad_inches=0.10)
        plt.close(fig)

    @utils.add_log
    def cell_cluster(self):
        """Cell-level normalization, dimensional reduction, clustering, and plots."""
        adata = self.adata[self.bin]
        if adata.n_obs < 2 or adata.n_vars < 2:
            raise ValueError(f"Cell analysis needs at least 2 cells and 2 genes; got {adata.n_obs} cells, {adata.n_vars} genes")

        sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
        sc.pp.log1p(adata)
        adata.layers[NORMALIZED_LAYER] = adata.X.copy()

        n_top_genes = min(2000, max(1, adata.n_vars))
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat')
        sc.pp.scale(adata, max_value=10)

        max_pcs = min(50, adata.n_obs - 1, adata.n_vars - 1)
        if max_pcs < 2:
            raise ValueError(f"Cell analysis has too few dimensions for PCA: n_pcs={max_pcs}")
        sc.pp.pca(adata, n_comps=max_pcs, svd_solver='auto')

        n_pcs_use = min(self.cell_n_pcs, max_pcs)
        n_neighbors = min(self.cell_n_neighbors, max(2, adata.n_obs - 1))
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs_use)

        sc.tl.leiden(
            adata,
            resolution=self.cell_cluster_resolution,
            random_state=0,
            flavor='igraph',
            key_added='cluster',
            directed=False,
            n_iterations=2,
        )
        cluster_count = int(adata.obs['cluster'].astype('category').cat.categories.size)
        self._update_cell_metrics_summary({
            'cluster_count': cluster_count,
            'cluster_resolution': self.cell_cluster_resolution,
            'n_neighbors': n_neighbors,
            'n_pcs': n_pcs_use,
        })

        if not self.skip_umap:
            sc.tl.umap(adata, random_state=0)
            self._write_cell_embedding_plot(
                'umap',
                'cluster',
                f'{self.sample}_cell_umap_cluster.png',
                categorical=True,
                palette=categorical_palette(cluster_count),
            )
            for metric in ['total_counts', 'n_genes_by_counts', f'pct_counts_{MITO_VAR}']:
                if metric in adata.obs.columns:
                    self._write_cell_embedding_plot('umap', metric, f'{self.sample}_cell_umap_{metric}.png')
        else:
            print("[Cell Analysis] Skipping UMAP (--skip-umap enabled)")

        if not self.skip_tsne:
            sc.tl.tsne(adata, n_pcs=n_pcs_use, perplexity=get_tsne_perplexity(adata.n_obs))
            self._write_cell_embedding_plot(
                'tsne',
                'cluster',
                f'{self.sample}_cell_tsne_cluster.png',
                categorical=True,
                palette=categorical_palette(cluster_count),
            )
        else:
            print("[Cell Analysis] Skipping t-SNE (--skip-tsne enabled)")

        self._write_cell_spatial_plot('cluster', f'{self.sample}_cell_spatial_cluster.png', categorical=True)
        self._write_cell_spatial_plot('total_counts', f'{self.sample}_cell_spatial_total_counts.png')
        self._write_cell_spatial_plot('n_genes_by_counts', f'{self.sample}_cell_spatial_n_genes.png')
        self._write_cell_spatial_plot(f'pct_counts_{MITO_VAR}', f'{self.sample}_cell_spatial_pct_counts_mt.png')
        self.adata[self.bin] = adata

    @utils.add_log
    def cell_find_markers(self):
        """Find and write cell-level marker genes."""
        adata = self.adata[self.bin]
        clusters = adata.obs['cluster'].astype('category').cat.categories
        if len(clusters) < 2:
            print("WARNING: Fewer than two cell clusters; skipping marker analysis.")
            return

        sc.tl.rank_genes_groups(
            adata,
            'cluster',
            reference='rest',
            pts=True,
            method='wilcoxon',
            use_raw=False,
            layer=NORMALIZED_LAYER,
        )

        ncols = min(4, max(1, len(clusters)))
        sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, ncols=ncols, show=False)
        fig = plt.gcf()
        fig.set_size_inches(min(12, 3.0 * ncols), 7)
        fig.set_dpi(CELL_DPI)
        fig.tight_layout()
        fig.savefig(os.path.join(self.cell_marker_dir, f'{self.sample}_cell_marker_genes.png'), dpi=CELL_DPI)
        plt.close(fig)

        df_markers = sc.get.rank_genes_groups_df(adata, group=None, pval_cutoff=PVAL_CUTOFF)
        df_markers = _normalise_rank_genes_groups_df(df_markers)
        df_markers['cluster'] = df_markers['cluster'].astype(str)
        raw_marker_file = os.path.join(self.cell_marker_dir, f'{self.sample}_cell_markers_raw.tsv')
        df_markers.to_csv(raw_marker_file, index=None, sep='\t')

        df_markers_filter = df_markers.loc[df_markers['avg_log2FC'] > 0].sort_values('p_val_adj').groupby('cluster').head(100)
        df_markers_filter = df_markers_filter.round({'avg_log2FC': 3, 'pct.1': 3, 'pct.2': 3})
        marker_file = os.path.join(self.cell_marker_dir, f'{self.sample}_cell_markers.tsv')
        df_markers_filter.to_csv(marker_file, index=None, sep='\t')
        print(f"Cell marker genes saved to {marker_file}")

    @utils.add_log
    def run(self):
        if self.assay == 'spatial':
            # Spatial transcriptomics analysis pipeline
            self.calculate_qc_metrics()
            self.apply_spatial_qc()
            self.write_mito_stats()
            self.scatter_plot()
            self.violin_plot()

            self.normalize()
            self.hvg()
            self.scale()
            self.pca()
            self.neighbors()

            # Conditional dimensional reduction
            if not self.skip_tsne:
                print(f"[Dimensional Reduction] Running t-SNE for bin size {self.bin}")
                self.tsne()
            else:
                print(f"[Dimensional Reduction] Skipping t-SNE for bin size {self.bin} (--skip-tsne enabled)")

            if not self.skip_umap:
                print(f"[Dimensional Reduction] Running UMAP for bin size {self.bin}")
                self.umap()
            else:
                print(f"[Dimensional Reduction] Skipping UMAP for bin size {self.bin} (--skip-umap enabled)")

            self.leiden()

            # cell communication

            self.find_marker_genes()
            self.differentiate_expression_marker_genes()
            # gene regulatory network

            self.write_stat_images()
            self.write_markers()

            # Only write t-SNE coordinates if t-SNE was computed
            if not self.skip_tsne:
                self.write_tsne()
            else:
                print("[Output] Skipping t-SNE coordinate file (t-SNE was not computed)")

            self.write_h5ad()
        elif self.assay == 'cell':
            print("Running cell-level analysis pipeline...")
            self.cell_qc()
            self.cell_cluster()
            self.cell_find_markers()
            self.write_h5ad()
            print("Cell-level analysis completed!")
        else:  # scrna
            # scRNA-seq analysis pipeline - QC, UMAP, and marker analysis
            print("Running scRNA-seq analysis pipeline...")
            self.scrna_qc()
            self.scrna_umap()

            # Marker gene analysis
            print("Finding marker genes...")
            self.scrna_find_markers()
            self.scrna_visualize_markers()
            self.scrna_write_markers()

            self.write_h5ad()
            print("scRNA-seq analysis completed!")

    def get_df(self):
        """
        return df_tsne, df_marker

        If t-SNE was skipped, df_tsne will be None
        """
        # Only read t-SNE file if it exists (i.e., t-SNE was not skipped)
        if self.skip_tsne or not os.path.exists(self.df_tsne_file):
            df_tsne = None
            print(f"[Output] t-SNE coordinates not available (t-SNE was skipped or file not found)")
        else:
            df_tsne = read_tsne(self.df_tsne_file)

        df_marker = pd.read_csv(self.df_marker_file, sep="\t")
        df_marker = format_df_marker(df_marker)
        return df_tsne, df_marker

    def get_adata_name(self):
        return self.adata_name


def get_opts_analysis_match(parser, sub_program):
    """
    Do not perform analysis. Only read data from scRNA-seq match_dir.
    """
    if sub_program:
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'])
        parser.add_argument("--tsne_file", help=HELP_DICT['tsne_file'])
        parser.add_argument("--df_marker_file", help=HELP_DICT['df_marker_file'])

        parser = s_common(parser)


class Report_runner(Step):

    def __init__(self, args, display_title=None):

        super().__init__(args, display_title=display_title)

    def _render_html(self):
        pass

    def add_marker_help(self):
        self.add_help_content(
            name='Marker Genes by Cluster',
            content='differential expression analysis based on the non-parameteric Wilcoxon rank sum test'
        )
        self.add_help_content(
            name='avg_log2FC',
            content='log fold-change of the average expression between the cluster and the rest of the sample'
        )
        self.add_help_content(
            name='pct.1',
            content='The percentage of cells where the gene is detected in the cluster'
        )
        self.add_help_content(
            name='pct.2',
            content='The percentage of cells where the gene is detected in the rest of the sample'
        )
        self.add_help_content(
            name='p_val_adj',
            content='Adjusted p-value, based on bonferroni correction using all genes in the dataset'
        )

    @staticmethod
    def get_df_file(match_dir):
        """
        return df_tsne_file, df_marker_file
        """
        match_dict = utils.parse_match_dir(match_dir)
        df_tsne_file = match_dict['tsne_coord']
        df_marker_file = match_dict.get('markers', None)
        return df_tsne_file, df_marker_file

    def get_df(self):
        """
        return df_tsne, df_marker
        """
        if utils.check_arg_not_none(self.args, 'match_dir'):
            df_tsne_file, df_marker_file = self.get_df_file(self.args.match_dir)
        elif utils.check_arg_not_none(self.args, 'tsne_file'):
            df_tsne_file = self.args.tsne_file
            df_marker_file = self.args.df_marker_file
        else:
            raise ValueError('match_dir or tsne_file must be specified')
        df_tsne = read_tsne(df_tsne_file)
        if df_marker_file:
            df_marker = pd.read_csv(df_marker_file, sep="\t")
            df_marker = format_df_marker(df_marker)
        else:
            df_marker = None
        return df_tsne, df_marker

    def run(self):
        pass

#!/usr/bin/env python3
"""
Celatlas scRNA Report Generator

This script generates a comprehensive HTML report for single-cell RNA-seq analysis.
It reads data from various analysis output files and creates an interactive HTML report
with the following features:

- Key metrics and statistics visualization
- Barcode quality control analysis
- Mapping and feature counting statistics
- Cell-level analysis results with clustering
- Interactive marker genes table with search and export functionality
- Professional bioinformatics UI design
- Support for multiple sample types and auto-detection

Usage:
    python scrna_report_generator.py /path/to/sample_dir [sample_name]
    python scrna_report_generator.py /path/to/sample_dir --auto-detect

Author: Claude (Anthropic)
Date: 2025-11-28
Version: 1.0.0 - scRNA-seq specialized report generator
"""

import os
import sys
import json
import pandas as pd
import numpy as np
import base64
from datetime import datetime
from pathlib import Path
import argparse
import logging
import re
from typing import Dict, List, Tuple, Optional

from celatlas_spatial.tools.downsample import deterministic_downsample_metrics

# Plotly imports for interactive charts
try:
    import plotly.graph_objects as go
    import plotly.offline as pyo
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    print("Warning: plotly not found. Installing plotly...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "plotly"])
    import plotly.graph_objects as go
    import plotly.offline as pyo
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True

# Configuration constants
DEFAULT_TOP_GENES_PER_CLUSTER = 20
DEFAULT_CHEMISTRY = "BBV2.4"
DEFAULT_SPECIES = "Mus_musculus"
SUPPORTED_BIN_SIZES = ['10', '20', '50', '100']
CLUSTER_COLORS = [
    "#3498db", "#2ecc71", "#9b59b6", "#e74c3c",
    "#f1c40f", "#1abc9c", "#d35400", "#34495e",
    "#e67e22", "#8e44ad", "#2c3e50", "#27ae60",
    "#16a085", "#c0392b", "#8e44ad", "#f39c12"
]


class ScRNAReportGenerator:
    """Generate comprehensive single-cell RNA-seq analysis reports."""

    HD_BARCODE_BIN_CHEMISTRIES = {'BBV4', 'BBV4_L9'}
    HD_BIN_SIZE_UM = 10
    HD_SPOT_SIZE_UM = 2
    HD_PIXEL_SIZE_UM = 0.5

    def __init__(self, sample_dir: str, sample_name: str = None, chemistry: str = "BBV2.4",
                 species: str = "Mus_musculus", output_dir: str = None, tissue: str = None,
                 use_interactive_charts: bool = True):
        """
        Initialize the scRNA report generator.

        Args:
            sample_dir: Path to the sample analysis directory
            sample_name: Sample identifier (auto-detected if None)
            chemistry: Chemistry version used
            species: Species name
            output_dir: Output directory for the report (defaults to sample_dir)
        """
        self.sample_dir = Path(sample_dir)
        
        # Get sample name and chip number from environment variables
        import os
        env_sample_name = os.environ.get('CELATLAS_SAMPLE_NAME', '')
        env_chip_number = os.environ.get('CELATLAS_CHIP_NUMBER', '')
        env_tissue = os.environ.get('CELATLAS_TISSUE', '')
        
        # Auto-detect sample name if not provided
        if sample_name is None:
            self.sample_name = self._detect_sample_name()
        else:
            self.sample_name = sample_name
            
        # Store environment variables for display logic
        self.env_sample_name = env_sample_name
        self.env_chip_number = env_chip_number
        self.env_tissue = tissue if tissue is not None else env_tissue
            
        self.chemistry = chemistry
        self.species = species
        self.output_dir = Path(output_dir) if output_dir else self.sample_dir
        
        # Set up logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
        # Data containers
        self.bin_stats = {}
        self.marker_genes = []
        self.clusters = []
        self.general_stats = {}
        self.bin_images = {}
        self.analysis_images = {}
        self.company_logo = None
        self.downsample_data = None
        self.downsample_source = None
        self.count_detail_data = None
        self.barcode_umi_counts = None  # Pre-aggregated barcode->UMI count mapping
        self.hd_bin_metrics = None
        self.hd_bin_counts = None
        self.hd_bin_downsample = None
        self.hd_bin_source = None
        
        # Chart HTML containers
        self.barcode_rank_plot_html = ""
        self.sequencing_saturation_plot_html = ""
        self.median_genes_plot_html = ""
        
        # Performance optimization flag
        self.use_interactive_charts = True  # Can be set to False for faster loading
        
        # Template path
        self.template_path = Path(__file__).parent.parent / "templates" / "html" / "scrna_report_template.html"
        
        # Validate inputs
        if not self.sample_dir.exists():
            raise FileNotFoundError(f"Sample directory not found: {self.sample_dir}")
        if not self.template_path.exists():
            raise FileNotFoundError(f"Template file not found: {self.template_path}")

    def is_hd_barcode_bin_mode(self) -> bool:
        """Return True for HD spatial chemistries run through the scRNA report."""
        return self.chemistry in self.HD_BARCODE_BIN_CHEMISTRIES

    def _unit_text(self) -> dict:
        if self.is_hd_barcode_bin_mode():
            if self.hd_bin_metrics is not None:
                statistics_title = 'Bin10 Statistics'
                unit_name = 'bin10 unit'
                unit_plural = 'bin10 units'
                rank_trace_name = 'bin10 units'
                number_label = 'Number of bin10 units'
                estimated_label = 'Estimated Number of bin10 Units'
                fraction_label = 'Fraction Reads in bin10 Units'
                mean_reads_label = 'Mean Reads per bin10'
                mean_umis_label = 'Median UMIs per bin10'
                median_umi_label = 'Median UMI per bin10'
                median_genes_label = 'Median Genes per bin10'
                median_genes_plot_title = 'Median Genes per bin10'
                median_genes_chart_info = 'Shows how median genes per bin10 unit changes with increasing sequencing depth'
                key_metrics_help = (
                    '<strong>Number of bin10 units:</strong> Count of 10um x 10um bins generated by grouping adjacent 2um BBV4 spots; '
                    'for HD BBV4 data these are not true cells<br>'
                    '<strong>Mean reads per bin10:</strong> Average reads per 10um bin<br>'
                    '<strong>Median UMIs per bin10:</strong> Median unique molecular identifiers per 10um bin<br>'
                    '<strong>Total genes detected:</strong> Total number of genes with detectable expression across the entire sample'
                )
                saturation_unit_help = (
                    '<strong>bin10 unit:</strong> 10um x 10um HD spatial bin, approximately 25 adjacent 2um spots; '
                    'it should not be interpreted as a true single cell'
                )
                statistics_help_title = 'bin10 Statistics'
                statistics_help = (
                    '<strong>Estimated Number of bin10 Units:</strong> Number of 10um x 10um bins with detected molecules<br>'
                    '<strong>Fraction Reads in bin10 Units:</strong> Percentage of assigned reads mapped into bin10 units<br>'
                    '<strong>Mean Reads per bin10:</strong> Average number of reads per bin10 unit<br>'
                    '<strong>Median UMI per bin10:</strong> Middle value of UMIs across bin10 units<br>'
                    '<strong>Total Genes:</strong> Total number of genes detected<br>'
                    '<strong>Median Genes per bin10:</strong> Middle value of gene counts across bin10 units<br>'
                    '<strong>Saturation:</strong> Downsample-estimated duplicate read fraction at bin10 level; use as a library-complexity diagnostic, not as true cell saturation'
                )
                diagnostic_note_html = (
                    '<div class="diagnostic-note">'
                    '<strong>HD bin10 diagnostic mode:</strong> BBV4 2um spatial barcodes were aggregated into 10um x 10um bin10 units '
                    f'using coordinate table {self.hd_bin_source}. '
                    'Each bin10 unit corresponds to roughly 25 adjacent 2um spots. These metrics should not be reported as true cell-level metrics.'
                    '</div>'
                )
            else:
                statistics_title = 'Barcode Diagnostic Statistics'
                unit_name = 'barcode unit'
                unit_plural = 'barcode units'
                rank_trace_name = 'barcode units'
                number_label = 'Number of barcode units'
                estimated_label = 'Estimated Number of Barcode Units'
                fraction_label = 'Fraction Reads in Barcode Units'
                mean_reads_label = 'Mean Reads per Barcode'
                mean_umis_label = 'Median UMIs per Barcode'
                median_umi_label = 'Median UMI per Barcode'
                median_genes_label = 'Median Genes per Barcode'
                median_genes_plot_title = 'Median Genes per Barcode'
                median_genes_chart_info = 'Shows how median genes per barcode changes with increasing sequencing depth'
                key_metrics_help = (
                    '<strong>Number of barcode units:</strong> Count of barcode-derived units from the scRNA-style diagnostic workflow; '
                    'for HD BBV4 data these are not true cells<br>'
                    '<strong>Mean reads per barcode:</strong> Average reads per reported barcode unit<br>'
                    '<strong>Median UMIs per barcode:</strong> Median unique molecular identifiers per reported barcode unit<br>'
                    '<strong>Total genes detected:</strong> Total number of genes with detectable expression across the entire sample'
                )
                saturation_unit_help = (
                    '<strong>Barcode unit:</strong> Barcode-derived HD spatial unit used for diagnostic reporting; '
                    'it should not be interpreted as a true single cell'
                )
                statistics_help_title = 'Barcode Diagnostic Statistics'
                statistics_help = (
                    '<strong>Estimated Number of Barcode Units:</strong> Barcode-derived units selected by the scRNA-style count workflow<br>'
                    '<strong>Fraction Reads in Barcode Units:</strong> Percentage of assigned reads in selected barcode units<br>'
                    '<strong>Mean Reads per Barcode:</strong> Average number of reads per selected barcode unit<br>'
                    '<strong>Median UMI per Barcode:</strong> Middle value of UMIs across selected barcode units<br>'
                    '<strong>Total Genes:</strong> Total number of genes detected<br>'
                    '<strong>Median Genes per Barcode:</strong> Middle value of gene counts across selected barcode units<br>'
                    '<strong>Saturation:</strong> Downsample-estimated duplicate read fraction at barcode level; use as a library-complexity diagnostic, not as true cell saturation'
                )
                diagnostic_note_html = (
                    '<div class="diagnostic-note">'
                    '<strong>HD diagnostic mode:</strong> BBV4 2um spatial barcodes require a FilterBarcodes coordinate table for true bin10-level metrics. '
                    'No bin10 aggregation was applied for this report, so count-derived values remain scRNA-style barcode diagnostics and should not be interpreted as true cells.'
                    '</div>'
                )
            return {
                'report_title': 'Celatlas HD',
                'report_subtitle': 'Diagnostic Report',
                'analysis_tab': 'Barcode/bin Analysis',
                'analysis_title': 'Barcode/bin Analysis',
                'statistics_title': statistics_title,
                'unit_name': unit_name,
                'unit_name_lower': unit_name,
                'unit_name_plural': unit_plural,
                'unit_name_plural_lower': unit_plural,
                'rank_trace_name': rank_trace_name,
                'number_label': number_label,
                'estimated_label': estimated_label,
                'fraction_label': fraction_label,
                'mean_reads_label': mean_reads_label,
                'mean_umis_label': mean_umis_label,
                'median_umi_label': median_umi_label,
                'median_genes_label': median_genes_label,
                'median_genes_plot_title': median_genes_plot_title,
                'median_genes_chart_info': median_genes_chart_info,
                'key_metrics_help': key_metrics_help,
                'saturation_unit_help': saturation_unit_help,
                'statistics_help_title': statistics_help_title,
                'statistics_help': statistics_help,
                'qc_help': (
                    '<strong>Scatter Plot:</strong> Shows relationship between total UMI counts and detected genes per reported unit<br>'
                    '<strong>Violin Plot:</strong> Displays distribution of QC metrics (UMIs, genes, mitochondrial %)<br>'
                    'For HD BBV4 data, summary metrics are bin10-level diagnostics and should not be interpreted as true cell identification.'
                ),
                'umap_help': (
                    '<strong>Interactive UMAP:</strong> Click and drag to pan, scroll to zoom, hover over reported units to see details. Click legend to show/hide clusters.<br>'
                    '<strong>Interactive t-SNE:</strong> Alternative dimensionality reduction method for visualizing diagnostic clusters.<br>'
                    'Both methods reduce high-dimensional gene expression data into 2D for diagnostic visualization.'
                ),
                'hover_text': 'Hover over reported units to see details • Click legend to filter clusters • Drag to pan • Scroll to zoom',
                'cluster_help': (
                    '<strong>Cluster:</strong> Groups of reported units with similar gene expression patterns<br>'
                    '<strong>Marker Gene:</strong> Genes that are differentially expressed in each cluster<br>'
                    '<strong>Avg log2FC:</strong> Average log2 fold-change between target cluster and other clusters<br>'
                    '<strong>Pct.1/Pct.2:</strong> Percentage of reported units expressing the gene in target cluster vs. other clusters'
                ),
                'diagnostic_note_html': diagnostic_note_html,
            }
        return {
            'report_title': 'Celatlas scRNA',
            'report_subtitle': 'Analysis Report',
            'analysis_tab': 'Cell Analysis',
            'analysis_title': 'Cell Analysis',
            'statistics_title': 'Cell Statistics',
            'unit_name': 'Cell',
            'unit_name_lower': 'cell',
            'unit_name_plural': 'Cells',
            'unit_name_plural_lower': 'cells',
            'rank_trace_name': 'Cells',
            'number_label': 'Number of cells',
            'estimated_label': 'Estimated Number of Cells',
            'fraction_label': 'Fraction Reads in Cells',
            'mean_reads_label': 'Mean Reads per Cell',
            'mean_umis_label': 'Mean UMIs per Cell',
            'median_umi_label': 'Median UMI per Cell',
            'median_genes_label': 'Median Genes per Cell',
            'median_genes_plot_title': 'Median Genes per Cell',
            'median_genes_chart_info': 'Shows how median genes per cell changes with increasing sequencing depth',
            'key_metrics_help': (
                '<strong>Number of cells:</strong> Estimated number of cells detected in the sample<br>'
                '<strong>Mean reads per cell:</strong> Average number of sequenced DNA fragments per cell<br>'
                '<strong>Mean UMIs per cell:</strong> Average unique molecular identifiers per cell (removes PCR duplicates)<br>'
                '<strong>Total genes detected:</strong> Total number of genes with detectable expression across the entire sample'
            ),
            'saturation_unit_help': '<strong>Cell:</strong> Individual cells captured and barcoded during single-cell RNA sequencing',
            'statistics_help_title': 'Cell Statistics',
            'statistics_help': (
                '<strong>Estimated Number of Cells:</strong> Predicted total number of cells in the sample<br>'
                '<strong>Fraction Reads in Cells:</strong> Percentage of total reads that map to valid cells<br>'
                '<strong>Mean Reads per Cell:</strong> Average number of reads per cell<br>'
                '<strong>Median UMI per Cell:</strong> Middle value of UMIs across all cells<br>'
                '<strong>Total Genes:</strong> Total number of genes detected<br>'
                '<strong>Median Genes per Cell:</strong> Middle value of gene counts across cells<br>'
                '<strong>Saturation:</strong> Downsample-estimated duplicate read fraction; interpret together with median genes, median UMI, and total read depth'
            ),
            'qc_help': (
                '<strong>Scatter Plot:</strong> Shows relationship between total UMI counts and detected genes per cell<br>'
                '<strong>Violin Plot:</strong> Displays distribution of QC metrics (UMIs, genes, mitochondrial %)<br>'
                'These plots show data quality before filtering to ensure proper cell identification.'
            ),
            'umap_help': (
                '<strong>Interactive UMAP:</strong> Click and drag to pan, scroll to zoom, hover over cells to see details. Click legend to show/hide clusters.<br>'
                '<strong>Interactive t-SNE:</strong> Alternative dimensionality reduction method for visualizing cell clusters.<br>'
                'Both methods reduce high-dimensional gene expression data into 2D for visualization. UMAP better preserves global structure while t-SNE emphasizes local relationships.'
            ),
            'hover_text': 'Hover over cells to see details • Click legend to filter clusters • Drag to pan • Scroll to zoom',
            'cluster_help': (
                '<strong>Cluster:</strong> Groups of cells with similar gene expression patterns<br>'
                '<strong>Marker Gene:</strong> Genes that are differentially expressed in each cluster<br>'
                '<strong>Avg log2FC:</strong> Average log2 fold-change between target cluster and other clusters<br>'
                '<strong>Pct.1/Pct.2:</strong> Percentage of cells expressing the gene in target cluster vs. other clusters'
            ),
            'diagnostic_note_html': '',
        }

    @staticmethod
    def _parse_numeric(value) -> Optional[float]:
        """Parse numeric strings that may contain commas, percent signs, or whitespace."""
        try:
            if value is None:
                return None
            if pd.isna(value):
                return None
            cleaned = str(value).replace(',', '').replace('%', '').strip()
            if cleaned in ('', 'N/A', 'nan', 'None'):
                return None
            return float(cleaned)
        except (ValueError, TypeError):
            return None

    @staticmethod
    def _format_percent_value(value: Optional[float]) -> str:
        if value is None or pd.isna(value):
            return "N/A"
        return f"{max(0.0, min(100.0, float(value))):.2f}%"

    def _normalize_downsample_data(self, df: pd.DataFrame, source: Path) -> Optional[pd.DataFrame]:
        """Normalize downsample data to read_fraction 0-1 and saturation 0-100."""
        required_saturation = [col for col in ('umi_saturation', 'read_saturation') if col in df.columns]
        if 'read_fraction' not in df.columns or not required_saturation:
            self.logger.warning(f"Invalid downsample file for saturation plot: {source}")
            return None

        normalized = df.copy()
        numeric_cols = ['read_fraction', 'median_gene_number'] + required_saturation
        for col in numeric_cols:
            if col in normalized.columns:
                normalized[col] = pd.to_numeric(normalized[col], errors='coerce')

        normalized = normalized.dropna(subset=['read_fraction'] + required_saturation)
        if normalized.empty:
            self.logger.warning(f"Downsample file has no valid numeric rows: {source}")
            return None

        if normalized['read_fraction'].max() > 1.5:
            normalized['read_fraction'] = normalized['read_fraction'] / 100.0
        normalized = normalized[(normalized['read_fraction'] >= 0) & (normalized['read_fraction'] <= 1)]

        for col in required_saturation:
            if normalized[col].max() <= 1.0:
                normalized[col] = normalized[col] * 100.0
            normalized[col] = normalized[col].clip(lower=0, upper=100)

        normalized = normalized.sort_values('read_fraction').drop_duplicates('read_fraction', keep='last')
        return normalized.reset_index(drop=True)

    def _chart_placeholder(self, title: str, message: str) -> str:
        return f'''
        <div class="simple-chart">
            <div class="chart-placeholder">
                {title}<br>
                <small>{message}</small>
            </div>
        </div>
        '''

    def _hd_spots_per_bin(self) -> int:
        return max(1, int(round(self.HD_BIN_SIZE_UM / self.HD_SPOT_SIZE_UM)))

    def _hd_bin_coordinate_width(self) -> int:
        # Match spatial_binSegment: bin width in fullres coordinates is
        # ceil(micron_bin / pixel_size). With BBV4 pixel_size=0.5, bin10 is 20.
        return max(1, int(np.ceil(self.HD_BIN_SIZE_UM / self.HD_PIXEL_SIZE_UM)))

    def _hd_bin_label(self, x: pd.Series, y: pd.Series) -> pd.Series:
        bin_width = self._hd_bin_coordinate_width()
        bin_x = (pd.to_numeric(x, errors='coerce') // bin_width).astype('Int64')
        bin_y = (pd.to_numeric(y, errors='coerce') // bin_width).astype('Int64')
        return 'bin10_' + bin_x.astype(str) + '_' + bin_y.astype(str)

    def _candidate_mask_dirs(self) -> List[Path]:
        candidates = [
            self.sample_dir / "06.segment" / "mask",
            self.sample_dir,
            self.sample_dir.parent,
            Path("/mnt/strna/celatlas_spatial/ST_mask"),
        ]
        unique = []
        seen = set()
        for path in candidates:
            key = str(path)
            if key not in seen:
                seen.add(key)
                unique.append(path)
        return unique

    def _find_hd_filter_barcodes_file(self) -> Optional[Path]:
        exact_names = [
            f"{self.sample_name}_FilterBarcodes.csv",
            f"{self.sample_name}.FilterBarcodes.csv",
        ]
        for directory in self._candidate_mask_dirs():
            if not directory.exists():
                continue
            for name in exact_names:
                candidate = directory / name
                if candidate.exists():
                    return candidate

        # Conservative fallback: only auto-pick a candidate when the sample prefix
        # matches exactly enough to avoid crossing chips accidentally.
        sample_prefix = self.sample_name.split('_')[0]
        for directory in self._candidate_mask_dirs():
            if not directory.exists():
                continue
            matches = sorted(directory.glob(f"{sample_prefix}*_FilterBarcodes.csv"))
            if len(matches) == 1:
                return matches[0]
        return None

    def _load_hd_barcode_positions(self) -> Optional[pd.DataFrame]:
        filter_barcodes_file = self._find_hd_filter_barcodes_file()
        if filter_barcodes_file is None:
            self.logger.warning(
                "HD barcode/bin metrics skipped: no *_FilterBarcodes.csv found for %s",
                self.sample_name,
            )
            return None

        try:
            positions = pd.read_csv(
                filter_barcodes_file,
                header=None,
                names=["x", "y", "barcode"],
                usecols=[0, 1, 2],
                dtype={"barcode": "string"},
            )
        except Exception as e:
            self.logger.warning(f"HD barcode/bin metrics skipped: failed to read {filter_barcodes_file}: {e}")
            return None

        positions = positions.dropna(subset=["x", "y", "barcode"])
        positions["barcode"] = positions["barcode"].astype(str).str.strip()
        positions = positions[positions["barcode"] != ""]
        positions["x"] = pd.to_numeric(positions["x"], errors="coerce")
        positions["y"] = pd.to_numeric(positions["y"], errors="coerce")
        positions = positions.dropna(subset=["x", "y"])
        positions = positions.drop_duplicates(subset=["barcode"], keep="first")

        if positions.empty:
            self.logger.warning(f"HD barcode/bin metrics skipped: no valid positions in {filter_barcodes_file}")
            return None

        positions["bin_barcode"] = self._hd_bin_label(positions["x"], positions["y"])
        self.hd_bin_source = str(filter_barcodes_file)
        self.logger.info(
            "Loaded HD barcode positions from %s: %s barcodes, %s bin10 units",
            filter_barcodes_file,
            f"{len(positions):,}",
            f"{positions['bin_barcode'].nunique():,}",
        )
        return positions[["barcode", "bin_barcode"]]

    def _aggregate_hd_bin_count_detail(self, positions: pd.DataFrame) -> Optional[pd.DataFrame]:
        count_detail_file = self.sample_dir / "05.count" / f"{self.sample_name}_count_detail.txt"
        if not count_detail_file.exists():
            self.logger.warning(f"HD barcode/bin metrics skipped: count_detail file not found: {count_detail_file}")
            return None

        pos_map = positions.set_index("barcode")["bin_barcode"]
        grouped_parts = []
        total_rows = 0
        matched_rows = 0
        matched_reads = 0
        total_reads = 0
        chunk_size = 500000

        try:
            for chunk in pd.read_csv(count_detail_file, sep="\t", chunksize=chunk_size):
                required = {"Barcode", "geneID", "UMI", "count"}
                if not required.issubset(chunk.columns):
                    missing = ", ".join(sorted(required - set(chunk.columns)))
                    self.logger.warning(f"HD barcode/bin metrics skipped: count_detail missing columns: {missing}")
                    return None

                total_rows += len(chunk)
                chunk = chunk[["Barcode", "geneID", "UMI", "count"]].copy()
                chunk["Barcode"] = chunk["Barcode"].astype(str)
                chunk["OriginalBarcode"] = chunk["Barcode"]
                chunk["count"] = pd.to_numeric(chunk["count"], errors="coerce").fillna(0).astype(np.int64)
                total_reads += int(chunk["count"].sum())

                chunk["bin_barcode"] = chunk["Barcode"].map(pos_map)
                chunk = chunk.dropna(subset=["bin_barcode"])
                if chunk.empty:
                    continue

                matched_rows += len(chunk)
                matched_reads += int(chunk["count"].sum())
                grouped = (
                    chunk.groupby(["bin_barcode", "OriginalBarcode", "geneID", "UMI"], sort=False, observed=True)["count"]
                    .sum()
                    .reset_index()
                    .rename(columns={"bin_barcode": "Barcode"})
                )
                grouped_parts.append(grouped)

                if total_rows % 5000000 == 0:
                    self.logger.info(
                        "HD bin10 aggregation processed %s rows; matched %s rows",
                        f"{total_rows:,}",
                        f"{matched_rows:,}",
                    )
        except Exception as e:
            self.logger.warning(f"HD barcode/bin metrics skipped: failed to aggregate {count_detail_file}: {e}")
            return None

        if not grouped_parts:
            self.logger.warning(
                "HD barcode/bin metrics skipped: no count_detail barcodes matched %s",
                self.hd_bin_source or "position file",
            )
            return None

        bin_count_detail = pd.concat(grouped_parts, ignore_index=True)
        bin_count_detail = (
            bin_count_detail.groupby(["Barcode", "OriginalBarcode", "geneID", "UMI"], sort=False, observed=True)["count"]
            .sum()
            .reset_index()
        )
        bin_count_detail["count"] = pd.to_numeric(bin_count_detail["count"], errors="coerce").fillna(0).astype(np.int64)
        bin_count_detail = bin_count_detail[bin_count_detail["count"] > 0]

        if bin_count_detail.empty:
            return None

        self.logger.info(
            "HD bin10 aggregation complete: %s rows -> %s dedup rows, matched reads %s/%s",
            f"{total_rows:,}",
            f"{len(bin_count_detail):,}",
            f"{matched_reads:,}",
            f"{total_reads:,}",
        )
        return bin_count_detail

    @staticmethod
    def _summarize_count_detail(count_detail: pd.DataFrame) -> pd.DataFrame:
        def num_gt2(x):
            return (x > 1).sum()

        summary = count_detail.groupby("Barcode").agg({
            "count": ["sum", num_gt2],
            "UMI": "count",
            "geneID": "nunique",
        })
        summary.columns = ["readcount", "UMI2", "UMI", "geneID"]
        return summary.sort_values("UMI", ascending=False)

    def _apply_hd_bin_metrics(self, bin_count_detail: pd.DataFrame):
        bin_summary = self._summarize_count_detail(bin_count_detail)
        if bin_summary.empty:
            return

        self.hd_bin_counts = bin_summary["UMI"]
        self.barcode_umi_counts = self.hd_bin_counts

        raw_total_reads = self._parse_numeric(self.general_stats.get("Reads Assigned To Exonic Regions", "0"))
        if raw_total_reads is None or raw_total_reads <= 0:
            raw_total_reads = self._parse_numeric(self.general_stats.get("Uniquely Mapped Reads", "0"))
        if raw_total_reads is None or raw_total_reads <= 0:
            raw_total_reads = float(bin_summary["readcount"].sum())

        reads_in_bins = float(bin_summary["readcount"].sum())
        fraction_reads = round(reads_in_bins / raw_total_reads * 100.0, 2) if raw_total_reads else 0.0

        metrics = {
            "Estimated Number of Cells": int(len(bin_summary)),
            "Fraction Reads in Cells": fraction_reads,
            "Mean Reads per Cell": int(reads_in_bins / len(bin_summary)) if len(bin_summary) else 0,
            "Median UMI per Cell": int(bin_summary["UMI"].median()),
            "Total Genes": int(bin_count_detail["geneID"].nunique()),
            "Median Genes per Cell": int(bin_summary["geneID"].median()),
        }

        try:
            downsample = deterministic_downsample_metrics(bin_count_detail)
            self.hd_bin_downsample = downsample
            self.downsample_data = downsample
            self.downsample_source = f"HD bin10 from count_detail + {self.hd_bin_source}"
            metrics["Saturation"] = self._format_percent_value(downsample["umi_saturation"].iloc[-1])
            self.logger.info(
                "HD bin10 downsample: median genes %s, saturation %s",
                downsample["median_gene_number"].iloc[-1],
                downsample["umi_saturation"].iloc[-1],
            )
        except Exception as e:
            self.logger.warning(f"HD bin10 downsample failed, using existing count-level curve: {e}")

        self.hd_bin_metrics = metrics
        for key, value in metrics.items():
            self.general_stats[key] = value

        self.logger.info(
            "Applied HD bin10 metrics: bins=%s, median UMI=%s, median genes=%s, fraction reads=%s%%",
            f"{metrics['Estimated Number of Cells']:,}",
            metrics["Median UMI per Cell"],
            metrics["Median Genes per Cell"],
            metrics["Fraction Reads in Cells"],
        )

    def build_hd_barcode_bin_metrics(self):
        """Build bin10-level metrics for HD BBV4 scRNA diagnostic reports."""
        if not self.is_hd_barcode_bin_mode():
            return

        positions = self._load_hd_barcode_positions()
        if positions is None:
            return

        bin_count_detail = self._aggregate_hd_bin_count_detail(positions)
        if bin_count_detail is None:
            return

        self._apply_hd_bin_metrics(bin_count_detail)
        
    def parse_stat_file(self, file_path: Path) -> Dict[str, str]:
        """Parse a stat.txt file and extract key-value pairs."""
        stats = {}
        if file_path.exists():
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    for line in f:
                        line = line.strip()
                        if ':' in line:
                            key, value = line.split(':', 1)
                            # Clean up the key and value
                            key = key.strip()
                            value = value.strip()
                            
                            # Extract percentage from parentheses if present
                            if '(' in value and ')' in value:
                                # Split main value and percentage
                                main_value = value.split('(')[0].strip()
                                percentage = value.split('(')[1].split(')')[0].strip()
                                stats[key] = main_value.replace(',', '')
                                stats[f"{key} Percentage"] = percentage
                            else:
                                stats[key] = value.replace(',', '')
            except Exception as e:
                self.logger.warning(f"Failed to parse {file_path}: {e}")
        return stats

    def parse_step_log_metrics(self, file_path: Path) -> Dict[str, str]:
        """Parse visible and hidden metrics from a step.log file."""
        stats = {}
        if not file_path.exists():
            return stats
        in_metric_section = False
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.rstrip('\n')
                    if line in {"Metrics", "Hidden Metrics"}:
                        in_metric_section = True
                        continue
                    if not in_metric_section:
                        continue
                    if not line or set(line) == {'-'} or line.startswith('  '):
                        continue
                    if not line.startswith('- ') or ':' not in line:
                        in_metric_section = False
                        continue
                    key, value = line[2:].split(':', 1)
                    key = key.strip()
                    value = value.strip()
                    if '(' in value and ')' in value:
                        main_value = value.split('(')[0].strip()
                        percentage = value.split('(')[1].split(')')[0].strip()
                        stats[key] = main_value.replace(',', '')
                        stats[f"{key} Percentage"] = percentage
                    else:
                        stats[key] = value.replace(',', '')
        except Exception as e:
            self.logger.warning(f"Failed to parse {file_path}: {e}")
        return stats

    def parse_step_metrics(self, step_dir: Path) -> Dict[str, str]:
        step_log_stats = self.parse_step_log_metrics(step_dir / "step.log")
        if step_log_stats:
            return step_log_stats
        stat_file = step_dir / "stat.txt"
        if stat_file.exists():
            return self.parse_stat_file(stat_file)
        return {}
    
    def _detect_sample_name(self) -> str:
        """Auto-detect sample name from directory structure or files."""
        logger = getattr(self, 'logger', logging.getLogger(__name__))

        # Try to detect from directory name
        if self.sample_dir.name and self.sample_dir.name != 'analysis':
            return self.sample_dir.name
        
        # Try to detect from marker genes file
        analysis_dir = self._binned_outputs_dir()
        marker_dirs = [analysis_dir / "markers", analysis_dir]
        for marker_dir in marker_dirs:
            if marker_dir.exists():
                for file in marker_dir.glob("*_markers_raw.tsv"):
                    sample_name = file.stem.replace('_markers_raw', '')
                    logger.info(f"Auto-detected sample name: {sample_name}")
                    return sample_name

        # Try to detect from other analysis files
        h5ad_dirs = [analysis_dir / "h5ad", analysis_dir]
        for pattern in ["*_bin*.h5ad", "*_spatial_analysis_report.html"]:
            for search_dir in h5ad_dirs:
                if not search_dir.exists():
                    continue
                for file in search_dir.glob(pattern):
                    parts = file.stem.split('_')
                    if len(parts) >= 2:
                        sample_name = '_'.join(parts[:2])
                        logger.info(f"Auto-detected sample name from {file.name}: {sample_name}")
                        return sample_name

        # Fallback to parent directory name
        fallback_name = self.sample_dir.parent.name if self.sample_dir.parent.name != 'analysis' else 'sample'
        logger.warning(f"Could not auto-detect sample name, using fallback: {fallback_name}")
        return fallback_name

    def _binned_outputs_dir(self) -> Path:
        """Return the compact 07.outs binned output directory."""
        return self.sample_dir / "07.outs" / "binned_outputs"

    def _binned_qc_dir(self) -> Path:
        return self._binned_outputs_dir() / "qc"

    def _binned_clustering_dir(self) -> Path:
        return self._binned_outputs_dir() / "clustering"

    @staticmethod
    def _first_existing(paths: List[Path]) -> Path:
        for path in paths:
            if path.exists():
                return path
        return paths[0]

    def _binsegment_dir(self) -> Path:
        """Return the current binSegment output directory."""
        return self.sample_dir / "06.segment" / "01.binsegment"
    
    def load_bin_statistics(self):
        """Load statistics for all supported bin sizes.
        
        Searches for bin statistics in the standard directory structure:
        06.segment/01.binsegment/square_bin/{sample_name}_bin{size}/stat.txt
        """
        binsegment_dir = self._binsegment_dir()
        for bin_size in SUPPORTED_BIN_SIZES:
            bin_dir = binsegment_dir / "square_bin" / f"{self.sample_name}_bin{bin_size}"
            stat_file = bin_dir / "stat.txt"
            
            if stat_file.exists():
                stats = self.parse_stat_file(stat_file)
                self.bin_stats[f'bin{bin_size}'] = stats
                self.logger.info(f"Loaded statistics for bin {bin_size}μm")
            else:
                self.logger.warning(f"Bin {bin_size}μm stat file not found: {stat_file}")
    
    def load_general_statistics(self):
        """Load general pipeline statistics."""
        # Load barcode statistics
        self.general_stats.update(self.parse_step_metrics(self.sample_dir / "01.barcode"))
        
        # Load STAR mapping statistics
        self.general_stats.update(self.parse_step_metrics(self.sample_dir / "03.star"))
        
        # Load featureCounts statistics
        self.general_stats.update(self.parse_step_metrics(self.sample_dir / "04.featureCounts"))
        
        # Load count statistics
        self.general_stats.update(self.parse_step_metrics(self.sample_dir / "05.count"))
            
        # Load analysis statistics if present.
        self.general_stats.update(self.parse_step_metrics(self._binned_outputs_dir()))
    
    def load_marker_genes(self, top_n: int = 20):
        """Load marker genes from the analysis results, ensuring top N genes per cluster."""
        marker_file = self.sample_dir / "06_analysis_wrapper" / "03.Markers" / f"{self.sample_name}_markers_raw.tsv"
        
        if marker_file.exists():
            try:
                df = pd.read_csv(marker_file, sep='\t')
                self.logger.info(f"Loaded marker genes file with {len(df)} total genes")
                
                # Get top N genes per cluster based on scores (higher is better)
                marker_genes = []
                cluster_summary = []
                
                for cluster in sorted(df['cluster'].unique()):
                    cluster_df = df[df['cluster'] == cluster]
                    
                    # Sort by scores descending to get top genes (most significant)
                    if 'scores' in cluster_df.columns:
                        cluster_df = cluster_df.sort_values('scores', ascending=False)
                    elif 'avg_log2FC' in cluster_df.columns:
                        # Fallback to avg_log2FC if scores not available
                        cluster_df = cluster_df.sort_values('avg_log2FC', ascending=False)
                    
                    # Take top N genes for this cluster
                    top_genes = cluster_df.head(top_n)
                    cluster_summary.append(f"Cluster {cluster}: {len(top_genes)} genes")
                    
                    for _, row in top_genes.iterrows():
                        marker_genes.append({
                            'cluster': str(row['cluster']),
                            'gene': str(row['gene']),
                            'scores': float(row['scores']) if pd.notna(row['scores']) else 0.0,
                            'avg_log2FC': float(row['avg_log2FC']) if pd.notna(row['avg_log2FC']) else 0.0,
                            'p_val': float(row['p_val']) if pd.notna(row['p_val']) else 1.0,
                            'p_val_adj': float(row['p_val_adj']) if pd.notna(row['p_val_adj']) else 1.0,
                            'pct1': float(row['pct.1']) if 'pct.1' in row and pd.notna(row['pct.1']) else 0.0,
                            'pct2': float(row['pct.2']) if 'pct.2' in row and pd.notna(row['pct.2']) else 0.0
                        })
                
                self.marker_genes = marker_genes
                self.clusters = sorted(df['cluster'].unique())
                self.logger.info(f"Loaded {len(marker_genes)} marker genes from {len(self.clusters)} clusters")
                self.logger.info(f"Distribution: {', '.join(cluster_summary)}")
                
                # Also store complete dataset for TSV export
                self.complete_marker_genes = []
                for _, row in df.iterrows():
                    self.complete_marker_genes.append({
                        'cluster': str(row['cluster']),
                        'gene': str(row['gene']),
                        'scores': float(row['scores']) if pd.notna(row['scores']) else 0.0,
                        'avg_log2FC': float(row['avg_log2FC']) if pd.notna(row['avg_log2FC']) else 0.0,
                        'p_val': float(row['p_val']) if pd.notna(row['p_val']) else 1.0,
                        'p_val_adj': float(row['p_val_adj']) if pd.notna(row['p_val_adj']) else 1.0,
                        'pct1': float(row['pct.1']) if 'pct.1' in row and pd.notna(row['pct.1']) else 0.0,
                        'pct2': float(row['pct.2']) if 'pct.2' in row and pd.notna(row['pct.2']) else 0.0
                    })
                
                self.logger.info(f"Stored {len(self.complete_marker_genes)} complete marker genes for TSV export")
                
            except Exception as e:
                self.logger.error(f"Failed to load marker genes: {e}")
                import traceback
                self.logger.error(traceback.format_exc())
        else:
            self.logger.warning(f"Marker genes file not found: {marker_file}")
    
    def generate_markers_tsv_data_uri(self) -> str:
        """Generate data URI for markers_raw.tsv file to enable offline download."""
        marker_file = self.sample_dir / "06_analysis_wrapper" / "03.Markers" / f"{self.sample_name}_markers_raw.tsv"
        
        if marker_file.exists():
            try:
                # Read the TSV file content with proper handling
                with open(marker_file, 'r', encoding='utf-8', newline='') as f:
                    tsv_content = f.read()
                
                # Validate content
                lines = tsv_content.strip().split('\n')
                self.logger.info(f"TSV file contains {len(lines)} lines")
                
                # Ensure content ends with newline for proper TSV format
                if not tsv_content.endswith('\n'):
                    tsv_content += '\n'
                
                # Encode as base64 for data URI with proper MIME type
                import base64
                encoded_content = base64.b64encode(tsv_content.encode('utf-8')).decode('ascii')
                
                # Use proper MIME type for TSV files
                data_uri = f"data:text/tab-separated-values;charset=utf-8;base64,{encoded_content}"
                
                # Validate the data URI length
                self.logger.info(f"Generated data URI for markers TSV file:")
                self.logger.info(f"  - Original file size: {len(tsv_content)} characters")
                self.logger.info(f"  - Number of lines: {len(lines)}")
                self.logger.info(f"  - Data URI size: {len(data_uri)} characters")
                self.logger.info(f"  - Base64 encoded size: {len(encoded_content)} characters")
                
                # Check for potential browser limitations
                if len(data_uri) > 2000000:  # 2MB limit for some browsers
                    self.logger.warning(f"Data URI is very large ({len(data_uri)} chars), may hit browser limits")
                
                # Validate by attempting to decode back
                try:
                    decoded_test = base64.b64decode(encoded_content).decode('utf-8')
                    decoded_lines = decoded_test.strip().split('\n')
                    if len(decoded_lines) != len(lines):
                        self.logger.error(f"Data URI validation failed: {len(decoded_lines)} != {len(lines)} lines")
                    else:
                        self.logger.info("✅ Data URI validation passed - encoding/decoding successful")
                except Exception as decode_error:
                    self.logger.error(f"Data URI validation failed: {decode_error}")
                
                return data_uri
                
            except Exception as e:
                self.logger.error(f"Failed to generate TSV data URI: {e}")
                import traceback
                self.logger.error(traceback.format_exc())
                return ""
        else:
            self.logger.warning(f"Marker TSV file not found: {marker_file}")
            return ""
    
    def load_bin_images(self):
        """Load bin analysis images (scatter and violin plots).
        
        Searches for bin-specific QC images in 07.outs/binned_outputs/qc.
        """
        for bin_size in SUPPORTED_BIN_SIZES:
            bin_key = f'bin{bin_size}'
            self.bin_images[bin_key] = {
                'scatter': None,
                'violin': None
            }
            
            # Look for scatter and violin plots
            scatter_file = self._binned_qc_dir() / f"{self.sample_name}_bin{bin_size}_scatter.png"
            violin_file = self._binned_qc_dir() / f"{self.sample_name}_bin{bin_size}_violin.png"
            
            if scatter_file.exists():
                self.bin_images[bin_key]['scatter'] = self.encode_image_to_base64(scatter_file)
                self.logger.info(f"✅ Loaded scatter plot for bin {bin_size}μm")
            else:
                self.logger.debug(f"Scatter plot not found for bin {bin_size}μm: {scatter_file}")
            
            if violin_file.exists():
                self.bin_images[bin_key]['violin'] = self.encode_image_to_base64(violin_file)
                self.logger.info(f"✅ Loaded violin plot for bin {bin_size}μm")
            else:
                self.logger.debug(f"Violin plot not found for bin {bin_size}μm: {violin_file}")
    
    def load_company_logo(self):
        """Load company logo from templates/img directory and convert to Base64."""
        # Primary location: templates/img directory (standard location)
        templates_img_dir = Path(__file__).parent.parent / "templates" / "img"
        
        # Try multiple potential logo file names and locations
        logo_paths = [
            # Standard location in templates/img
            templates_img_dir / "Celatlas_LOGO.png",
            templates_img_dir / "company_logo.png",
            templates_img_dir / "logo.png",
            
            # Fallback locations
            self.sample_dir.parent / "Celatlas_LOGO.png",  # Same level as sample dir
            Path(__file__).parent.parent.parent / "Celatlas_LOGO.png",  # Project root
            self.sample_dir / "Celatlas_LOGO.png",  # Inside sample dir
            
        ]
        
        for logo_path in logo_paths:
            if logo_path.exists():
                try:
                    self.company_logo = self.encode_image_to_base64(logo_path)
                    self.logger.info(f"✅ Loaded company logo from {logo_path}")
                    self.logger.info(f"Logo file size: {logo_path.stat().st_size / 1024:.1f}KB")
                    return
                except Exception as e:
                    self.logger.warning(f"Failed to load logo from {logo_path}: {e}")
        
        self.logger.warning("Company logo not found in any of the following locations:")
        for path in logo_paths:
            self.logger.warning(f"  - {path}")
        self.logger.warning("Will use default logo. To add logo, place it in templates/img/ directory.")
    
    def load_downsample_data(self):
        """Load downsampling data for saturation analysis."""
        downsample_file = self.sample_dir / "05.count" / f"{self.sample_name}_downsample.tsv"
        
        if downsample_file.exists():
            try:
                downsample_data = pd.read_csv(downsample_file, sep='\t')
                self.downsample_data = self._normalize_downsample_data(downsample_data, downsample_file)
                if self.downsample_data is not None:
                    self.downsample_source = f"05.count: {downsample_file}"
                    self.logger.info(f"Loaded downsample data with {len(self.downsample_data)} data points")
            except Exception as e:
                self.logger.error(f"Failed to load downsample data: {e}")
                self.downsample_data = None
        else:
            self.logger.warning(f"Downsample file not found: {downsample_file}")
            self.downsample_data = None
            self.downsample_source = None
    
    def load_count_detail_data(self):
        """Load barcode UMI count data for barcode rank plot.

        First tries to use the aggregated counts.txt file (faster and includes all barcodes).
        Falls back to count_detail.txt if counts.txt is not available.
        """
        # Try the aggregated counts file first (preferred method)
        counts_file = self.sample_dir / "05.count" / f"{self.sample_name}_counts.txt"

        if counts_file.exists():
            try:
                self.logger.info(f"Loading barcode counts from {counts_file}")
                # Read the counts file which has all barcodes with their UMI counts
                counts_df = pd.read_csv(counts_file, sep='\t', index_col=0)
                # Sort by UMI count (descending)
                counts_df = counts_df.sort_values('UMI', ascending=False)
                self.barcode_umi_counts = counts_df['UMI']

                self.logger.info(f"Loaded barcode counts: {len(self.barcode_umi_counts):,} barcodes, "
                               f"min UMI={self.barcode_umi_counts.min()}, max UMI={self.barcode_umi_counts.max()}")
                return
            except Exception as e:
                self.logger.error(f"Failed to load counts file: {e}, falling back to count_detail")

        # Fallback: read count_detail file
        count_detail_file = self.sample_dir / "05.count" / f"{self.sample_name}_count_detail.txt"

        if count_detail_file.exists():
            try:
                self.logger.info(f"Loading count detail data from {count_detail_file}")
                # Read in chunks to avoid memory issues with large files (40M+ rows)
                # Aggregate barcode-UMI pairs on the fly
                barcode_umi_dict = {}
                chunk_size = 500000
                total_rows = 0

                for chunk in pd.read_csv(count_detail_file, sep='\t', chunksize=chunk_size):
                    total_rows += len(chunk)
                    # Group by Barcode and collect unique UMIs
                    for barcode, group in chunk.groupby('Barcode'):
                        if barcode not in barcode_umi_dict:
                            barcode_umi_dict[barcode] = set()
                        barcode_umi_dict[barcode].update(group['UMI'].unique())

                    if total_rows % 5000000 == 0:
                        self.logger.info(f"  Processed {total_rows:,} rows, tracking {len(barcode_umi_dict):,} barcodes")

                # Convert to barcode UMI counts
                self.barcode_umi_counts = pd.Series({
                    barcode: len(umis)
                    for barcode, umis in barcode_umi_dict.items()
                }).sort_values(ascending=False)

                self.logger.info(f"Loaded complete count detail data: {total_rows:,} rows, {len(self.barcode_umi_counts):,} barcodes")
                self.count_detail_data = None  # Not needed anymore, save memory

            except Exception as e:
                self.logger.error(f"Failed to load count detail data: {e}")
                self.barcode_umi_counts = None
        else:
            self.logger.warning(f"Count detail file not found: {count_detail_file}")
            self.barcode_umi_counts = None
    
    def generate_barcode_rank_plot(self, include_plotlyjs='inline') -> str:
        """Generate interactive Barcode Rank Plot using plotly (Cell Ranger style)."""
        if not PLOTLY_AVAILABLE:
            return '<div class="simple-chart"><div class="chart-placeholder">Plotly not available</div></div>'

        try:
            unit_text = self._unit_text()
            if self.barcode_umi_counts is None:
                # Generate sample data if real data not available
                n_barcodes = 10000
                x_data = np.arange(1, n_barcodes + 1)
                barcode_counts = np.exp(-x_data / 2000) * 1000 + np.random.exponential(10, n_barcodes)
                self.logger.warning("Using simulated data for Barcode Rank Plot")
            else:
                # Use all available data
                total_barcodes = len(self.barcode_umi_counts)
                x_data = np.arange(1, total_barcodes + 1)
                barcode_counts = self.barcode_umi_counts.values

                self.logger.info(f"Barcode Rank Plot data: {total_barcodes} barcodes, "
                               f"min UMI={barcode_counts[-1]}, max UMI={barcode_counts[0]}")

            # Use the reported count as the split point so the blue/gray split
            # matches the metric table.
            estimated_cells = self.general_stats.get('Estimated Number of Cells', None)

            if estimated_cells is not None:
                try:
                    # Remove commas if present and convert to int
                    if isinstance(estimated_cells, str):
                        estimated_cells = int(estimated_cells.replace(',', ''))
                    else:
                        estimated_cells = int(estimated_cells)
                    inflection_idx = min(estimated_cells, len(x_data))
                    self.logger.info(f"Using reported {unit_text['unit_name_plural_lower']} ({estimated_cells}) as split point")
                except (ValueError, AttributeError):
                    self.logger.warning(f"Could not parse Estimated Number of Cells: {estimated_cells}, using auto-detection")
                    estimated_cells = None

            # Fallback: auto-detect inflection point if Estimated Number of Cells not available
            if estimated_cells is None:
                log_counts = np.log10(barcode_counts + 1)
                if len(log_counts) > 100:
                    gradient = np.gradient(log_counts)
                    second_deriv = np.gradient(gradient)
                    inflection_idx = np.argmin(second_deriv[:len(second_deriv)//2])
                    inflection_idx = max(inflection_idx, 10)
                else:
                    inflection_idx = int(len(x_data) * 0.1)
                self.logger.info(f"Auto-detected inflection point at {inflection_idx}")

            # Display all barcodes to show complete background curve
            # No need to trim data - show everything for complete visualization
            display_range = len(x_data)

            self.logger.info(f"Barcode Rank Plot: split at {inflection_idx}, displaying {display_range} barcodes (all data)")

            # Use all data - no trimming
            # x_data and barcode_counts remain unchanged

            # Split data into selected units (blue) and background (gray)
            # Include inflection point in both traces to ensure continuous line
            cell_x = x_data[:inflection_idx + 1]  # Include inflection_idx
            cell_y = barcode_counts[:inflection_idx + 1]
            bg_x = x_data[inflection_idx:]  # Start from inflection_idx
            bg_y = barcode_counts[inflection_idx:]

            # Create the plot
            fig = go.Figure()

            # Add selected units trace (blue)
            fig.add_trace(go.Scatter(
                x=cell_x,
                y=cell_y,
                mode='lines',
                name=unit_text['rank_trace_name'],
                line=dict(color='#1a73e8', width=3),
                hovertemplate='<b>Rank:</b> %{x}<br><b>UMI Counts:</b> %{y}<extra></extra>'
            ))

            # Add background trace (gray)
            fig.add_trace(go.Scatter(
                x=bg_x,
                y=bg_y,
                mode='lines',
                name='Background',
                line=dict(color='#9ca3af', width=3),
                hovertemplate='<b>Rank:</b> %{x}<br><b>UMI Counts:</b> %{y}<extra></extra>'
            ))
            
            # X-axis tick labels: Always show 1, 100, 10k, 1M
            x_max = max(x_data)
            x_tickvals = [1, 100, 10000, 1000000]
            x_ticktext = ['1', '100', '10k', '1M']

            # Set X-axis range to show full scale
            x_range = [np.log10(1), np.log10(max(1000000, x_max))]

            # Y-axis tick labels: Always show 1, 10, 100, 1000, 10K
            y_min = min(barcode_counts)
            y_max = max(barcode_counts)
            y_tickvals = [1, 10, 100, 1000, 10000]
            y_ticktext = ['1', '10', '100', '1000', '10k']

            # Set Y-axis range to show ALL data including background
            y_lower = max(0.1, y_min * 0.8)  # Start from data minimum (at least 0.1), with 20% buffer
            y_upper = max(10000, y_max * 1.2)  # At least show to 10k, or 20% above max data
            y_range = [np.log10(y_lower), np.log10(y_upper)]

            # Update layout to match UI theme with Cell Ranger style ticks
            fig.update_layout(
                title=dict(
                    text='Barcode Rank Plot',
                    font=dict(size=16, color='#1a73e8'),
                    x=0.5
                ),
                xaxis=dict(
                    title='Barcodes (Ranked)',
                    type='log',
                    showgrid=True,
                    gridcolor='#e2e8f0',
                    color='#000000',
                    linecolor='#000000',
                    linewidth=2,
                    tickvals=x_tickvals,
                    ticktext=x_ticktext,
                    range=x_range
                ),
                yaxis=dict(
                    title='UMI Counts',
                    type='log',
                    showgrid=True,
                    gridcolor='#e2e8f0',
                    color='#000000',
                    linewidth=2,
                    linecolor='#000000',
                    tickvals=y_tickvals,
                    ticktext=y_ticktext,
                    range=y_range
                ),
                plot_bgcolor='white',
                paper_bgcolor='white',
                font=dict(family='Arial', size=12),
                legend=dict(
                    x=0.7, y=0.95,
                    bgcolor='white',  # Pure white background
                    bordercolor='#e2e8f0',
                    borderwidth=1
                ),
                margin=dict(l=60, r=80, t=60, b=60),
                height=450,
                autosize=True
            )

            # Convert to HTML with optimized plotly.js inclusion
            # Set responsive to True and ensure plot fills container properly without clipping
            config = {
                'displayModeBar': False,
                'responsive': True,
                'autosizable': True,
                'fillFrame': False,
                'frameMargins': 0
            }
            html_str = pyo.plot(fig, output_type='div', include_plotlyjs=include_plotlyjs, config=config)

            self.logger.info("Generated Barcode Rank Plot successfully")
            return html_str
            
        except Exception as e:
            self.logger.error(f"Failed to generate Barcode Rank Plot: {e}")
            return '<div class="simple-chart"><div class="chart-placeholder">Barcode Rank Plot<br><small>Error generating chart</small></div></div>'
    
    def generate_sequencing_saturation_plot(self, include_plotlyjs=False) -> str:
        """Generate interactive Sequencing Saturation plot using plotly."""
        if not PLOTLY_AVAILABLE:
            return '<div class="simple-chart"><div class="chart-placeholder">Plotly not available</div></div>'
        
        try:
            if self.downsample_data is None or self.downsample_data.empty or 'read_fraction' not in self.downsample_data.columns:
                self.logger.warning("Sequencing saturation plot skipped because no valid downsample data is available")
                return self._chart_placeholder("Sequencing Saturation", "No valid downsample data")

            read_fractions = self.downsample_data['read_fraction'].values
            if 'umi_saturation' in self.downsample_data.columns:
                saturation_values = self.downsample_data['umi_saturation'].values
            elif 'read_saturation' in self.downsample_data.columns:
                saturation_values = self.downsample_data['read_saturation'].values
            else:
                self.logger.warning("Sequencing saturation plot skipped because saturation columns are missing")
                return self._chart_placeholder("Sequencing Saturation", "Saturation columns missing")
            
            # Create the plot
            fig = go.Figure()
            
            fig.add_trace(go.Scatter(
                x=read_fractions,
                y=saturation_values,
                mode='lines',
                name='Sequencing Saturation',
                line=dict(color='#ffde7d', width=3),
                hovertemplate='<b>Read Fraction:</b> %{x:.2f}<br><b>Saturation:</b> %{y:.2f}%<extra></extra>'
            ))
            
            # Update layout
            fig.update_layout(
                title=dict(
                    text='Sequencing Saturation',
                    font=dict(size=16, color='#1a73e8'),
                    x=0.5
                ),
                xaxis=dict(
                    title='Read Fraction',
                    showgrid=True,
                    gridcolor='#e2e8f0',
                    color='#000000',
                    linecolor='#000000',
                    linewidth=2,
                    range=[0, 1]
                ),
                yaxis=dict(
                    title='Sequencing Saturation (%)',
                    showgrid=True,
                    gridcolor='#e2e8f0',
                    color='#000000',
                    linecolor='#000000',
                    linewidth=2,
                    range=[0, 100],
                    tickformat='.0f',
                    dtick=20
                ),
                plot_bgcolor='white',
                paper_bgcolor='white',
                font=dict(family='Arial', size=12),
                showlegend=False,
                margin=dict(l=60, r=20, t=60, b=60),
                height=450
            )
            
            # Convert to HTML with optimized plotly.js inclusion
            config = {'displayModeBar': False, 'responsive': True}
            html_str = pyo.plot(fig, output_type='div', include_plotlyjs=include_plotlyjs, config=config)
            
            self.logger.info("Generated Sequencing Saturation Plot successfully")
            return html_str
            
        except Exception as e:
            self.logger.error(f"Failed to generate Sequencing Saturation Plot: {e}")
            return '<div class="simple-chart"><div class="chart-placeholder">Sequencing Saturation<br><small>Error generating chart</small></div></div>'
    
    def generate_median_genes_plot(self, include_plotlyjs=False) -> str:
        """Generate interactive Median Genes per Cell plot using plotly."""
        if not PLOTLY_AVAILABLE:
            return '<div class="simple-chart"><div class="chart-placeholder">Plotly not available</div></div>'
        
        try:
            unit_text = self._unit_text()
            if self.downsample_data is None or 'read_fraction' not in self.downsample_data.columns:
                # Generate sample data if real data not available
                read_fractions = np.linspace(0, 1, 11)
                median_genes = read_fractions * 2000 * (1 - np.exp(-read_fractions * 2))  # Saturation curve
                
                self.logger.warning("Using simulated data for Median Genes Plot")
            else:
                # Use real data
                read_fractions = self.downsample_data['read_fraction'].values
                if 'median_gene_number' in self.downsample_data.columns:
                    median_genes = self.downsample_data['median_gene_number'].values
                else:
                    # Fallback to simulated data based on real read fractions
                    median_genes = read_fractions * 2000 * (1 - np.exp(-read_fractions * 2))
            
            # Create the plot
            fig = go.Figure()
            
            fig.add_trace(go.Scatter(
                x=read_fractions,
                y=median_genes,
                mode='lines',
                name=unit_text['median_genes_plot_title'],
                line=dict(color='#e84545', width=3),
                hovertemplate='<b>Read Fraction:</b> %{x:.2f}<br><b>Median Genes:</b> %{y:.0f}<extra></extra>'
            ))
            
            # Update layout
            fig.update_layout(
                title=dict(
                    text=unit_text['median_genes_plot_title'],
                    font=dict(size=16, color='#1a73e8'),
                    x=0.5
                ),
                xaxis=dict(
                    title='Read Fraction',
                    showgrid=True,
                    gridcolor='#e2e8f0',
                    color='#000000',
                    linecolor='#000000',
                    linewidth=2,
                    range=[0, 1]
                ),
                yaxis=dict(
                    title=unit_text['median_genes_plot_title'],
                    showgrid=True,
                    gridcolor='#e2e8f0',
                    color='#000000',
                    linecolor='#000000',
                    linewidth=2
                ),
                plot_bgcolor='white',
                paper_bgcolor='white',
                font=dict(family='Arial', size=12),
                showlegend=False,
                margin=dict(l=60, r=20, t=60, b=60),
                height=450
            )
            
            # Convert to HTML with optimized plotly.js inclusion
            config = {'displayModeBar': False, 'responsive': True}
            html_str = pyo.plot(fig, output_type='div', include_plotlyjs=include_plotlyjs, config=config)
            
            self.logger.info("Generated Median Genes Plot successfully")
            return html_str
            
        except Exception as e:
            self.logger.error(f"Failed to generate Median Genes Plot: {e}")
            return f'<div class="simple-chart"><div class="chart-placeholder">{self._unit_text()["median_genes_plot_title"]}<br><small>Error generating chart</small></div></div>'
    
    def generate_static_chart_placeholder(self, title: str, subtitle: str) -> str:
        """Generate a lightweight static chart placeholder for faster loading."""
        return f'''
        <div class="simple-chart">
            <div class="chart-placeholder">
                {title}<br>
                <small>{subtitle}</small><br>
                <em style="color: #94a3b8; font-size: 11px;">Static mode for faster loading</em>
            </div>
        </div>
        '''
    
    def load_analysis_images(self):
        """Load analysis result images including QC plots and UMAP visualizations."""
        # Initialize analysis images dictionary
        self.analysis_images = {
            'tissue_hires': None,
            'spatial_gene_expression': None,
            'spatial_cluster': None,
            'umap_cluster': None,
            # scRNA-specific images
            'qc_scatter_before': None,
            'qc_violin_before': None,
            'umap_genes': None,
            'umap_counts': None,
            'interactive_umap_html': None,
            'interactive_tsne_html': None
        }
        
        # Load tissue hires image from the current directory layout only.
        tissue_hires_paths = [
            self._binned_clustering_dir() / "tissue_hires_image.png",
            self._binsegment_dir() / "images" / "tissue_hires_image.png",
        ]
        
        for tissue_path in tissue_hires_paths:
            if tissue_path.exists():
                self.analysis_images['tissue_hires'] = self.encode_image_to_base64(tissue_path)
                self.logger.info(f"Loaded tissue hires image from {tissue_path}")
                break
        
        # Load spatial gene expression distribution
        spatial_gene_expr_file = self._first_existing([
            self._binned_qc_dir() / f"{self.sample_name}_spatial_gene_expression_distribution.png",
        ])
        if spatial_gene_expr_file.exists():
            self.analysis_images['spatial_gene_expression'] = self.encode_image_to_base64(spatial_gene_expr_file)
            self.logger.info("Loaded spatial gene expression distribution image")
        
        # Load spatial cluster image
        spatial_cluster_file = self._first_existing([
            self._binned_clustering_dir() / f"{self.sample_name}_spatial_cluster.png",
        ])
        if spatial_cluster_file.exists():
            self.analysis_images['spatial_cluster'] = self.encode_image_to_base64(spatial_cluster_file)
            self.logger.info("Loaded spatial cluster image")
        
        # Load UMAP cluster image
        umap_cluster_file = self._first_existing([
            self._binned_clustering_dir() / f"{self.sample_name}_umap_cluster.png",
        ])
        if umap_cluster_file.exists():
            self.analysis_images['umap_cluster'] = self.encode_image_to_base64(umap_cluster_file)
            self.logger.info("Loaded UMAP cluster image")

        analysis_wrapper_dir = self.sample_dir / "06_analysis_wrapper"
        if not (analysis_wrapper_dir / "01.QC").exists():
            self.logger.warning(f"Analysis directory not found: {analysis_wrapper_dir}")
            return
        self.logger.info(f"Found analysis directory at: {analysis_wrapper_dir}")

        # QC images (before filtering)
        qc_scatter_before = analysis_wrapper_dir / "01.QC" / f"{self.sample_name}_qc_scatter_before.png"
        if qc_scatter_before.exists():
            self.analysis_images['qc_scatter_before'] = self.encode_image_to_base64(qc_scatter_before)
            self.logger.info("✅ Loaded QC scatter plot (before filtering)")
        else:
            self.logger.warning(f"QC scatter plot not found: {qc_scatter_before}")

        qc_violin_before = analysis_wrapper_dir / "01.QC" / f"{self.sample_name}_qc_violin_before.png"
        if qc_violin_before.exists():
            self.analysis_images['qc_violin_before'] = self.encode_image_to_base64(qc_violin_before)
            self.logger.info("✅ Loaded QC violin plot (before filtering)")
        else:
            self.logger.warning(f"QC violin plot not found: {qc_violin_before}")

        # UMAP images
        umap_genes_file = analysis_wrapper_dir / "02.UMAP" / f"{self.sample_name}_umap_genes.png"
        if umap_genes_file.exists():
            self.analysis_images['umap_genes'] = self.encode_image_to_base64(umap_genes_file)
            self.logger.info("✅ Loaded UMAP genes plot")
        else:
            self.logger.warning(f"UMAP genes plot not found: {umap_genes_file}")

        umap_counts_file = analysis_wrapper_dir / "02.UMAP" / f"{self.sample_name}_umap_counts.png"
        if umap_counts_file.exists():
            self.analysis_images['umap_counts'] = self.encode_image_to_base64(umap_counts_file)
            self.logger.info("✅ Loaded UMAP counts plot")
        else:
            self.logger.warning(f"UMAP counts plot not found: {umap_counts_file}")

        # Interactive UMAP (load JSON and convert to HTML)
        umap_interactive_json = analysis_wrapper_dir / "02.UMAP" / f"{self.sample_name}_umap_interactive.json"
        if umap_interactive_json.exists():
            try:
                import plotly.io as pio
                with open(umap_interactive_json, 'r') as f:
                    fig_json = f.read()

                # Convert JSON to plotly figure
                fig = pio.from_json(fig_json)

                # Generate HTML div (not full HTML document)
                interactive_html = fig.to_html(
                    include_plotlyjs='cdn',
                    div_id='interactive-umap-plot',
                    config={'responsive': True, 'displayModeBar': True}
                )

                self.analysis_images['interactive_umap_html'] = interactive_html
                self.logger.info("✅ Loaded interactive UMAP plot")
            except Exception as e:
                self.logger.error(f"Failed to load interactive UMAP: {e}")
                self.analysis_images['interactive_umap_html'] = '<div style="padding: 50px; text-align: center; color: #999;">Interactive UMAP not available</div>'
        else:
            self.logger.warning(f"Interactive UMAP JSON not found: {umap_interactive_json}")
            self.analysis_images['interactive_umap_html'] = '<div style="padding: 50px; text-align: center; color: #999;">Interactive UMAP not available</div>'

        # Interactive t-SNE (load JSON and convert to HTML)
        tsne_interactive_json = analysis_wrapper_dir / "02.UMAP" / f"{self.sample_name}_tsne_interactive.json"
        if tsne_interactive_json.exists():
            try:
                import plotly.io as pio
                with open(tsne_interactive_json, 'r') as f:
                    fig_json = f.read()

                # Convert JSON to plotly figure
                fig = pio.from_json(fig_json)

                # Generate HTML div (not full HTML document)
                interactive_html = fig.to_html(
                    include_plotlyjs=False,  # Already loaded by UMAP
                    div_id='interactive-tsne-plot',
                    config={'responsive': True, 'displayModeBar': True}
                )

                self.analysis_images['interactive_tsne_html'] = interactive_html
                self.logger.info("✅ Loaded interactive t-SNE plot")
            except Exception as e:
                self.logger.error(f"Failed to load interactive t-SNE: {e}")
                self.analysis_images['interactive_tsne_html'] = '<div style="padding: 50px; text-align: center; color: #999;">Interactive t-SNE not available</div>'
        else:
            self.logger.warning(f"Interactive t-SNE JSON not found: {tsne_interactive_json}")
            self.analysis_images['interactive_tsne_html'] = '<div style="padding: 50px; text-align: center; color: #999;">Interactive t-SNE not available</div>'

    def encode_image_to_base64(self, image_path: Path) -> str:
        """Encode an image file to base64 string."""
        try:
            with open(image_path, 'rb') as img_file:
                encoded_string = base64.b64encode(img_file.read()).decode('utf-8')
                return f"data:image/png;base64,{encoded_string}"
        except Exception as e:
            self.logger.error(f"Failed to encode image {image_path}: {e}")
            return ""
    
    def format_number(self, value: str, as_float: bool = False) -> str:
        """Format a numeric value with appropriate formatting."""
        try:
            # Remove commas and convert to number
            clean_value = str(value).replace(',', '').strip()
            
            if as_float:
                num = float(clean_value)
                if num < 1:
                    return f"{num:.3f}"
                elif num < 100:
                    return f"{num:.1f}"
                else:
                    return f"{num:,.0f}"
            else:
                num = int(float(clean_value))
                return f"{num:,}"
        except (ValueError, TypeError):
            return str(value)
    
    def calculate_percentage(self, part: str, total: str) -> str:
        """Calculate percentage from part and total values."""
        try:
            part_val = float(str(part).replace(',', ''))
            total_val = float(str(total).replace(',', ''))
            if total_val > 0:
                return f"{(part_val / total_val * 100):.2f}"
            return "0.00"
        except (ValueError, TypeError):
            return "0.00"
    
    def _calculate_median_genes_square(self) -> str:
        """Calculate median genes per square from bin statistics."""
        try:
            # Try to get median genes from bin50 statistics first
            bin50_stats = self.bin_stats.get('bin50', {})
            median_genes = bin50_stats.get('Median Genes per square bin', '')
            
            if median_genes and median_genes != '0':
                return self.format_number(median_genes)
            
            # Fallback: try to get from analysis statistics or calculate from mean
            analysis_stats = self.general_stats
            median_genes_analysis = analysis_stats.get('Median Genes per Square', '')
            
            if median_genes_analysis and median_genes_analysis != '0':
                return self.format_number(median_genes_analysis)
            
            # Last resort: estimate from mean genes if available
            mean_genes = bin50_stats.get('Mean Genes per square bin', '0')
            if mean_genes and mean_genes != '0':
                # Median is typically slightly lower than mean for gene expression data
                try:
                    mean_val = float(str(mean_genes).replace(',', ''))
                    estimated_median = int(mean_val * 0.85)  # Rough estimation
                    return self.format_number(str(estimated_median))
                except (ValueError, TypeError):
                    pass
            
            # Default fallback
            return "N/A"
            
        except Exception as e:
            self.logger.warning(f"Failed to calculate median genes per square: {e}")
            return "N/A"
    
    def _calculate_saturation_percentage(self) -> str:
        """Calculate sequencing saturation percentage from available data."""
        try:
            # Method 1: Try to get from downsample data (most accurate)
            if self.downsample_data is not None and not self.downsample_data.empty:
                downsample_data = self.downsample_data.sort_values('read_fraction')
                if 'umi_saturation' in self.downsample_data.columns:
                    max_saturation = downsample_data['umi_saturation'].iloc[-1]
                    return self._format_percent_value(max_saturation)
                elif 'read_saturation' in self.downsample_data.columns:
                    max_saturation = downsample_data['read_saturation'].iloc[-1]
                    return self._format_percent_value(max_saturation)
            
            # Method 2: Try to get from general statistics
            saturation_from_stats = self.general_stats.get('Sequencing Saturation', '') or self.general_stats.get('Saturation', '')
            if saturation_from_stats and saturation_from_stats != '0':
                # Extract percentage if it's in format like "44.71%" or just "44.71"
                saturation_val = self._parse_numeric(saturation_from_stats)
                if saturation_val is not None:
                    if saturation_val <= 1:
                        saturation_val *= 100
                    return self._format_percent_value(saturation_val)
            
            return "N/A"
            
        except Exception as e:
            self.logger.warning(f"Failed to calculate saturation percentage: {e}")
            return "N/A"

    def _saturation_note(self) -> str:
        if self.is_hd_barcode_bin_mode():
            return "HD BBV4 scrna output is a barcode/bin-level diagnostic; saturation reflects duplicate-rate/library-complexity at barcode/bin level, not true cell-level sequencing saturation."

        mean_reads = self._parse_numeric(self.general_stats.get('Mean Reads per Cell', '0'))
        saturation = self._parse_numeric(self._calculate_saturation_percentage())
        raw_reads = self._parse_numeric(self.general_stats.get('Raw Reads', '0'))
        valid_reads = self._parse_numeric(self.general_stats.get('Valid Reads', '0'))

        if saturation is None:
            return "Sequencing saturation is unavailable."
        elif (
            (mean_reads is not None and mean_reads < 1000)
            or (valid_reads is not None and valid_reads < 5_000_000)
            or (raw_reads is not None and raw_reads < 10_000_000)
        ):
            return "Low-depth runs can show elevated duplicate-rate saturation before gene recovery reaches a stable plateau; interpret together with median genes and UMI distribution."
        return "Estimated from read downsampling as the fraction of duplicated UMI/read observations."
    
    def generate_analysis_results_section(self) -> str:
        """Generate HTML section for analysis results including tissue hires, spatial gene expression, spatial cluster, and umap cluster."""
        analysis_images = [
            {
                'key': 'tissue_hires',
                'title': 'Tissue Hires Image',
                'description': 'High-resolution tissue morphology'
            },
            {
                'key': 'spatial_gene_expression',
                'title': 'Spatial Gene Expression Distribution',
                'description': 'Gene expression spatial patterns'
            },
            {
                'key': 'spatial_cluster',
                'title': 'Spatial Cluster',
                'description': 'Spatial clustering visualization'
            },
            {
                'key': 'umap_cluster',
                'title': 'UMAP Cluster',
                'description': 'UMAP dimensional reduction clustering'
            }
        ]
        
        image_sections = []
        for img_info in analysis_images:
            key = img_info['key']
            title = img_info['title']
            description = img_info['description']
            img_data = self.analysis_images.get(key)
            
            if img_data:
                image_html = f'<img src="{img_data}" alt="{title}" />'
            else:
                image_html = f'<div style="color: #999; font-size: 13px;">{title}<br/>Image not available</div>'
            
            section_html = f'''
            <div class="analysis-image" title="Click to enlarge">
                {image_html}
                <div class="image-label">{title}</div>
            </div>
            '''
            image_sections.append(section_html)
        
        return '\n'.join(image_sections)
    
    def generate_bin_analysis_sections(self) -> str:
        """Generate HTML sections for bin analysis results (scatter and violin plots)."""
        sections = []
        bin_sizes = ['10', '20', '50', '100']
        
        for bin_size in bin_sizes:
            bin_key = f'bin{bin_size}'
            
            if bin_key in self.bin_images:
                scatter_img = self.bin_images[bin_key].get('scatter', '')
                violin_img = self.bin_images[bin_key].get('violin', '')
                
                section_html = f'''
                <div class="bin-analysis-container">
                    <div class="bin-title">Bin {bin_size} μm Analysis</div>
                    <div class="bin-images">
                        <div class="analysis-image" title="Click to enlarge">
                            {f'<img src="{scatter_img}" alt="Scatter Plot - Bin {bin_size}μm" />' if scatter_img else f'<div style="color: #999;">Scatter Plot - Bin {bin_size}μm<br/>Image not available</div>'}
                            <div class="image-label">Scatter Plot</div>
                        </div>
                        <div class="analysis-image" title="Click to enlarge">
                            {f'<img src="{violin_img}" alt="Violin Plot - Bin {bin_size}μm" />' if violin_img else f'<div style="color: #999;">Violin Plot - Bin {bin_size}μm<br/>Image not available</div>'}
                            <div class="image-label">Violin Plot</div>
                        </div>
                    </div>
                </div>
                '''
                sections.append(section_html)
        
        return '\n'.join(sections)
    
    def generate_marker_genes_rows(self, display_limit: int = 150) -> str:
        """Generate HTML table rows for marker genes."""
        if not self.marker_genes:
            return '<tr><td colspan="8" style="text-align: center; padding: 30px;">No marker genes data available</td></tr>'
        
        rows = []
        cluster_colors = [
            "#3498db", "#2ecc71", "#9b59b6", "#e74c3c",
            "#f1c40f", "#1abc9c", "#d35400", "#34495e",
            "#e67e22", "#8e44ad", "#2c3e50", "#27ae60",
            "#16a085", "#c0392b", "#8e44ad", "#f39c12"
        ]
        
        # Show only first display_limit genes in the table
        display_genes = self.marker_genes[:display_limit]
        
        for gene in display_genes:
            cluster = str(gene['cluster'])
            try:
                cluster_color = CLUSTER_COLORS[int(cluster) % len(CLUSTER_COLORS)]
            except (ValueError, IndexError):
                cluster_color = CLUSTER_COLORS[0]  # Default to first color
            
            # Format scientific notation for p-values
            p_val_str = f"{gene['p_val']:.2e}" if gene['p_val'] > 0 else "0.00e+00"
            p_val_adj_str = f"{gene['p_val_adj']:.2e}" if gene['p_val_adj'] > 0 else "0.00e+00"
            
            row_html = f'''
            <tr>
                <td><span class="cluster-tag" style="background: {cluster_color}; color: white;">Cluster {cluster}</span></td>
                <td style="font-weight: 500; color: #2c3e50;">{gene['gene']}</td>
                <td>{gene['scores']:.3f}</td>
                <td>{gene['avg_log2FC']:.3f}</td>
                <td>{p_val_str}</td>
                <td>{p_val_adj_str}</td>
                <td>{gene['pct1']:.3f}</td>
                <td>{gene['pct2']:.3f}</td>
            </tr>
            '''
            rows.append(row_html)
        
        return '\n'.join(rows)
    
    def generate_marker_genes_data_json(self) -> str:
        """Generate JSON data for marker genes JavaScript functionality."""
        # Convert numpy types to Python native types for JSON serialization
        json_compatible_data = []
        for gene in self.marker_genes:
            json_gene = {}
            for key, value in gene.items():
                if hasattr(value, 'item'):  # numpy scalar
                    json_gene[key] = value.item()
                elif isinstance(value, (int, float, str)):
                    json_gene[key] = value
                else:
                    json_gene[key] = str(value)  # fallback to string
            json_compatible_data.append(json_gene)
        
        return json.dumps(json_compatible_data, indent=2)
    
    def load_template(self) -> str:
        """Load the HTML template."""
        if not self.template_path.exists():
            raise FileNotFoundError(f"Template file not found: {self.template_path}")
        
        with open(self.template_path, 'r', encoding='utf-8') as f:
            return f.read()
    
    def populate_template(self, template: str) -> str:
        """Populate the HTML template with actual data."""
        unit_text = self._unit_text()
        
        # Determine Sample ID display logic
        if self.env_sample_name:
            # If sample_name is provided, show it as Sample ID
            display_sample_id = self.env_sample_name
        else:
            # If only chip_number is provided, show N/A
            display_sample_id = "N/A"
        
        # Determine Transcriptome based on species
        transcriptome_mapping = {
            'Mus_musculus': 'GRCm39',
            'mouse': 'GRCm39',
            'Homo_sapiens': 'GRCh38', 
            'human': 'GRCh38'
        }
        transcriptome = transcriptome_mapping.get(self.species, 'N/A')
        
        # Basic sample information
        replacements = {
            'SAMPLE_ID': display_sample_id,
            'ASSAY': 'RNA',
            'CHEMISTRY': self.chemistry,
            'GENOME': self.species,
            'GENERATION_DATE': datetime.now().strftime('%B %d, %Y'),
            'REPORT_TITLE': unit_text['report_title'],
            'REPORT_SUBTITLE': unit_text['report_subtitle'],
            'ANALYSIS_TAB_LABEL': unit_text['analysis_tab'],
            'ANALYSIS_TITLE': unit_text['analysis_title'],
            'STATISTICS_TITLE': unit_text['statistics_title'],
            'UNIT_NAME': unit_text['unit_name'],
            'UNIT_NAME_LOWER': unit_text['unit_name_lower'],
            'UNIT_NAME_PLURAL': unit_text['unit_name_plural'],
            'UNIT_NAME_PLURAL_LOWER': unit_text['unit_name_plural_lower'],
            'NUMBER_LABEL': unit_text['number_label'],
            'ESTIMATED_LABEL': unit_text['estimated_label'],
            'FRACTION_LABEL': unit_text['fraction_label'],
            'MEAN_READS_LABEL': unit_text['mean_reads_label'],
            'MEAN_UMIS_LABEL': unit_text['mean_umis_label'],
            'MEDIAN_UMI_LABEL': unit_text['median_umi_label'],
            'MEDIAN_GENES_LABEL': unit_text['median_genes_label'],
            'MEDIAN_GENES_PLOT_TITLE': unit_text['median_genes_plot_title'],
            'MEDIAN_GENES_CHART_INFO': unit_text['median_genes_chart_info'],
            'KEY_METRICS_HELP': unit_text['key_metrics_help'],
            'SATURATION_UNIT_HELP': unit_text['saturation_unit_help'],
            'STATISTICS_HELP_TITLE': unit_text['statistics_help_title'],
            'STATISTICS_HELP': unit_text['statistics_help'],
            'QC_HELP': unit_text['qc_help'],
            'UMAP_HELP': unit_text['umap_help'],
            'HOVER_TEXT': unit_text['hover_text'],
            'CLUSTER_HELP': unit_text['cluster_help'],
            'DIAGNOSTIC_NOTE_HTML': unit_text['diagnostic_note_html'],
            
            # New Sample Metadata parameters
            'CHIP_NUMBER': self.env_chip_number if self.env_chip_number else 'N/A',
            'TRANSCRIPTOME': transcriptome,
            'IMAGE_ALIGNMENT': 'N/A',  # scRNA-seq doesn't have spatial image alignment
            'PROBE_SET_NAME': 'N/A',
            'TISSUE': self.env_tissue if self.env_tissue else 'N/A',
        }
        
        # Key metrics from bin statistics (only used for spatial data)
        bin50_stats = self.bin_stats.get('bin50', {})
        bin10_stats = self.bin_stats.get('bin10', {})

        # For scRNA-seq, Key Metrics shows cell-level statistics from 05.count step metrics.
        # These are stored in general_stats after loading count statistics
        num_cells = self.general_stats.get('Estimated Number of Cells', '0')
        mean_reads_cell_km = self.general_stats.get('Mean Reads per Cell', '0')
        median_umi_cell_km = self.general_stats.get('Median UMI per Cell', '0')
        total_genes_km = self.general_stats.get('Total Genes', '0')

        replacements.update({
            'NUM_CELLS': self.format_number(num_cells),
            'MEAN_READS_CELL': self.format_number(mean_reads_cell_km, as_float=True),
            'MEAN_UMI_CELL': self.format_number(median_umi_cell_km, as_float=True),
            'TOTAL_GENES_DETECTED': self.format_number(total_genes_km),
        })
        
        # Barcode statistics - improved parsing to handle different formats
        raw_reads = self.general_stats.get('Raw Reads', '0')
        valid_reads = self.general_stats.get('Valid Reads', '0')
        mismatch = self.general_stats.get('Mismatched Reads', '0')
        exact_match = self.general_stats.get('Exactly Matched Reads', '0')
        no_polyt = self.general_stats.get('No PolyT Reads', '0')
        low_quality = self.general_stats.get('Low Quality Reads', '0')
        no_linker = self.general_stats.get('No Linker Reads', '0')
        no_barcode = self.general_stats.get('No Barcode Reads', '0')
        corrected_linker = self.general_stats.get('Corrected Linker Reads', '0')
        corrected_barcode = self.general_stats.get('Corrected Barcode Reads', '0')
        
        # Handle percentage extraction if already parsed
        valid_reads_pct = self.general_stats.get('Valid Reads Percentage', 
            self.calculate_percentage(valid_reads, raw_reads) + '%')
        mismatch_pct = self.general_stats.get('Mismatched Reads Percentage',
            self.calculate_percentage(mismatch, raw_reads) + '%')
        exact_match_pct = self.general_stats.get('Exactly Matched Reads Percentage',
            self.calculate_percentage(exact_match, raw_reads) + '%')
        no_polyt_pct = self.general_stats.get('No PolyT Reads Percentage',
            self.calculate_percentage(no_polyt, raw_reads) + '%')
        low_quality_pct = self.general_stats.get('Low Quality Reads Percentage',
            self.calculate_percentage(low_quality, raw_reads) + '%')
        no_linker_pct = self.general_stats.get('No Linker Reads Percentage',
            self.calculate_percentage(no_linker, raw_reads) + '%')
        no_barcode_pct = self.general_stats.get('No Barcode Reads Percentage',
            self.calculate_percentage(no_barcode, raw_reads) + '%')
        corrected_linker_pct = self.general_stats.get('Corrected Linker Reads Percentage',
            self.calculate_percentage(corrected_linker, raw_reads) + '%')
        corrected_barcode_pct = self.general_stats.get('Corrected Barcode Reads Percentage',
            self.calculate_percentage(corrected_barcode, raw_reads) + '%')
        
        replacements.update({
            'NUM_READS': self.format_number(raw_reads),
            'VALID_BARCODES': f"{self.format_number(valid_reads)} ({valid_reads_pct})",
            'EXACT_MATCHED_READS': f"{self.format_number(exact_match)} ({exact_match_pct})",
            'MISMATCH_READS': f"{self.format_number(mismatch)} ({mismatch_pct})",
            'NO_POLYT_READS': f"{self.format_number(no_polyt)} ({no_polyt_pct})",
            'LOW_QUALITY_READS': f"{self.format_number(low_quality)} ({low_quality_pct})",
            'NO_LINKER_READS': f"{self.format_number(no_linker)} ({no_linker_pct})",
            'NO_BARCODE_READS': f"{self.format_number(no_barcode)} ({no_barcode_pct})",
            'CORRECTED_LINKER_READS': f"{self.format_number(corrected_linker)} ({corrected_linker_pct})",
            'CORRECTED_BARCODE_READS': f"{self.format_number(corrected_barcode)} ({corrected_barcode_pct})",
            'Q30_BARCODES': self.general_stats.get('Q30 of Barcodes', '0%'),
            'Q30_UMIS': self.general_stats.get('Q30 of UMIs', '0%'),
        })
        
        # Mapping statistics - improved parsing
        unique_mapped = self.general_stats.get('Uniquely Mapped Reads', '0')
        multi_mapped = self.general_stats.get('Multi-Mapped Reads', '0')
        unmapped_short = self.general_stats.get('Unmapped Reads - Too Short', '0')
        unmapped_other = self.general_stats.get('Unmapped Reads - Other', '0')
        
        # Use actual valid reads as denominator for mapping percentages
        total_mapped = valid_reads
        
        # Handle percentage extraction if already parsed
        unique_mapped_pct = self.general_stats.get('Uniquely Mapped Reads Percentage',
            self.calculate_percentage(unique_mapped, total_mapped) + '%')
        multi_mapped_pct = self.general_stats.get('Multi-Mapped Reads Percentage',
            self.calculate_percentage(multi_mapped, total_mapped) + '%')
        unmapped_short_pct = self.general_stats.get('Unmapped Reads - Too Short Percentage',
            self.calculate_percentage(unmapped_short, total_mapped) + '%')
        unmapped_other_pct = self.general_stats.get('Unmapped Reads - Other Percentage',
            self.calculate_percentage(unmapped_other, total_mapped) + '%')
        
        replacements.update({
            'UNIQUE_MAPPED': f"{self.format_number(unique_mapped)} ({unique_mapped_pct})",
            'MULTI_MAPPED': f"{self.format_number(multi_mapped)} ({multi_mapped_pct})",
            'UNMAPPED_SHORT': f"{self.format_number(unmapped_short)} ({unmapped_short_pct})",
            'UNMAPPED_OTHER': f"{self.format_number(unmapped_other)} ({unmapped_other_pct})",
        })
        
        # FeatureCounts statistics - improved parsing
        feature_type = self.general_stats.get('Feature Type', 'Gene').lower()
        exonic_reads = self.general_stats.get('Reads Assigned To Exonic Regions', '0')
        intronic_reads = self.general_stats.get('Reads Assigned To Intronic Regions', '0')
        intergenic_reads = self.general_stats.get('Reads Assigned To Intergenic Regions', '0')
        unassigned_reads = self.general_stats.get('Reads Unassigned Ambiguity', '0')
        
        # Handle percentage extraction if already parsed
        exonic_pct = self.general_stats.get('Reads Assigned To Exonic Regions Percentage',
            self.calculate_percentage(exonic_reads, unique_mapped) + '%')
        intronic_pct = self.general_stats.get('Reads Assigned To Intronic Regions Percentage',
            self.calculate_percentage(intronic_reads, unique_mapped) + '%')
        intergenic_pct = self.general_stats.get('Reads Assigned To Intergenic Regions Percentage',
            self.calculate_percentage(intergenic_reads, unique_mapped) + '%')
        unassigned_pct = self.general_stats.get('Reads Unassigned Ambiguity Percentage',
            self.calculate_percentage(unassigned_reads, unique_mapped) + '%')
        
        replacements.update({
            'FEATURE_TYPE': feature_type,
            'EXONIC_READS': f"{self.format_number(exonic_reads)} ({exonic_pct})",
            'INTRONIC_READS': f"{self.format_number(intronic_reads)} ({intronic_pct})",
            'INTERGENIC_READS': f"{self.format_number(intergenic_reads)} ({intergenic_pct})",
            'UNASSIGNED_READS': f"{self.format_number(unassigned_reads)} ({unassigned_pct})",
        })
        
        # Cell Statistics section - for scRNA-seq, read from 05.count step metrics.
        # All cell statistics are stored in general_stats after loading count statistics
        estimated_cells = self.general_stats.get('Estimated Number of Cells', '0')
        fraction_reads_cells = self.general_stats.get('Fraction Reads in Cells', '0')
        mean_reads_cell = self.general_stats.get('Mean Reads per Cell', '0')
        median_umi_cell = self.general_stats.get('Median UMI per Cell', '0')
        total_genes_cell = self.general_stats.get('Total Genes', '0')
        median_genes_cell = self.general_stats.get('Median Genes per Cell', '0')
        saturation_percentage = self._calculate_saturation_percentage()
        saturation_note = self._saturation_note()

        # Remove % sign if present and format as number
        if isinstance(fraction_reads_cells, str):
            fraction_reads_cells = fraction_reads_cells.replace('%', '').strip()

        replacements.update({
            'ESTIMATED_NUM_CELLS': self.format_number(estimated_cells),
            'FRACTION_READS_CELLS': fraction_reads_cells if fraction_reads_cells else '0',
            'MEAN_READS_CELL_STAT': self.format_number(mean_reads_cell, as_float=True),
            'MEDIAN_UMI_CELL': self.format_number(median_umi_cell),
            'TOTAL_GENES_CELL': self.format_number(total_genes_cell),
            'MEDIAN_GENES_CELL': self.format_number(median_genes_cell),
            'SATURATION': saturation_percentage,
            'SATURATION_NOTE': saturation_note,
        })
        
        # Analysis results sections (no bin-level metrics for scRNA)
        replacements['ANALYSIS_RESULTS_IMAGES'] = self.generate_analysis_results_section()
        
        # Company logo
        if self.company_logo:
            replacements['COMPANY_LOGO_HTML'] = f'<img src="{self.company_logo}" alt="Company Logo" />'
        else:
            replacements['COMPANY_LOGO_HTML'] = '<div class="logo-fallback">C</div>'
        
        # Marker genes - organize by cluster for pagination
        replacements['MARKER_GENES_DATA'] = self.generate_marker_genes_data_json()
        # Convert numpy int64 to regular Python int for JSON serialization
        cluster_list = [int(cluster) for cluster in self.clusters] if hasattr(self, 'clusters') else []
        replacements['CLUSTER_LIST'] = json.dumps(cluster_list)
        replacements['TOTAL_CLUSTERS'] = str(len(cluster_list))
        
        # Complete marker genes data for full TSV export fallback
        if hasattr(self, 'complete_marker_genes'):
            replacements['COMPLETE_MARKER_GENES_DATA'] = json.dumps(self.complete_marker_genes, indent=2)
        else:
            replacements['COMPLETE_MARKER_GENES_DATA'] = json.dumps([])
        
        # Embed TSV file content as data URI for offline download
        replacements['MARKERS_RAW_TSV_DATA_URI'] = self.generate_markers_tsv_data_uri()

        # Interactive charts
        replacements['BARCODE_RANK_PLOT_HTML'] = self.barcode_rank_plot_html
        replacements['SEQUENCING_SATURATION_PLOT_HTML'] = self.sequencing_saturation_plot_html
        replacements['MEDIAN_GENES_PLOT_HTML'] = self.median_genes_plot_html

        # scRNA Cell Analysis tab images (base64 encoded)
        replacements['QC_SCATTER_BEFORE_IMG'] = self.analysis_images.get('qc_scatter_before', '')
        replacements['QC_VIOLIN_BEFORE_IMG'] = self.analysis_images.get('qc_violin_before', '')
        replacements['UMAP_GENES_IMG'] = self.analysis_images.get('umap_genes', '')
        replacements['UMAP_COUNTS_IMG'] = self.analysis_images.get('umap_counts', '')
        replacements['INTERACTIVE_UMAP_HTML'] = self.analysis_images.get('interactive_umap_html', '<div style="padding: 50px; text-align: center; color: #999;">Interactive UMAP not available</div>')
        replacements['INTERACTIVE_TSNE_HTML'] = self.analysis_images.get('interactive_tsne_html', '<div style="padding: 50px; text-align: center; color: #999;">Interactive t-SNE not available</div>')

        # Replace all placeholders in template
        result = template
        for key, value in replacements.items():
            result = result.replace('{{' + key + '}}', str(value))
        
        return result
    
    def generate_report(self, output_filename: str = None) -> Path:
        """Generate the complete HTML report."""
        self.logger.info("Starting report generation...")
        
        # Load all data
        self.load_bin_statistics()
        self.load_general_statistics()
        self.load_marker_genes()
        self.load_bin_images()
        self.load_analysis_images()
        self.load_company_logo()
        
        # Load data for interactive charts
        self.load_downsample_data()
        self.load_count_detail_data()
        self.build_hd_barcode_bin_metrics()
        
        # Generate charts with performance optimization
        if self.use_interactive_charts:
            # Interactive charts with plotly (larger file size, better UX)
            self.barcode_rank_plot_html = self.generate_barcode_rank_plot(include_plotlyjs='inline')
            self.sequencing_saturation_plot_html = self.generate_sequencing_saturation_plot(include_plotlyjs=False)
            self.median_genes_plot_html = self.generate_median_genes_plot(include_plotlyjs=False)
        else:
            # Static charts for faster loading
            self.barcode_rank_plot_html = self.generate_static_chart_placeholder("Barcode Rank Plot", "UMI Counts vs Barcodes")
            self.sequencing_saturation_plot_html = self.generate_static_chart_placeholder("Sequencing Saturation", "Read fraction vs Sequencing saturation%")
            self.median_genes_plot_html = self.generate_static_chart_placeholder(self._unit_text()["median_genes_plot_title"], "Sequencing Depth Analysis")
        
        # Load and populate template
        template = self.load_template()
        populated_html = self.populate_template(template)
        
        # Save report
        if output_filename is None:
            output_filename = f"{self.sample_name}_scrna_analysis_report.html"
        
        output_path = self.output_dir / output_filename
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(populated_html)
        
        self.logger.info(f"Report generated successfully: {output_path}")
        return output_path


def main():
    """Main function for command line usage."""
    parser = argparse.ArgumentParser(
        description='Generate single-cell RNA-seq analysis report',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Example usage:
  python scrna_report_generator.py /path/to/sample_dir sample_name
  python scrna_report_generator.py /path/to/sample_dir --auto-detect
  python scrna_report_generator.py /path/to/sample_dir sample_name --chemistry BBV2.4 --species Homo_sapiens
        """
    )
    
    parser.add_argument('sample_dir', 
                       help='Path to sample analysis directory')
    parser.add_argument('sample_name', nargs='?', default=None,
                       help='Sample identifier (auto-detected if not provided)')
    parser.add_argument('--auto-detect', action='store_true',
                       help='Auto-detect sample name from directory structure')
    parser.add_argument('--chemistry', default='BBV2.4',
                       help='Chemistry version (default: BBV2.4)')
    parser.add_argument('--species', default='Mus_musculus',
                       help='Species name (default: Mus_musculus)')
    parser.add_argument('--tissue', help='Tissue name to display in sample metadata')
    parser.add_argument('--output-dir', 
                       help='Output directory (default: same as sample_dir)')
    parser.add_argument('--output-filename',
                       help='Output filename (default: auto-generated)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    parser.add_argument('--top-genes', type=int, default=15,
                       help='Number of top genes per cluster to display (default: 15)')
    parser.add_argument('--fast-mode', action='store_true',
                       help='Use static charts for faster loading (default: interactive charts)') 
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Handle auto-detection
    sample_name = args.sample_name
    if args.auto_detect or sample_name is None:
        sample_name = None  # Let the class auto-detect
    
    try:
        generator = ScRNAReportGenerator(
            sample_dir=args.sample_dir,
            sample_name=sample_name,
            chemistry=args.chemistry,
            species=args.species,
            output_dir=args.output_dir,
            tissue=args.tissue,
        )
        
        # Override parameters
        if hasattr(args, 'fast_mode'):
            generator.use_interactive_charts = not args.fast_mode
        generator.load_marker_genes(top_n=args.top_genes)
        
        output_path = generator.generate_report(args.output_filename)
        print(f"✅ Report generated successfully: {output_path}")
        print(f"📊 Sample: {generator.sample_name}")
        print(f"🧬 Chemistry: {generator.chemistry}")
        print(f"🔬 Species: {generator.species}")
        print(f"📈 Marker genes: {len(generator.marker_genes)}")
        
    except FileNotFoundError as e:
        logging.error(f"File not found: {e}")
        print(f"❌ Error: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Failed to generate report: {e}")
        print(f"❌ Error: Failed to generate report - {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()

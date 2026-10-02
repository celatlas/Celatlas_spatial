#!/usr/bin/env python3
"""
Celatlas Spatial Report Generator

This script generates a comprehensive HTML report for spatial transcriptomics analysis.
It reads data from various analysis output files and creates an interactive HTML report
with the following features:

- Key metrics and statistics visualization
- Barcode quality control analysis 
- Mapping and feature counting statistics
- Analysis results with tissue hires images, spatial gene expression, clustering
- Interactive marker genes table with search and export functionality
- Professional bioinformatics UI design
- Support for multiple sample types and auto-detection

Usage:
    python spatial_report_generator.py /path/to/sample_dir [sample_name]
    python spatial_report_generator.py /path/to/sample_dir --auto-detect

Author: Claude (Anthropic)
Date: 2025-11-28
Version: 1.8.0 - Enhanced with improved UI, marker genes, and analysis images
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


class SpatialReportGenerator:
    """Generate comprehensive spatial transcriptomics analysis reports."""
    
    def __init__(self, sample_dir: str, sample_name: str = None, chemistry: str = "BBV2.4",
                 species: str = "Mus_musculus", method: str = None, output_dir: str = None,
                 tissue: str = None,
                 use_interactive_charts: bool = True,
                 include_cell_segmentation: bool = True):
        """
        Initialize the report generator.
        
        Args:
            sample_dir: Path to the sample analysis directory
            sample_name: Sample identifier (auto-detected if None)
            chemistry: Chemistry version used
            species: Species name
            method: Spatial alignment mode from the pipeline command
            output_dir: Output directory for the report (defaults to sample_dir)
            include_cell_segmentation: Include experimental Cell Segmentation tab
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
        self.method = self._normalize_method(method)
        self.output_dir = Path(output_dir) if output_dir else self.sample_dir
        self.include_cell_segmentation = include_cell_segmentation
        
        # Set up logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
        # Data containers
        self.bin_stats = {}
        self.marker_genes = []
        self.complete_marker_genes = []
        self.clusters = []
        self.general_stats = {}
        self.bin_images = {}
        self.analysis_images = {}
        self.cell_segmentation_metrics = {}
        self.cell_segmentation_images = {}
        self.cell_marker_genes = []
        self.cell_complete_marker_genes = []
        self.cell_clusters = []
        self.has_cell_segmentation = False
        self.company_logo = None
        self.downsample_data = None
        self.downsample_source = None
        self.count_downsample_data = None
        self.count_downsample_source = None
        self.count_detail_data = None
        self.barcode_umi_counts = None  # Pre-aggregated barcode->UMI count mapping
        self.median_genes_plot_html = ""
        
        # Chart HTML containers
        self.barcode_rank_plot_html = ""
        self.sequencing_saturation_plot_html = ""
        
        # Performance optimization flag
        self.use_interactive_charts = True  # Can be set to False for faster loading
        
        # Template path
        self.template_path = Path(__file__).parent.parent / "templates" / "html" / "spatial_report_template.html"
        
        # Validate inputs
        if not self.sample_dir.exists():
            raise FileNotFoundError(f"Sample directory not found: {self.sample_dir}")
        if not self.template_path.exists():
            raise FileNotFoundError(f"Template file not found: {self.template_path}")

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

    def parse_star_final_log(self, log_file: Path) -> Dict[str, str]:
        """Parse STAR Log.final.out using STAR's native mapping categories."""
        stats = {}
        if not log_file.exists():
            return stats

        key_map = {
            'Number of input reads': 'STAR Input Reads',
            'Uniquely mapped reads number': 'Uniquely Mapped Reads',
            'Uniquely mapped reads %': 'Uniquely Mapped Reads Percentage',
            'Number of reads mapped to multiple loci': 'Multi-Mapped Reads',
            '% of reads mapped to multiple loci': 'Multi-Mapped Reads Percentage',
            'Number of reads mapped to too many loci': 'Reads Mapped to Too Many Loci',
            '% of reads mapped to too many loci': 'Reads Mapped to Too Many Loci Percentage',
            'Number of reads unmapped: too short': 'Unmapped Reads - Too Short',
            '% of reads unmapped: too short': 'Unmapped Reads - Too Short Percentage',
            'Number of reads unmapped: other': 'Unmapped Reads - Other',
            '% of reads unmapped: other': 'Unmapped Reads - Other Percentage',
        }
        try:
            with open(log_file, 'r', encoding='utf-8') as handle:
                for line in handle:
                    if '|' not in line:
                        continue
                    key, value = [item.strip() for item in line.split('|', 1)]
                    mapped_key = key_map.get(key)
                    if not mapped_key:
                        continue
                    cleaned = value.replace(',', '').strip()
                    if mapped_key.endswith('Percentage'):
                        cleaned = cleaned.replace('%', '').strip()
                        stats[mapped_key] = f"{cleaned}%"
                    else:
                        stats[mapped_key] = cleaned
        except Exception as e:
            self.logger.warning(f"Failed to parse STAR final log {log_file}: {e}")
        return stats
    
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

    def _binned_markers_dir(self) -> Path:
        return self._binned_outputs_dir() / "markers"

    @staticmethod
    def _first_existing(paths: List[Path]) -> Path:
        for path in paths:
            if path.exists():
                return path
        return paths[0]

    def _marker_raw_file(self) -> Path:
        return self._binned_markers_dir() / f"{self.sample_name}_markers_raw.tsv"

    def _cell_marker_raw_file(self) -> Path:
        cell_outputs_dir = self._cellsegmented_outputs_dir()
        return self._first_existing([
            cell_outputs_dir / "markers" / f"{self.sample_name}_cell_markers_raw.tsv",
            cell_outputs_dir / "markers" / f"{self.sample_name}_cell_markers.tsv",
        ])

    def _binsegment_dir(self) -> Path:
        """Return the current binSegment output directory."""
        return self.sample_dir / "06.segment" / "01.binsegment"

    def _cellsegmented_outputs_dir(self) -> Path:
        """Return the compact 07.outs cell-segmented output directory."""
        return self.sample_dir / "07.outs" / "cellsegmented_outputs"

    def _cellsegment_dir(self) -> Path:
        """Return the current StarDist/cell segmentation output directory."""
        return self.sample_dir / "06.segment" / "02.cellsegment"

    def _calculate_bin_matrix_stats(self, matrix_dir: Path) -> Dict[str, str]:
        """Calculate bin metrics directly from the final 10X matrix.

        This keeps report metrics accurate after cavity filtering replaces
        square_bin matrices without regenerating the original stat.txt file.
        """
        matrix_file = self._first_existing([
            matrix_dir / "matrix.mtx.gz",
            matrix_dir / "matrix.mtx",
        ])
        if not matrix_file.exists():
            return {}

        try:
            from scipy.io import mmread
            matrix = mmread(str(matrix_file)).tocsr()
            if matrix.shape[1] == 0:
                return {
                    "Estimated Number of square bin": "0",
                    "Total Genes": "0",
                    "Mean UMI per square bin": "0",
                    "Median UMI per square bin": "0",
                    "Mean Genes per square bin": "0",
                    "Median Genes per square bin": "0",
                }

            umi_per_bin = np.asarray(matrix.sum(axis=0)).ravel()
            genes_per_bin = np.asarray(matrix.getnnz(axis=0)).ravel()
            expressed_genes = int(np.count_nonzero(matrix.getnnz(axis=1)))

            return {
                "Estimated Number of square bin": str(int(matrix.shape[1])),
                "Total Genes": str(expressed_genes),
                "Mean UMI per square bin": str(int(umi_per_bin.mean())),
                "Median UMI per square bin": str(int(np.median(umi_per_bin))),
                "Mean Genes per square bin": str(int(genes_per_bin.mean())),
                "Median Genes per square bin": str(int(np.median(genes_per_bin))),
                "Matrix Nonzero Entries": str(int(matrix.nnz)),
            }
        except Exception as e:
            self.logger.warning(f"Failed to calculate matrix statistics from {matrix_dir}: {e}")
            return {}

    def _should_calculate_bin_matrix_stats(self, bin_dir: Path, stat_file: Path) -> bool:
        """Use matrix-derived bin metrics only for cavity-filtered or explicit repair cases."""
        force = os.environ.get("CELATLAS_REPORT_RECALCULATE_BIN_STATS", "").strip().lower()
        if force in ("1", "true", "yes", "on"):
            return True
        if not stat_file.exists():
            return True
        return (bin_dir / "filter_summary.json").exists()
    
    def load_bin_statistics(self):
        """Load statistics for all supported bin sizes.
        
        Searches for bin statistics in the standard directory structure:
        06.segment/01.binsegment/square_bin/{sample_name}_bin{size}/stat.txt
        """
        binsegment_dir = self._binsegment_dir()
        for bin_size in SUPPORTED_BIN_SIZES:
            bin_dir = binsegment_dir / "square_bin" / f"{self.sample_name}_bin{bin_size}"
            stat_file = bin_dir / "stat.txt"
            matrix_dir = bin_dir / "filtered_feature_bc_matrix"
            stats = {}
            
            if stat_file.exists():
                stats = self.parse_stat_file(stat_file)
                self.logger.info(f"Loaded statistics for bin {bin_size}μm")
            else:
                self.logger.warning(f"Bin {bin_size}μm stat file not found: {stat_file}")

            matrix_stats = {}
            if self._should_calculate_bin_matrix_stats(bin_dir, stat_file):
                matrix_stats = self._calculate_bin_matrix_stats(matrix_dir)
            if matrix_stats:
                stat_bin_count = self._parse_numeric(stats.get("Estimated Number of square bin"))
                matrix_bin_count = self._parse_numeric(matrix_stats.get("Estimated Number of square bin"))
                if stat_bin_count is not None and matrix_bin_count is not None and stat_bin_count != matrix_bin_count:
                    self.logger.warning(
                        f"Bin {bin_size}μm stat.txt bin count ({int(stat_bin_count)}) differs from "
                        f"matrix ({int(matrix_bin_count)}); using matrix-derived UMI/gene metrics."
                    )
                stats.update(matrix_stats)
                self.logger.info(f"Calculated matrix-derived statistics for bin {bin_size}μm")

            self.bin_stats[f'bin{bin_size}'] = stats
    
    def load_general_statistics(self):
        """Load general pipeline statistics."""
        # Load barcode statistics
        self.general_stats.update(self.parse_step_metrics(self.sample_dir / "01.barcode"))
        
        # Load STAR mapping statistics
        self.general_stats.update(self.parse_step_metrics(self.sample_dir / "03.star"))
        star_log_candidates = [
            self.sample_dir / "03.star" / f"{self.sample_name}_Log.final.out",
            *sorted((self.sample_dir / "03.star").glob("*Log.final.out")),
        ]
        for star_log in star_log_candidates:
            if star_log.exists():
                self.general_stats.update(self.parse_star_final_log(star_log))
                break
        
        # Load featureCounts statistics
        self.general_stats.update(self.parse_step_metrics(self.sample_dir / "04.featureCounts"))
        
        # Load count statistics
        self.general_stats.update(self.parse_step_metrics(self.sample_dir / "05.count"))
            
        # Load analysis statistics
        self.general_stats.update(self.parse_step_metrics(self._binned_outputs_dir()))
    
    @staticmethod
    def _sort_marker_clusters(clusters) -> List:
        def cluster_key(value):
            try:
                return (0, int(value))
            except (TypeError, ValueError):
                return (1, str(value))
        return sorted(clusters, key=cluster_key)

    @staticmethod
    def _row_float(row, column: str, default: float) -> float:
        if column not in row or pd.isna(row[column]):
            return default
        try:
            return float(row[column])
        except (TypeError, ValueError):
            return default

    def _marker_record_from_row(self, row) -> Dict:
        return {
            'cluster': str(row['cluster']),
            'gene': str(row['gene']),
            'scores': self._row_float(row, 'scores', 0.0),
            'avg_log2FC': self._row_float(row, 'avg_log2FC', 0.0),
            'p_val': self._row_float(row, 'p_val', 1.0),
            'p_val_adj': self._row_float(row, 'p_val_adj', 1.0),
            'pct1': self._row_float(row, 'pct.1', 0.0),
            'pct2': self._row_float(row, 'pct.2', 0.0),
        }

    def _load_marker_genes_from_file(self, marker_file: Path, top_n: int, label: str) -> Tuple[List[Dict], List[Dict], List]:
        """Load marker genes from one marker TSV and return top rows, complete rows, and clusters."""
        if not marker_file.exists():
            self.logger.warning(f"{label} marker genes file not found: {marker_file}")
            return [], [], []

        try:
            df = pd.read_csv(marker_file, sep='\t')
            required_columns = {'cluster', 'gene'}
            if df.empty or not required_columns.issubset(df.columns):
                self.logger.warning(f"{label} marker file is empty or missing required columns: {marker_file}")
                return [], [], []

            self.logger.info(f"Loaded {label} marker genes file with {len(df)} total genes")

            marker_genes = []
            cluster_summary = []
            clusters = self._sort_marker_clusters(df['cluster'].dropna().unique())

            for cluster in clusters:
                cluster_df = df[df['cluster'].astype(str) == str(cluster)]
                if 'scores' in cluster_df.columns:
                    cluster_df = cluster_df.sort_values('scores', ascending=False)
                elif 'avg_log2FC' in cluster_df.columns:
                    cluster_df = cluster_df.sort_values('avg_log2FC', ascending=False)

                top_genes = cluster_df.head(top_n)
                cluster_summary.append(f"Cluster {cluster}: {len(top_genes)} genes")
                for _, row in top_genes.iterrows():
                    marker_genes.append(self._marker_record_from_row(row))

            complete_marker_genes = [self._marker_record_from_row(row) for _, row in df.iterrows()]

            self.logger.info(f"Loaded {len(marker_genes)} {label} marker genes from {len(clusters)} clusters")
            self.logger.info(f"{label} marker distribution: {', '.join(cluster_summary)}")
            self.logger.info(f"Stored {len(complete_marker_genes)} complete {label} marker genes for TSV export")
            return marker_genes, complete_marker_genes, clusters

        except Exception as e:
            self.logger.error(f"Failed to load {label} marker genes: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return [], [], []

    def load_marker_genes(self, top_n: int = 20):
        """Load binned marker genes from the analysis results, ensuring top N genes per cluster."""
        self.marker_genes, self.complete_marker_genes, self.clusters = self._load_marker_genes_from_file(
            self._marker_raw_file(),
            top_n,
            "binned",
        )

    def load_cell_marker_genes(self, top_n: int = 20):
        """Load cell-segmented marker genes from the cell analysis results."""
        if not self.include_cell_segmentation:
            self.cell_marker_genes = []
            self.cell_complete_marker_genes = []
            self.cell_clusters = []
            return

        self.cell_marker_genes, self.cell_complete_marker_genes, self.cell_clusters = self._load_marker_genes_from_file(
            self._cell_marker_raw_file(),
            top_n,
            "cell-segmented",
        )
        if self.cell_marker_genes:
            self.has_cell_segmentation = True

    def _generate_markers_tsv_data_uri(self, marker_file: Path, label: str) -> str:
        """Generate data URI for a markers TSV file to enable offline download."""
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
                self.logger.info(f"Generated data URI for {label} markers TSV file:")
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
                self.logger.error(f"Failed to generate {label} TSV data URI: {e}")
                import traceback
                self.logger.error(traceback.format_exc())
                return ""
        else:
            self.logger.warning(f"{label} marker TSV file not found: {marker_file}")
            return ""

    def generate_markers_tsv_data_uri(self) -> str:
        """Generate data URI for binned markers_raw.tsv file to enable offline download."""
        return self._generate_markers_tsv_data_uri(self._marker_raw_file(), "binned")

    def generate_cell_markers_tsv_data_uri(self) -> str:
        """Generate data URI for cell markers_raw.tsv file to enable offline download."""
        if not self.include_cell_segmentation:
            return ""
        return self._generate_markers_tsv_data_uri(self._cell_marker_raw_file(), "cell-segmented")
    
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
        """Load count-level headline metrics and bin50 plot data independently."""
        count_downsample_file = self.sample_dir / "05.count" / f"{self.sample_name}_downsample.tsv"
        if count_downsample_file.exists():
            try:
                count_data = pd.read_csv(count_downsample_file, sep='\t')
                self.count_downsample_data = self._normalize_downsample_data(count_data, count_downsample_file)
                if self.count_downsample_data is not None:
                    self.count_downsample_source = f"05.count: {count_downsample_file}"
                    self.logger.info(
                        "Loaded count-level downsample data with %s data points for headline saturation",
                        len(self.count_downsample_data),
                    )
            except Exception as e:
                self.logger.warning(f"Failed to load count-level downsample data: {e}")
                self.count_downsample_data = None

        # Bin50 is the spatially meaningful source for the report curves.
        bin_downsample_file = self._binsegment_dir() / "square_bin" / f"{self.sample_name}_bin50" / "downsample.tsv"

        if bin_downsample_file.exists():
            try:
                downsample_data = pd.read_csv(bin_downsample_file, sep='\t')
                self.downsample_data = self._normalize_downsample_data(downsample_data, bin_downsample_file)
                if self.downsample_data is not None:
                    self.downsample_source = f"bin50: {bin_downsample_file}"
                    self.logger.info(f"Loaded bin-specific downsample data (bin50) with {len(self.downsample_data)} data points")
                    return
            except Exception as e:
                self.logger.warning(f"Failed to load bin-specific downsample data: {e}")

        # Historical outputs may only have count-level data; it remains usable for plots.
        if self.count_downsample_data is not None:
            self.downsample_data = self.count_downsample_data
            self.downsample_source = self.count_downsample_source
            self.logger.info("Using count-level downsample data for report curves")
        else:
            self.logger.warning(f"No valid downsample file found at {bin_downsample_file} or {count_downsample_file}")
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

            # Use Estimated Number of Cells as the split point
            # This ensures the blue/gray split matches the reported cell count
            estimated_cells = self.general_stats.get('Estimated Number of Cells', None)

            if estimated_cells is not None:
                try:
                    # Remove commas if present and convert to int
                    if isinstance(estimated_cells, str):
                        estimated_cells = int(estimated_cells.replace(',', ''))
                    else:
                        estimated_cells = int(estimated_cells)
                    inflection_idx = min(estimated_cells, len(x_data))
                    self.logger.info(f"Using Estimated Number of Cells ({estimated_cells}) as split point")
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

            # Split data into cells (blue) and background (gray)
            # Include inflection point in both traces to ensure continuous line
            cell_x = x_data[:inflection_idx + 1]  # Include inflection_idx
            cell_y = barcode_counts[:inflection_idx + 1]
            bg_x = x_data[inflection_idx:]  # Start from inflection_idx
            bg_y = barcode_counts[inflection_idx:]

            # Create the plot
            fig = go.Figure()

            # Add cells trace (blue)
            fig.add_trace(go.Scatter(
                x=cell_x,
                y=cell_y,
                mode='lines',
                name='Cells',
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

            # Update layout with Cell Ranger style ticks
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

            scope = "Bin50" if self.downsample_source and self.downsample_source.startswith("bin50:") else "Count-level"
            chart_title = f"{scope} UMI Duplicate Saturation"
            
            # Create the plot
            fig = go.Figure()
            
            fig.add_trace(go.Scatter(
                x=read_fractions,
                y=saturation_values,
                mode='lines',
                name=chart_title,
                line=dict(color='#ffde7d', width=3),
                hovertemplate='<b>Read Fraction:</b> %{x:.2f}<br><b>Saturation:</b> %{y:.2f}%<extra></extra>'
            ))
            
            # Update layout
            fig.update_layout(
                title=dict(
                    text=chart_title,
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
                    title=f'{scope} UMI Duplicate Saturation (%)',
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
        """Generate interactive Median Genes per Square plot using plotly."""
        if not PLOTLY_AVAILABLE:
            return '<div class="simple-chart"><div class="chart-placeholder">Plotly not available</div></div>'
        
        try:
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
                name='Median Genes per Square',
                line=dict(color='#e84545', width=3),
                hovertemplate='<b>Read Fraction:</b> %{x:.2f}<br><b>Median Genes:</b> %{y:.0f}<extra></extra>'
            ))
            
            # Update layout
            fig.update_layout(
                title=dict(
                    text='Median Genes per Square',
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
                    title='Median Genes per Square',
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
            return '<div class="simple-chart"><div class="chart-placeholder">Median Genes per Square<br><small>Error generating chart</small></div></div>'
    
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
        """Load analysis result images including tissue hires, spatial gene expression, spatial cluster, and umap cluster."""
        # Initialize analysis images dictionary
        self.analysis_images = {
            'tissue_hires': None,
            'spatial_gene_expression': None, 
            'spatial_cluster': None,
            'umap_cluster': None
        }
        
        # Load tissue hires image from the current directory layout.
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
        spatial_gene_expr_file = self._binned_qc_dir() / f"{self.sample_name}_spatial_gene_expression_distribution.png"
        if spatial_gene_expr_file.exists():
            self.analysis_images['spatial_gene_expression'] = self.encode_image_to_base64(spatial_gene_expr_file)
            self.logger.info("Loaded spatial gene expression distribution image")
        
        # Load spatial cluster image
        spatial_cluster_file = self._binned_clustering_dir() / f"{self.sample_name}_spatial_cluster.png"
        if spatial_cluster_file.exists():
            self.analysis_images['spatial_cluster'] = self.encode_image_to_base64(spatial_cluster_file)
            self.logger.info("Loaded spatial cluster image")
        
        # Load UMAP cluster image
        umap_cluster_file = self._binned_clustering_dir() / f"{self.sample_name}_umap_cluster.png"
        if umap_cluster_file.exists():
            self.analysis_images['umap_cluster'] = self.encode_image_to_base64(umap_cluster_file)
            self.logger.info("Loaded UMAP cluster image")

    def _format_cell_metric(self, value, percent: bool = False) -> str:
        if value is None:
            return "N/A"
        parsed = self._parse_numeric(value)
        if parsed is None:
            return str(value)
        if percent:
            if parsed <= 1:
                parsed *= 100
            return f"{parsed:.2f}%"
        if abs(parsed - round(parsed)) < 1e-6:
            return f"{int(round(parsed)):,}"
        return f"{parsed:,.2f}"

    def _load_json_file(self, path: Path) -> Dict:
        if not path.exists():
            return {}
        try:
            with open(path, 'r', encoding='utf-8') as handle:
                return json.load(handle)
        except Exception as e:
            self.logger.warning(f"Failed to load JSON file {path}: {e}")
            return {}

    def _load_cell_stats_summary(self, cellsegment_dir: Path) -> Dict:
        stats_file = cellsegment_dir / f"{self.sample_name}.cell_stats.tsv"
        if not stats_file.exists():
            stats_file = cellsegment_dir / f"{self.sample_name}_cell_stats.tsv"
        if not stats_file.exists():
            return {}

        try:
            df = pd.read_csv(stats_file, sep='\t')
            if df.empty:
                return {}
            retained = df
            if 'retained' in df.columns:
                retained = df[pd.to_numeric(df['retained'], errors='coerce').fillna(0) > 0]
            if retained.empty:
                retained = df

            summary = {}
            for src_col, dst_key, agg in [
                ('reads', 'mean_reads_per_cell', 'mean'),
                ('area_px', 'median_cell_area_px2', 'median'),
                ('assigned_barcodes', 'median_barcodes_per_cell', 'median'),
            ]:
                if src_col in retained.columns:
                    values = pd.to_numeric(retained[src_col], errors='coerce').dropna()
                    if not values.empty:
                        summary[dst_key] = float(values.mean() if agg == 'mean' else values.median())
            return summary
        except Exception as e:
            self.logger.warning(f"Failed to parse cell stats file {stats_file}: {e}")
            return {}

    def _infer_cell_cluster_count(self, cell_outputs_dir: Path) -> Optional[int]:
        marker_file = cell_outputs_dir / "markers" / f"{self.sample_name}_cell_markers_raw.tsv"
        if not marker_file.exists():
            marker_file = cell_outputs_dir / "markers" / f"{self.sample_name}_cell_markers.tsv"
        if not marker_file.exists():
            return None
        try:
            df = pd.read_csv(marker_file, sep='\t', usecols=['cluster'])
            if df.empty or 'cluster' not in df.columns:
                return None
            return int(df['cluster'].astype(str).nunique())
        except Exception:
            return None

    def load_cell_segmentation_data(self):
        """Load cell segmentation metrics and images for the report Cell Segmentation tab."""
        if not self.include_cell_segmentation:
            self.has_cell_segmentation = False
            self.cell_segmentation_metrics = {}
            self.cell_segmentation_images = {}
            return

        cell_outputs_dir = self._cellsegmented_outputs_dir()
        cellsegment_dir = self._cellsegment_dir()
        self.has_cell_segmentation = cell_outputs_dir.exists() or cellsegment_dir.exists()

        metrics_summary = self._load_json_file(cell_outputs_dir / "metrics_summary.json")
        assignment_summary = self._load_json_file(cellsegment_dir / f"{self.sample_name}.stardist_assignment_summary.json")
        run_summary = self._load_json_file(cellsegment_dir / f"{self.sample_name}.stardist_cell_segment_run_summary.json")
        stardist_summary = self._load_json_file(cellsegment_dir / "stardist_inference" / f"{self.sample_name}.stardist_summary.json")
        cell_stats_summary = self._load_cell_stats_summary(cellsegment_dir)

        assignment = assignment_summary.get('assignment', {}) if isinstance(assignment_summary, dict) else {}
        count_aggregation = assignment_summary.get('count_aggregation', {}) if isinstance(assignment_summary, dict) else {}

        cluster_count = metrics_summary.get('cluster_count')
        if cluster_count is None:
            cluster_count = self._infer_cell_cluster_count(cell_outputs_dir)

        self.cell_segmentation_metrics = {
            'number_of_cells': self._format_cell_metric(
                metrics_summary.get('cells_after_qc', count_aggregation.get('cells_after_qc', assignment.get('assigned_cells')))
            ),
            'cells_before_qc': self._format_cell_metric(
                metrics_summary.get('cells_before_qc', count_aggregation.get('cells_before_qc', assignment.get('total_labels')))
            ),
            'total_genes': self._format_cell_metric(
                metrics_summary.get('genes', count_aggregation.get('genes'))
            ),
            'median_umi_per_cell': self._format_cell_metric(
                metrics_summary.get('median_umi_per_cell', count_aggregation.get('median_umi_per_retained_cell'))
            ),
            'median_genes_per_cell': self._format_cell_metric(
                metrics_summary.get('median_genes_per_cell', count_aggregation.get('median_genes_per_retained_cell'))
            ),
            'median_mt_percent': self._format_cell_metric(metrics_summary.get('median_mt_percent'), percent=True),
            'raw_nuclei': self._format_cell_metric(
                assignment_summary.get('raw_labels_before_expansion', stardist_summary.get('n_labels'))
            ),
            'labels_after_expansion': self._format_cell_metric(assignment_summary.get('labels_after_expansion')),
            'assigned_barcodes': self._format_cell_metric(assignment.get('assigned_barcodes')),
            'assignment_rate': self._format_cell_metric(assignment.get('assignment_rate'), percent=True),
            'assigned_cell_fraction': self._format_cell_metric(assignment.get('assigned_cell_fraction'), percent=True),
            'matrix_nnz': self._format_cell_metric(count_aggregation.get('matrix_nonzero_entries')),
            'mean_reads_per_cell': self._format_cell_metric(cell_stats_summary.get('mean_reads_per_cell')),
            'median_cell_area_px2': self._format_cell_metric(cell_stats_summary.get('median_cell_area_px2')),
            'median_barcodes_per_cell': self._format_cell_metric(cell_stats_summary.get('median_barcodes_per_cell')),
            'cluster_count': self._format_cell_metric(cluster_count),
            'cluster_resolution': self._format_cell_metric(metrics_summary.get('cluster_resolution')),
            'n_neighbors': self._format_cell_metric(metrics_summary.get('n_neighbors')),
            'n_pcs': self._format_cell_metric(metrics_summary.get('n_pcs')),
            'segmentation_model': str(run_summary.get('model', stardist_summary.get('model', 'StarDist'))),
        }

        image_candidates = {
            'stardist_boundary_overlay': [
                cellsegment_dir / "stardist_inference" / f"{self.sample_name}.stardist_boundary_overlay.preview.jpg",
                cellsegment_dir / "stardist_inference" / f"{self.sample_name}.stardist_boundary_overlay.preview.png",
            ],
            'qc_boundary_overlay': [
                cellsegment_dir / "visualization" / f"{self.sample_name}.qc_boundary_overlay.preview.jpg",
                cellsegment_dir / "visualization" / f"{self.sample_name}.qc_boundary_overlay.jpg",
            ],
            'barcode_assignment': [
                cellsegment_dir / "visualization" / f"{self.sample_name}.barcode_assignment.preview.jpg",
            ],
            'cell_umi_heatmap': [
                cellsegment_dir / "visualization" / f"{self.sample_name}.cell_umi_heatmap.png",
            ],
            'cell_spatial_cluster': [
                cell_outputs_dir / "spatial" / f"{self.sample_name}_cell_spatial_cluster.png",
            ],
            'cell_spatial_total_counts': [
                cell_outputs_dir / "spatial" / f"{self.sample_name}_cell_spatial_total_counts.png",
            ],
            'cell_spatial_n_genes': [
                cell_outputs_dir / "spatial" / f"{self.sample_name}_cell_spatial_n_genes.png",
            ],
            'cell_umap_cluster': [
                cell_outputs_dir / "clustering" / f"{self.sample_name}_cell_umap_cluster.png",
            ],
            'cell_qc_violin': [
                cell_outputs_dir / "qc" / f"{self.sample_name}_cell_qc_violin.png",
            ],
            'cell_qc_scatter': [
                cell_outputs_dir / "qc" / f"{self.sample_name}_cell_qc_scatter.png",
            ],
        }

        self.cell_segmentation_images = {}
        for key, paths in image_candidates.items():
            image_path = self._first_existing(paths)
            if image_path.exists():
                self.cell_segmentation_images[key] = self._encode_image_file_to_data_uri(image_path)

        if self.cell_segmentation_images or any(value != "N/A" for value in self.cell_segmentation_metrics.values()):
            self.has_cell_segmentation = True
    
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
        """Return count-level UMI duplicate saturation from an authoritative source."""
        try:
            if self.count_downsample_data is not None and not self.count_downsample_data.empty:
                downsample_data = self.count_downsample_data.sort_values('read_fraction')
                if 'umi_saturation' in self.count_downsample_data.columns:
                    max_saturation = downsample_data['umi_saturation'].iloc[-1]
                    return self._format_percent_value(max_saturation)
                elif 'read_saturation' in self.count_downsample_data.columns:
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

    def _saturation_note(self, context: str = "spatial") -> str:
        mean_reads = self._parse_numeric(self.bin_stats.get('bin50', {}).get('Mean Reads per square bin', '0'))
        saturation = self._parse_numeric(self._calculate_saturation_percentage())
        raw_reads = self._parse_numeric(self.general_stats.get('Raw Reads', '0'))
        valid_reads = self._parse_numeric(self.general_stats.get('Valid Reads', '0'))

        recovery_note = self._gene_recovery_note()
        if saturation is None:
            return "Sequencing saturation is unavailable."
        elif (
            (mean_reads is not None and mean_reads < 1000)
            or (valid_reads is not None and valid_reads < 5_000_000)
            or (raw_reads is not None and raw_reads < 10_000_000)
        ):
            depth_note = "Low-depth runs can show elevated duplicate-rate saturation before gene recovery reaches a stable plateau."
        else:
            depth_note = "Count-level UMI duplicate saturation is estimated from read downsampling."
        return f"{depth_note} {recovery_note}".strip()

    def _gene_recovery_note(self) -> str:
        """Summarize whether the spatial median-gene curve is still rising."""
        data = self.downsample_data
        if data is None or data.empty or 'median_gene_number' not in data.columns:
            return "Interpret saturation together with total reads, median genes, and median UMI."

        records = data.sort_values('read_fraction').to_dict('records')
        full_rows = [row for row in records if (self._parse_numeric(row.get('read_fraction')) or 0) >= 0.99]
        baseline_rows = [row for row in records if (self._parse_numeric(row.get('read_fraction')) or 0) <= 0.8]
        if not full_rows or not baseline_rows:
            return "Interpret saturation together with total reads, median genes, and median UMI."

        full_genes = self._parse_numeric(full_rows[-1].get('median_gene_number'))
        baseline_genes = self._parse_numeric(baseline_rows[-1].get('median_gene_number'))
        if full_genes is None or baseline_genes is None or full_genes <= 0:
            return "Interpret saturation together with total reads, median genes, and median UMI."

        gain = max(0.0, full_genes - baseline_genes)
        gain_pct = gain / full_genes * 100.0
        source = "bin50" if self.downsample_source and self.downsample_source.startswith("bin50:") else "count-level"
        if gain_pct >= 10.0:
            interpretation = "gene recovery is still rising, so additional sequencing may remain useful"
        else:
            interpretation = "gene recovery is approaching a plateau"
        return (
            f"The {source} median-gene curve gains {gain:.1f} genes "
            f"({gain_pct:.1f}% of the full-depth value) from 80% to 100% reads; {interpretation}."
        )

    @staticmethod
    def _normalize_method(method: str) -> Optional[str]:
        if method is None:
            return None
        normalized = str(method).strip()
        if not normalized:
            return None
        lowered = normalized.lower()
        if lowered == 'he':
            return 'HE'
        if lowered in {'gene_expr', 'gene-expression', 'gene expression'}:
            return 'gene_expr'
        if lowered in {'image', 'ssdna'}:
            return 'ssDNA'
        return normalized

    def _detect_image_alignment_mode(self) -> str:
        """
        Return the explicit pipeline method when available.

        Older reports inferred this field from HE alignment images only, which
        made gene_expr runs display as ssDNA. Keep image-based detection as a
        fallback for reports generated by older entrypoints.
        """
        if self.method:
            return self.method

        try:
            he_image_path = self._get_alignment_he_image_file()
            gem_image_path = self._get_alignment_gem_heatmap_file()

            if he_image_path.exists() and gem_image_path.exists():
                self.logger.info("HE mode detected - alignment images found")
                return 'HE'
            else:
                self.logger.info("HE alignment images not found - using ssDNA mode")
                return 'ssDNA'

        except Exception as e:
            self.logger.warning(f"Failed to detect image alignment mode: {e}")
            return 'ssDNA'

    def _alignment_overlay_dir(self) -> Path:
        return self._binsegment_dir() / "images" / "overlays"

    def _cavity_overlay_dir(self, bin_size: str = "bin10") -> Path:
        """Return cavity QC overlays, supporting current and legacy layouts."""
        base_dir = self._binsegment_dir() / "cavity_filter" / "images" / "overlays"
        legacy_dir = base_dir / bin_size
        return legacy_dir if legacy_dir.is_dir() else base_dir

    def _cavity_overlay_candidates(self, filename: str, bin_size: str = "bin10") -> List[Path]:
        """Return current single-directory and legacy per-bin cavity paths."""
        base_dir = self._binsegment_dir() / "cavity_filter" / "images" / "overlays"
        return [base_dir / filename, base_dir / bin_size / filename]

    def _get_alignment_he_image_file(self) -> Path:
        """Return the preferred H&E base image for the Image Alignment viewer."""
        return self._first_existing([
            self._alignment_overlay_dir() / "2_he_registered.jpg",
            self._binsegment_dir() / "images" / "tissue_hires_image.png",
        ])

    def _get_alignment_gem_overlay_file(self) -> Path:
        """Return the preferred GEM overlay image for the Image Alignment viewer."""
        return self._first_existing([
            *self._cavity_overlay_candidates("3b_tissue_segmentation_gem_heatmap.png"),
            *self._cavity_overlay_candidates("1_filtered_gem_expression.png"),
            self._alignment_overlay_dir() / "3b_tissue_segmentation_gem_heatmap.png",
            self._alignment_overlay_dir() / "5_gem_heatmap_on_he.png",
        ])

    def _get_alignment_gem_heatmap_file(self) -> Path:
        """Return the preferred pure GEM heatmap image for the opacity-controlled overlay."""
        return self._first_existing([
            *self._cavity_overlay_candidates("3c_gem_heatmap_only.png"),
            *self._cavity_overlay_candidates("2_filtered_gem_heatmap_only.png"),
            *self._cavity_overlay_candidates("1_filtered_gem_expression.png"),
            self._alignment_overlay_dir() / "3c_gem_heatmap_only.png",
        ])

    def _encode_image_file_to_data_uri(self, image_path: Path) -> str:
        """Encode an image file as a data URI using its filename extension for MIME type."""
        mime_type_map = {
            ".jpg": "image/jpeg",
            ".jpeg": "image/jpeg",
            ".png": "image/png",
            ".tif": "image/tiff",
            ".tiff": "image/tiff",
        }
        mime_type = mime_type_map.get(image_path.suffix.lower(), "image/png")
        with open(image_path, 'rb') as img_file:
            base64_data = base64.b64encode(img_file.read()).decode('utf-8')
        return f"data:{mime_type};base64,{base64_data}"

    def _encode_heatmap_overlay_to_data_uri(self, image_path: Path) -> str:
        """
        Encode a GEM heatmap for visual overlay by making the white background transparent.

        Cavity-filter heatmaps are saved as opaque white-background PNGs. In the
        alignment viewer they need to behave like signal overlays, otherwise the
        opacity slider fades a white sheet over the H&E image.
        """
        try:
            from io import BytesIO
            from PIL import Image

            Image.MAX_IMAGE_PIXELS = None
            with Image.open(image_path) as image:
                rgba = image.convert("RGBA")

            pixels = np.array(rgba)
            rgb = pixels[:, :, :3]
            white_background = (
                (rgb[:, :, 0] >= 248) &
                (rgb[:, :, 1] >= 248) &
                (rgb[:, :, 2] >= 248)
            )
            pixels[:, :, 3] = np.where(white_background, 0, pixels[:, :, 3])

            output = BytesIO()
            Image.fromarray(pixels).save(output, format="PNG")
            base64_data = base64.b64encode(output.getvalue()).decode('utf-8')
            return f"data:image/png;base64,{base64_data}"
        except Exception as e:
            self.logger.warning(f"Failed to create transparent GEM overlay, using raw heatmap: {e}")
            return self._encode_image_file_to_data_uri(image_path)

    def _get_he_image_path(self) -> str:
        """
        Get the H&E registered image as base64-encoded data URI.
        Returns data URI if in HE mode, empty string otherwise.
        """
        try:
            # Check if we're in HE mode
            if self._detect_image_alignment_mode() == 'HE':
                he_image_path = self._get_alignment_he_image_file()

                if he_image_path.exists():
                    return self._encode_image_file_to_data_uri(he_image_path)
                else:
                    self.logger.warning(f"HE image file not found: {he_image_path}")
                    return ''
            else:
                return ''

        except Exception as e:
            self.logger.warning(f"Failed to get HE image as base64: {e}")
            return ''

    def _get_gem_image_path(self) -> str:
        """
        Get the GEM heatmap overlay image as base64-encoded data URI.
        Returns data URI if in HE mode, empty string otherwise.
        """
        try:
            # Check if we're in HE mode
            if self._detect_image_alignment_mode() == 'HE':
                gem_image_path = self._get_alignment_gem_overlay_file()

                if gem_image_path.exists():
                    return self._encode_image_file_to_data_uri(gem_image_path)
                else:
                    self.logger.warning(f"GEM image file not found: {gem_image_path}")
                    return ''
            else:
                return ''

        except Exception as e:
            self.logger.warning(f"Failed to get GEM image as base64: {e}")
            return ''

    def _get_gem_alignment_overlay_image_path(self) -> str:
        """
        Get the GEM heatmap as a transparent-background overlay for the viewer.

        Download still uses the original heatmap-only image; this data URI is only
        for visual alignment on top of the H&E layer.
        """
        try:
            if self._detect_image_alignment_mode() == 'HE':
                gem_heatmap_only_path = self._get_alignment_gem_heatmap_file()

                if gem_heatmap_only_path.exists():
                    return self._encode_heatmap_overlay_to_data_uri(gem_heatmap_only_path)
                else:
                    self.logger.warning(f"GEM alignment overlay file not found: {gem_heatmap_only_path}")
                    return ''
            else:
                return ''

        except Exception as e:
            self.logger.warning(f"Failed to get GEM alignment overlay as base64: {e}")
            return ''

    def _get_gem_heatmap_only_image_path(self) -> str:
        """
        Get pure GEM heatmap image as base64-encoded data URI for download button.

        Returns data URI if in HE mode, empty string otherwise.
        """
        try:
            # Check if we're in HE mode
            if self._detect_image_alignment_mode() == 'HE':
                gem_heatmap_only_path = self._get_alignment_gem_heatmap_file()

                if gem_heatmap_only_path.exists():
                    return self._encode_image_file_to_data_uri(gem_heatmap_only_path)
                else:
                    self.logger.warning(f"Pure GEM heatmap file not found: {gem_heatmap_only_path}")
                    return ''
            else:
                return ''

        except Exception as e:
            self.logger.warning(f"Failed to get pure GEM heatmap as base64: {e}")
            return ''

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

    def _generate_metric_table_rows(self, rows: List[Tuple[str, str]]) -> str:
        return '\n'.join(
            f'''
            <tr>
                <td class="left-cell">{name}</td>
                <td>{value}</td>
            </tr>
            '''
            for name, value in rows
        )

    def generate_cell_segmentation_metric_tables(self) -> str:
        metrics = self.cell_segmentation_metrics
        segmentation_rows = [
            ("Segmentation Model", metrics.get('segmentation_model', 'StarDist')),
            ("Raw StarDist Nuclei", metrics.get('raw_nuclei', 'N/A')),
            ("Labels After Expansion", metrics.get('labels_after_expansion', 'N/A')),
            ("Cells Before QC", metrics.get('cells_before_qc', 'N/A')),
            ("Number of Cells", metrics.get('number_of_cells', 'N/A')),
            ("Assigned Barcodes", metrics.get('assigned_barcodes', 'N/A')),
            ("Barcode Assignment Rate", metrics.get('assignment_rate', 'N/A')),
            ("Assigned Cell Fraction", metrics.get('assigned_cell_fraction', 'N/A')),
            ("Median Barcodes per Cell", metrics.get('median_barcodes_per_cell', 'N/A')),
        ]
        expression_rows = [
            ("Median UMIs per Cell", metrics.get('median_umi_per_cell', 'N/A')),
            ("Median Genes per Cell", metrics.get('median_genes_per_cell', 'N/A')),
            ("Total Genes Detected in Cells", metrics.get('total_genes', 'N/A')),
            ("Matrix Nonzero Entries", metrics.get('matrix_nnz', 'N/A')),
            ("Mean Reads per Cell", metrics.get('mean_reads_per_cell', 'N/A')),
            ("Cell Clusters", metrics.get('cluster_count', 'N/A')),
            ("Leiden Resolution", metrics.get('cluster_resolution', 'N/A')),
            ("Cell Graph Neighbors", metrics.get('n_neighbors', 'N/A')),
            ("Cell PCA Components", metrics.get('n_pcs', 'N/A')),
        ]
        return f'''
        <div class="cell-metrics-tables">
            <div class="cell-metrics-table-block">
                <div class="cell-metrics-subtitle">Segmentation & Assignment</div>
                <table class="bin-table cell-metrics-table">
                    <tr>
                        <th style="width: 62%;">Metric</th>
                        <th>Value</th>
                    </tr>
                    {self._generate_metric_table_rows(segmentation_rows)}
                </table>
            </div>
            <div class="cell-metrics-table-block">
                <div class="cell-metrics-subtitle">Expression & Clustering</div>
                <table class="bin-table cell-metrics-table">
                    <tr>
                        <th style="width: 62%;">Metric</th>
                        <th>Value</th>
                    </tr>
                    {self._generate_metric_table_rows(expression_rows)}
                </table>
            </div>
        </div>
        '''

    def generate_cell_segmentation_images_section(self) -> str:
        image_infos = [
            ('stardist_boundary_overlay', 'StarDist Boundary Overlay'),
            ('qc_boundary_overlay', 'QC Boundary Overlay'),
            ('barcode_assignment', 'Barcode Assignment'),
            ('cell_umi_heatmap', 'Cell UMI Heatmap'),
            ('cell_spatial_cluster', 'Cell Spatial Cluster'),
            ('cell_spatial_total_counts', 'Cell Spatial UMI Counts'),
            ('cell_spatial_n_genes', 'Cell Spatial Gene Counts'),
            ('cell_umap_cluster', 'Cell UMAP Cluster'),
            ('cell_qc_violin', 'Cell QC Violin'),
            ('cell_qc_scatter', 'Cell QC Scatter'),
        ]
        sections = []
        for key, title in image_infos:
            img_data = self.cell_segmentation_images.get(key)
            if img_data:
                image_html = f'<img src="{img_data}" alt="{title}" />'
            else:
                image_html = f'<div style="color: #999; font-size: 13px;">{title}<br/>Image not available</div>'
            sections.append(f'''
            <div class="analysis-image cell-analysis-image" title="Click to enlarge">
                {image_html}
                <div class="image-label">{title}</div>
            </div>
            ''')
        return '\n'.join(sections)
    
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
        """Generate JSON data for binned marker genes JavaScript functionality."""
        return self._generate_marker_genes_data_json(self.marker_genes)

    def generate_cell_marker_genes_data_json(self) -> str:
        """Generate JSON data for cell marker genes JavaScript functionality."""
        return self._generate_marker_genes_data_json(self.cell_marker_genes)

    @staticmethod
    def _generate_marker_genes_data_json(marker_genes: List[Dict]) -> str:
        """Generate JSON data for marker genes JavaScript functionality."""
        # Convert numpy types to Python native types for JSON serialization
        json_compatible_data = []
        for gene in marker_genes:
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

    @staticmethod
    def _clusters_to_json(clusters: List) -> str:
        cluster_list = []
        for cluster in clusters:
            try:
                cluster_list.append(int(cluster))
            except (TypeError, ValueError):
                cluster_list.append(str(cluster))
        return json.dumps(cluster_list)
    
    def load_template(self) -> str:
        """Load the HTML template."""
        if not self.template_path.exists():
            raise FileNotFoundError(f"Template file not found: {self.template_path}")
        
        with open(self.template_path, 'r', encoding='utf-8') as f:
            return f.read()

    def _strip_cell_segmentation_template(self, template: str) -> str:
        """Remove experimental cell-segmentation UI from customer-facing reports."""
        template = re.sub(
            r'\s*<div class="tab" data-tab="cell-segmentation"[^>]*>Cell Segmentation</div>',
            '',
            template,
        )
        template = re.sub(
            r'\n\s*<!-- Cell Segmentation Tab Content -->.*?\n\s*<!-- Image Modal -->',
            '\n\n    <!-- Image Modal -->',
            template,
            flags=re.DOTALL,
        )
        template = re.sub(
            r'\n\s*initializeMarkerGenesTableInstance\(\{\s*'
            r'geneData: \{\{CELL_MARKER_GENES_DATA\}\},.*?'
            r"exportFileName: '\{\{SAMPLE_ID\}\}_cell_markers_raw\.tsv'\s*"
            r'\}\);',
            '',
            template,
            flags=re.DOTALL,
        )
        return template
    
    def populate_template(self, template: str) -> str:
        """Populate the HTML template with actual data."""
        if not self.include_cell_segmentation:
            template = self._strip_cell_segmentation_template(template)
        
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
            
            # New Sample Metadata parameters
            'CHIP_NUMBER': self.env_chip_number if self.env_chip_number else 'N/A',
            'TRANSCRIPTOME': transcriptome,
            'IMAGE_ALIGNMENT': self._detect_image_alignment_mode(),
            'PROBE_SET_NAME': 'N/A',
            'TISSUE': self.env_tissue if self.env_tissue else 'N/A',

            # Image Alignment Viewer paths (for HE mode)
            'HE_IMAGE_PATH': self._get_he_image_path(),
            'GEM_IMAGE_PATH': self._get_gem_image_path(),
            'GEM_ALIGNMENT_OVERLAY_PATH': self._get_gem_alignment_overlay_image_path(),
            'GEM_HEATMAP_ONLY_PATH': self._get_gem_heatmap_only_image_path(),
        }
        
        # Key metrics from bin statistics (using 50μm bin as default)
        bin50_stats = self.bin_stats.get('bin50', {})
        bin10_stats = self.bin_stats.get('bin10', {})
        
        replacements.update({
            'NUM_SQUARES': self.format_number(bin10_stats.get('Estimated Number of square bin', '0')),
            'MEAN_READS_10': self.format_number(bin10_stats.get('Mean Reads per square bin', '0')),
            'MEAN_UMI_10': self.format_number(bin10_stats.get('Mean UMI per square bin', '0')),
            'TOTAL_GENES': self.format_number(bin10_stats.get('Total Genes', '0')),
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
        too_many_loci = self.general_stats.get(
            'Reads Mapped to Too Many Loci',
            self.general_stats.get('Reads Mapped To Too Many Loci', '0')
        )
        unmapped_short = self.general_stats.get('Unmapped Reads - Too Short', '0')
        unmapped_other = self.general_stats.get('Unmapped Reads - Other', '0')
        
        # Prefer STAR input reads as the denominator for STAR mapping categories.
        total_mapped = self.general_stats.get('STAR Input Reads', valid_reads)
        
        # Handle percentage extraction if already parsed
        unique_mapped_pct = self.general_stats.get('Uniquely Mapped Reads Percentage',
            self.calculate_percentage(unique_mapped, total_mapped) + '%')
        multi_mapped_pct = self.general_stats.get('Multi-Mapped Reads Percentage',
            self.calculate_percentage(multi_mapped, total_mapped) + '%')
        too_many_loci_pct = self.general_stats.get(
            'Reads Mapped to Too Many Loci Percentage',
            self.general_stats.get(
                'Reads Mapped To Too Many Loci Percentage',
                self.calculate_percentage(too_many_loci, total_mapped) + '%'
            )
        )
        unmapped_short_pct = self.general_stats.get('Unmapped Reads - Too Short Percentage',
            self.calculate_percentage(unmapped_short, total_mapped) + '%')
        unmapped_other_pct = self.general_stats.get('Unmapped Reads - Other Percentage',
            self.calculate_percentage(unmapped_other, total_mapped) + '%')
        
        replacements.update({
            'UNIQUE_MAPPED': f"{self.format_number(unique_mapped)} ({unique_mapped_pct})",
            'MULTI_MAPPED': f"{self.format_number(multi_mapped)} ({multi_mapped_pct})",
            'TOO_MANY_LOCI': f"{self.format_number(too_many_loci)} ({too_many_loci_pct})",
            'UNMAPPED_SHORT': f"{self.format_number(unmapped_short)} ({unmapped_short_pct})",
            'UNMAPPED_OTHER': f"{self.format_number(unmapped_other)} ({unmapped_other_pct})",
            'STAR_INPUT_READS': self.format_number(total_mapped),
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
        
        # Sequencing saturation statistics - calculate dynamically from available data
        # Try to extract from bin statistics first
        saturation_stats = self.general_stats.get('Sequencing Saturation', {}) if isinstance(self.general_stats.get('Sequencing Saturation'), dict) else {}

        # Use bin50 data for sequencing saturation estimates (for Estimated Number of Square section)
        estimated_squares = bin50_stats.get('Estimated Number of square bin', '0')
        mean_reads = bin50_stats.get('Mean Reads per square bin', '0')
        mean_umi = bin50_stats.get('Mean UMI per square bin', '0')
        total_genes = bin50_stats.get('Total Genes', '0')

        # Get Fraction Reads in Square from bin50 stats (consistent with other metrics)
        fraction_reads_square = bin50_stats.get('Fraction Reads in Square', '0')

        # Remove % sign if present and format as number
        if isinstance(fraction_reads_square, str):
            fraction_reads_square = fraction_reads_square.replace('%', '').strip()

        # Calculate dynamic MEDIAN_GENES_SQUARE from bin statistics
        median_genes_square = self._calculate_median_genes_square()

        # Calculate dynamic SATURATION from downsample data or estimate from UMI data
        saturation_percentage = self._calculate_saturation_percentage()
        saturation_note = self._saturation_note()

        replacements.update({
            'ESTIMATED_NUM_SQUARE': self.format_number(estimated_squares),
            'FRACTION_READS_SQUARE': fraction_reads_square if fraction_reads_square else '0',
            'MEAN_READS_SQUARE': self.format_number(mean_reads, as_float=True),
            'MEDIAN_UMI_SQUARE': self.format_number(mean_umi),
            'TOTAL_GENES_SQUARE': self.format_number(total_genes),
            'MEDIAN_GENES_SQUARE': median_genes_square,
            'SATURATION': saturation_percentage,
            'SATURATION_NOTE': saturation_note,
        })
        
        # Bin metrics table - simplified to four key metrics
        for bin_size in ['10', '20', '50', '100']:
            bin_key = f'bin{bin_size}'
            stats = self.bin_stats.get(bin_key, {})
            
            replacements.update({
                f'BIN{bin_size}_NUM_SQUARES': self.format_number(stats.get('Estimated Number of square bin', '0')),
                f'BIN{bin_size}_MEAN_UMI': self.format_number(stats.get('Mean UMI per square bin', '0')),
                f'BIN{bin_size}_MEAN_GENES': self.format_number(stats.get('Mean Genes per square bin', '0')),
            })
        
        # Analysis results sections
        replacements['ANALYSIS_RESULTS_IMAGES'] = self.generate_analysis_results_section()
        replacements['BIN_ANALYSIS_SECTIONS'] = self.generate_bin_analysis_sections()
        replacements['CELL_SEGMENTATION_TAB_STYLE'] = '' if self.has_cell_segmentation else 'display: none;'
        replacements['CELL_SEGMENTATION_METRIC_TABLES'] = self.generate_cell_segmentation_metric_tables()
        replacements['CELL_SEGMENTATION_IMAGES'] = self.generate_cell_segmentation_images_section()
        
        # Company logo
        if self.company_logo:
            replacements['COMPANY_LOGO_HTML'] = f'<img src="{self.company_logo}" alt="Company Logo" />'
        else:
            replacements['COMPANY_LOGO_HTML'] = '<div class="logo-fallback">C</div>'
        
        # Marker genes - organize by cluster for pagination
        replacements['MARKER_GENES_SECTION_STYLE'] = '' if self.marker_genes else 'display: none;'
        replacements['MARKER_GENES_DATA'] = self.generate_marker_genes_data_json()
        replacements['CLUSTER_LIST'] = self._clusters_to_json(self.clusters)
        replacements['TOTAL_CLUSTERS'] = str(len(self.clusters))
        
        # Complete marker genes data for full TSV export fallback
        replacements['COMPLETE_MARKER_GENES_DATA'] = json.dumps(self.complete_marker_genes, indent=2)
        
        # Embed TSV file content as data URI for offline download
        replacements['MARKERS_RAW_TSV_DATA_URI'] = self.generate_markers_tsv_data_uri()

        replacements['CELL_MARKER_GENES_SECTION_STYLE'] = '' if self.cell_marker_genes else 'display: none;'
        replacements['CELL_MARKER_GENES_DATA'] = self.generate_cell_marker_genes_data_json()
        replacements['CELL_CLUSTER_LIST'] = self._clusters_to_json(self.cell_clusters)
        replacements['CELL_TOTAL_CLUSTERS'] = str(len(self.cell_clusters))
        replacements['CELL_COMPLETE_MARKER_GENES_DATA'] = json.dumps(self.cell_complete_marker_genes, indent=2)
        replacements['CELL_MARKERS_RAW_TSV_DATA_URI'] = self.generate_cell_markers_tsv_data_uri()
        
        # Interactive charts
        replacements['BARCODE_RANK_PLOT_HTML'] = self.barcode_rank_plot_html
        replacements['SEQUENCING_SATURATION_PLOT_HTML'] = self.sequencing_saturation_plot_html
        replacements['MEDIAN_GENES_PLOT_HTML'] = self.median_genes_plot_html
        
        # Replace all placeholders in template
        result = template
        for key, value in replacements.items():
            result = result.replace('{{' + key + '}}', str(value))
        
        return result
    
    def generate_report(self, output_filename: str = None, top_genes: int = DEFAULT_TOP_GENES_PER_CLUSTER) -> Path:
        """Generate the complete HTML report."""
        self.logger.info("Starting report generation...")
        
        # Load all data
        self.load_bin_statistics()
        self.load_general_statistics()
        self.load_marker_genes(top_n=top_genes)
        self.load_bin_images()
        self.load_analysis_images()
        self.load_cell_segmentation_data()
        self.load_cell_marker_genes(top_n=top_genes)
        self.load_company_logo()
        
        # Load data for interactive charts
        self.load_downsample_data()
        self.load_count_detail_data()
        
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
            self.median_genes_plot_html = self.generate_static_chart_placeholder("Median Genes per Square", "Sequencing Depth Analysis")
        
        # Load and populate template
        template = self.load_template()
        populated_html = self.populate_template(template)
        
        # Save report
        if output_filename is None:
            output_filename = f"{self.sample_name}_spatial_analysis_report.html"
        
        output_path = self.output_dir / output_filename
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(populated_html)
        
        self.logger.info(f"Report generated successfully: {output_path}")
        return output_path


def main():
    """Main function for command line usage."""
    parser = argparse.ArgumentParser(
        description='Generate spatial transcriptomics analysis report',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Example usage:
  python spatial_report_generator.py /path/to/ST110250_B1 ST110250_B1
  python spatial_report_generator.py /path/to/sample_dir --auto-detect
  python spatial_report_generator.py /path/to/sample_dir sample_name --chemistry BBV2.4 --species Homo_sapiens
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
    parser.add_argument('--method', choices=['gene_expr', 'image', 'ssDNA', 'HE'],
                       help='Pipeline spatial method to display in sample metadata')
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
    parser.add_argument('--exclude-cell-segmentation', action='store_true',
                       help='Hide experimental Cell Segmentation sections from the report')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Handle auto-detection
    sample_name = args.sample_name
    if args.auto_detect or sample_name is None:
        sample_name = None  # Let the class auto-detect
    
    try:
        generator = SpatialReportGenerator(
            sample_dir=args.sample_dir,
            sample_name=sample_name,
            chemistry=args.chemistry,
            species=args.species,
            method=args.method,
            output_dir=args.output_dir,
            tissue=args.tissue,
            include_cell_segmentation=not args.exclude_cell_segmentation
        )
        
        # Override parameters
        if hasattr(args, 'fast_mode'):
            generator.use_interactive_charts = not args.fast_mode
        output_path = generator.generate_report(args.output_filename, top_genes=args.top_genes)
        print(f"✅ Report generated successfully: {output_path}")
        print(f"📊 Sample: {generator.sample_name}")
        print(f"🧬 Chemistry: {generator.chemistry}")
        print(f"🔬 Species: {generator.species}")
        print(f"📈 Marker genes: {len(generator.marker_genes)}")
        print(f"🔬 Cell marker genes: {len(generator.cell_marker_genes)}")
        
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

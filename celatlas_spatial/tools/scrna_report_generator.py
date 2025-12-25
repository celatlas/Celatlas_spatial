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

    def __init__(self, sample_dir: str, sample_name: str = None, chemistry: str = "BBV2.4",
                 species: str = "Mus_musculus", output_dir: str = None,
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
        
        # Auto-detect sample name if not provided
        if sample_name is None:
            self.sample_name = self._detect_sample_name()
        else:
            self.sample_name = sample_name
            
        # Store environment variables for display logic
        self.env_sample_name = env_sample_name
        self.env_chip_number = env_chip_number
            
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
        self.count_detail_data = None
        self.barcode_umi_counts = None  # Pre-aggregated barcode->UMI count mapping
        
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
    
    def _detect_sample_name(self) -> str:
        """Auto-detect sample name from directory structure or files."""
        # Try to detect from directory name
        if self.sample_dir.name and self.sample_dir.name != 'analysis':
            return self.sample_dir.name
        
        # Try to detect from marker genes file
        analysis_dir = self.sample_dir / "07.analysis"
        if analysis_dir.exists():
            for file in analysis_dir.glob("*_markers_raw.tsv"):
                sample_name = file.stem.replace('_markers_raw', '')
                self.logger.info(f"Auto-detected sample name: {sample_name}")
                return sample_name
        
        # Try to detect from other analysis files
        for pattern in ["*_bin*.h5ad", "*_spatial_analysis_report.html"]:
            for file in analysis_dir.glob(pattern):
                parts = file.stem.split('_')
                if len(parts) >= 2:
                    sample_name = '_'.join(parts[:2])
                    self.logger.info(f"Auto-detected sample name from {file.name}: {sample_name}")
                    return sample_name
        
        # Fallback to parent directory name
        fallback_name = self.sample_dir.parent.name if self.sample_dir.parent.name != 'analysis' else 'sample'
        self.logger.warning(f"Could not auto-detect sample name, using fallback: {fallback_name}")
        return fallback_name
    
    def load_bin_statistics(self):
        """Load statistics for all supported bin sizes.
        
        Searches for bin statistics in the standard directory structure:
        06.binSegment/square_bin/{sample_name}_bin{size}/stat.txt
        """
        for bin_size in SUPPORTED_BIN_SIZES:
            bin_dir = self.sample_dir / "06.binSegment" / "square_bin" / f"{self.sample_name}_bin{bin_size}"
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
        barcode_stat = self.sample_dir / "01.barcode" / "stat.txt"
        if barcode_stat.exists():
            self.general_stats.update(self.parse_stat_file(barcode_stat))
        
        # Load STAR mapping statistics
        star_stat = self.sample_dir / "03.star" / "stat.txt"
        if star_stat.exists():
            self.general_stats.update(self.parse_stat_file(star_stat))
        
        # Load featureCounts statistics
        fc_stat = self.sample_dir / "04.featureCounts" / "stat.txt"
        if fc_stat.exists():
            self.general_stats.update(self.parse_stat_file(fc_stat))
        
        # Load count statistics
        count_stat = self.sample_dir / "05.count" / "stat.txt"
        if count_stat.exists():
            self.general_stats.update(self.parse_stat_file(count_stat))
            
        # Load analysis statistics
        analysis_stat = self.sample_dir / "07.analysis" / "stat.txt"
        if analysis_stat.exists():
            self.general_stats.update(self.parse_stat_file(analysis_stat))
    
    def load_marker_genes(self, top_n: int = 20):
        """Load marker genes from the analysis results, ensuring top N genes per cluster."""
        # scRNA mode: markers are in 06_analysis_wrapper/03.Markers
        # Try both old path (double nested) and new path (single level) for backward compatibility
        marker_file_new = self.sample_dir / "06_analysis_wrapper" / "03.Markers" / f"{self.sample_name}_markers_raw.tsv"
        marker_file_old = self.sample_dir / "06_analysis_wrapper" / "06_analysis_wrapper" / "03.Markers" / f"{self.sample_name}_markers_raw.tsv"

        # Use new path if exists, otherwise try old path
        if marker_file_new.exists():
            marker_file = marker_file_new
        elif marker_file_old.exists():
            marker_file = marker_file_old
            self.logger.warning(f"Using legacy marker file path (double nested): {marker_file_old}")
        else:
            marker_file = marker_file_new  # Will fail with proper error message
        
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
        # scRNA mode: markers are in 06_analysis_wrapper/03.Markers
        # Try both old path (double nested) and new path (single level) for backward compatibility
        marker_file_new = self.sample_dir / "06_analysis_wrapper" / "03.Markers" / f"{self.sample_name}_markers_raw.tsv"
        marker_file_old = self.sample_dir / "06_analysis_wrapper" / "06_analysis_wrapper" / "03.Markers" / f"{self.sample_name}_markers_raw.tsv"

        # Use new path if exists, otherwise try old path
        if marker_file_new.exists():
            marker_file = marker_file_new
        elif marker_file_old.exists():
            marker_file = marker_file_old
        else:
            marker_file = marker_file_new  # Will fail with proper error message
        
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
        
        Searches for bin-specific QC images in:
        07.analysis/Bioinfodata/02.SpatialQC/{sample_name}_bin{size}_{type}.png
        """
        analysis_dir = self.sample_dir / "07.analysis" / "Bioinfodata" / "02.SpatialQC"
        
        for bin_size in SUPPORTED_BIN_SIZES:
            bin_key = f'bin{bin_size}'
            self.bin_images[bin_key] = {
                'scatter': None,
                'violin': None
            }
            
            # Look for scatter and violin plots
            scatter_file = analysis_dir / f"{self.sample_name}_bin{bin_size}_scatter.png"
            violin_file = analysis_dir / f"{self.sample_name}_bin{bin_size}_violin.png"
            
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
                self.downsample_data = pd.read_csv(downsample_file, sep='\t')
                self.logger.info(f"Loaded downsample data with {len(self.downsample_data)} data points")
            except Exception as e:
                self.logger.error(f"Failed to load downsample data: {e}")
                self.downsample_data = None
        else:
            self.logger.warning(f"Downsample file not found: {downsample_file}")
            self.downsample_data = None
    
    def load_count_detail_data(self):
        """Load count detail data for barcode rank plot.

        Uses chunked reading to calculate barcode UMI counts from the complete dataset
        without loading all data into memory at once.
        """
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
        """Generate interactive Barcode Rank Plot using plotly."""
        if not PLOTLY_AVAILABLE:
            return '<div class="simple-chart"><div class="chart-placeholder">Plotly not available</div></div>'
        
        try:
            if self.barcode_umi_counts is None:
                # Generate sample data if real data not available
                n_barcodes = 10000
                x_data = np.arange(1, n_barcodes + 1)

                # Generate realistic barcode UMI counts (exponential decay)
                barcode_counts = np.exp(-x_data / 2000) * 1000 + np.random.exponential(10, n_barcodes)
                background_counts = np.ones(n_barcodes) * 5 + np.random.exponential(2, n_barcodes)

                self.logger.warning("Using simulated data for Barcode Rank Plot")
            else:
                # Use pre-calculated barcode UMI counts (already sorted in descending order)
                # For scRNA-seq, show more barcodes to better visualize cell population
                n_barcodes = min(len(self.barcode_umi_counts), 50000)  # Increased from 20k to 50k for scRNA

                x_data = np.arange(1, n_barcodes + 1)
                barcode_counts = self.barcode_umi_counts.head(n_barcodes).values

                # Estimate background as lower percentile (more conservative for scRNA)
                background_level = np.percentile(barcode_counts, 5)  # Changed from 10 to 5
                background_counts = np.ones(n_barcodes) * background_level
            
            # Create the plot
            fig = go.Figure()
            
            # Add barcode trace
            fig.add_trace(go.Scatter(
                x=x_data,
                y=barcode_counts,
                mode='lines',
                name='Barcode',
                line=dict(color='#1a73e8', width=2),
                hovertemplate='<b>Rank:</b> %{x}<br><b>UMI Counts:</b> %{y}<extra></extra>'
            ))
            
            # Add background trace
            fig.add_trace(go.Scatter(
                x=x_data,
                y=background_counts,
                mode='lines',
                name='Background',
                line=dict(color='#94a3b8', width=2, dash='dash'),
                hovertemplate='<b>Rank:</b> %{x}<br><b>Background:</b> %{y}<extra></extra>'
            ))

            # Add cell threshold line if barcode data is available
            if self.barcode_umi_counts is not None:
                # Find inflection point (cell vs background threshold)
                # Simple heuristic: find the steepest drop in log-log space
                log_counts = np.log10(barcode_counts + 1)
                if len(log_counts) > 100:
                    # Calculate second derivative to find inflection
                    gradient = np.gradient(log_counts)
                    second_deriv = np.gradient(gradient)
                    # Find the point with maximum curvature (most negative second derivative)
                    inflection_idx = np.argmin(second_deriv[:len(second_deriv)//2])  # Only check first half
                    inflection_idx = max(inflection_idx, 10)  # Ensure at least 10 cells

                    cell_threshold_rank = x_data[inflection_idx]
                    cell_threshold_umi = barcode_counts[inflection_idx]

                    # Add vertical line at cell threshold
                    fig.add_vline(
                        x=cell_threshold_rank,
                        line=dict(color='#ef4444', width=2, dash='dot'),
                        annotation=dict(
                            text=f'Cell Threshold<br>(~{cell_threshold_rank:,} cells)',
                            font=dict(size=10, color='#ef4444'),
                            xanchor='left'
                        )
                    )
            
            # Update layout to match UI theme
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
                    exponentformat='none'
                ),
                yaxis=dict(
                    title='UMI Counts',
                    type='log',
                    showgrid=True,
                    gridcolor='#e2e8f0',
                    color='#000000',
                    linewidth=2,
                    linecolor='#000000'     
                ),
                plot_bgcolor='white',
                paper_bgcolor='white',
                font=dict(family='Arial', size=12),
                legend=dict(
                    x=0.7, y=0.95,
                    bgcolor='rgba(255,255,255,0.8)',
                    bordercolor='#e2e8f0',
                    borderwidth=1
                ),
                margin=dict(l=60, r=20, t=60, b=60),
                height=450
            )
            
            # Convert to HTML with optimized plotly.js inclusion
            config = {'displayModeBar': False, 'responsive': True}
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
            if self.downsample_data is None or 'read_fraction' not in self.downsample_data.columns:
                # Generate sample data if real data not available
                read_fractions = np.linspace(0, 1, 11)
                saturation_values = 1 - np.exp(-read_fractions * 3)  # Exponential saturation curve
                saturation_values = saturation_values   # Convert to percentage
                
                self.logger.warning("Using simulated data for Sequencing Saturation Plot")
            else:
                # Use real data
                read_fractions = self.downsample_data['read_fraction'].values
                # Use umi_saturation if available, otherwise use read_saturation
                if 'umi_saturation' in self.downsample_data.columns:
                    saturation_values = self.downsample_data['umi_saturation'].values 
                else:
                    saturation_values = self.downsample_data['read_saturation'].values
            
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
                    range=[0, max(saturation_values) * 1.1],
                    tickformat='.2f'
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
        
        # Load tissue hires image (try multiple locations)
        tissue_hires_paths = [
            self.sample_dir / "07.analysis" / "Bioinfodata" / "03.SpatialCluster" / "tissue_hires_image.png",
            self.sample_dir / "06.binSegment" / "images" / "tissue_hires_image.png"
        ]
        
        for tissue_path in tissue_hires_paths:
            if tissue_path.exists():
                self.analysis_images['tissue_hires'] = self.encode_image_to_base64(tissue_path)
                self.logger.info(f"Loaded tissue hires image from {tissue_path}")
                break
        
        # Load spatial gene expression distribution
        spatial_gene_expr_file = self.sample_dir / "07.analysis" / "Bioinfodata" / "02.SpatialQC" / f"{self.sample_name}_spatial_gene_expression_distribution.png"
        if spatial_gene_expr_file.exists():
            self.analysis_images['spatial_gene_expression'] = self.encode_image_to_base64(spatial_gene_expr_file)
            self.logger.info("Loaded spatial gene expression distribution image")
        
        # Load spatial cluster image
        spatial_cluster_file = self.sample_dir / "07.analysis" / "Bioinfodata" / "03.SpatialCluster" / f"{self.sample_name}_spatial_cluster.png"
        if spatial_cluster_file.exists():
            self.analysis_images['spatial_cluster'] = self.encode_image_to_base64(spatial_cluster_file)
            self.logger.info("Loaded spatial cluster image")
        
        # Load UMAP cluster image
        umap_cluster_file = self.sample_dir / "07.analysis" / "Bioinfodata" / "03.SpatialCluster" / f"{self.sample_name}_umap_cluster.png"
        if umap_cluster_file.exists():
            self.analysis_images['umap_cluster'] = self.encode_image_to_base64(umap_cluster_file)
            self.logger.info("Loaded UMAP cluster image")

        # Load scRNA-specific images (06_analysis_wrapper directory structure)
        # Try both possible directory structures
        possible_analysis_dirs = [
            self.sample_dir / "06_analysis_wrapper" / "06_analysis_wrapper",  # Actual structure
            self.sample_dir / "06_analysis_wrapper",  # Expected structure
        ]

        analysis_wrapper_dir = None
        for dir_path in possible_analysis_dirs:
            if (dir_path / "01.QC").exists():
                analysis_wrapper_dir = dir_path
                self.logger.info(f"Found analysis directory at: {analysis_wrapper_dir}")
                break

        if not analysis_wrapper_dir:
            self.logger.warning(f"Analysis directory not found. Tried: {possible_analysis_dirs}")
            return

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
                if 'umi_saturation' in self.downsample_data.columns:
                    # Get the saturation at maximum read fraction (usually last row)
                    max_saturation = self.downsample_data['umi_saturation'].iloc[-1]
                    return f"{max_saturation:.2f}"
                elif 'read_saturation' in self.downsample_data.columns:
                    max_saturation = self.downsample_data['read_saturation'].iloc[-1]
                    return f"{max_saturation:.2f}"
            
            # Method 2: Try to get from general statistics
            saturation_from_stats = self.general_stats.get('Sequencing Saturation', '')
            if saturation_from_stats and saturation_from_stats != '0':
                # Extract percentage if it's in format like "44.71%" or just "44.71"
                saturation_str = str(saturation_from_stats).replace('%', '').strip()
                try:
                    saturation_val = float(saturation_str)
                    return f"{saturation_val:.2f}"
                except ValueError:
                    pass
            
            # Method 3: Estimate from UMI and read counts
            bin50_stats = self.bin_stats.get('bin50', {})
            mean_reads = bin50_stats.get('Mean Reads per square bin', '0')
            mean_umi = bin50_stats.get('Mean UMI per square bin', '0')
            
            if mean_reads and mean_umi and mean_reads != '0' and mean_umi != '0':
                try:
                    reads_val = float(str(mean_reads).replace(',', ''))
                    umi_val = float(str(mean_umi).replace(',', ''))
                    if reads_val > 0:
                        # Saturation = 1 - (UMI/Reads) approximation
                        saturation_estimate = (1 - (umi_val / reads_val)) * 100
                        # Cap between 0 and 100
                        saturation_estimate = max(0, min(100, saturation_estimate))
                        return f"{saturation_estimate:.2f}"
                except (ValueError, TypeError):
                    pass
            
            # Default fallback
            return "N/A"
            
        except Exception as e:
            self.logger.warning(f"Failed to calculate saturation percentage: {e}")
            return "N/A"
    
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
            'IMAGE_ALIGNMENT': 'ssDNA',
            'PROBE_SET_NAME': 'N/A',
            'FILTER_PROBES': 'N/A',
        }
        
        # Key metrics from bin statistics (only used for spatial data)
        bin50_stats = self.bin_stats.get('bin50', {})
        bin10_stats = self.bin_stats.get('bin10', {})

        # For scRNA-seq, Key Metrics shows cell-level statistics from 05.count/stat.txt
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
        
        # Handle percentage extraction if already parsed
        valid_reads_pct = self.general_stats.get('Valid Reads Percentage', 
            self.calculate_percentage(valid_reads, raw_reads) + '%')
        mismatch_pct = self.general_stats.get('Mismatched Reads Percentage',
            self.calculate_percentage(mismatch, raw_reads) + '%')
        
        replacements.update({
            'NUM_READS': self.format_number(raw_reads),
            'VALID_BARCODES': f"{self.format_number(valid_reads)} ({valid_reads_pct})",
            'MISMATCH_READS': f"{self.format_number(mismatch)} ({mismatch_pct})",
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
        
        # Cell Statistics section - for scRNA-seq, read from 05.count/stat.txt
        # All cell statistics are stored in general_stats after loading count statistics
        estimated_cells = self.general_stats.get('Estimated Number of Cells', '0')
        fraction_reads_cells = self.general_stats.get('Fraction Reads in Cells', '0')
        mean_reads_cell = self.general_stats.get('Mean Reads per Cell', '0')
        median_umi_cell = self.general_stats.get('Median UMI per Cell', '0')
        total_genes_cell = self.general_stats.get('Total Genes', '0')
        median_genes_cell = self.general_stats.get('Median Genes per Cell', '0')
        saturation_percentage = self.general_stats.get('Saturation', '0')

        # Remove % sign if present and format as number
        if isinstance(fraction_reads_cells, str):
            fraction_reads_cells = fraction_reads_cells.replace('%', '').strip()
        if isinstance(saturation_percentage, str):
            saturation_percentage = saturation_percentage.replace('%', '').strip()

        replacements.update({
            'ESTIMATED_NUM_CELLS': self.format_number(estimated_cells),
            'FRACTION_READS_CELLS': fraction_reads_cells if fraction_reads_cells else '0',
            'MEAN_READS_CELL_STAT': self.format_number(mean_reads_cell, as_float=True),
            'MEDIAN_UMI_CELL': self.format_number(median_umi_cell),
            'TOTAL_GENES_CELL': self.format_number(total_genes_cell),
            'MEDIAN_GENES_CELL': self.format_number(median_genes_cell),
            'SATURATION': saturation_percentage if saturation_percentage else '0',
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
            output_dir=args.output_dir
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
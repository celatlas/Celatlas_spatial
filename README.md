# Celatlas Spatial - Spatial Transcriptomics Analysis Pipeline

[![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Celatlas Spatial is a comprehensive spatial transcriptomics analysis pipeline that provides end-to-end processing from raw sequencing data to publication-ready results.

## Installation & Setup

### Installation Methods

Choose the installation method that best fits your situation:

```bash
# Option 1: From distribution package (if available)
pip install dist/celatlas_spatial-1.6.0-py3-none-any.whl

# Option 2: From PyPI (recommended for most users)
pip install celatlas-spatial

# Option 3: From source (for development)
git clone 
cd celatlas-spatial
pip install -e .
```

### Verification

Verify installation success:

```bash
celatlas_spatial --help
```

### Requirements

- **Python**: = 3.9
- **External tools**: STAR, SAMtools (install via conda)
- **Python packages**: Automatically installed (PyTorch, Pandas, NumPy, Plotly, Jinja2)

## Prerequisites

### Must install...

**Important:** You must install the STAR aligner and other necessary software separately before using Celatlas Spatial:

```bash
conda install --file conda_pkgs.txt
```

### Reference Genome Setup

Before running analysis, you need to build reference genome indices:

**1. Create reference directory structure:**
```bash
mkdir -p reference/Homo_sapiens
mkdir -p reference/Mus_musculus
```

**2. Download genome files (example for human):** 

Note:​​ For detailed parameters and the latest guidelines, please refer to the official STAR documentation at https://github.com/alexdobin/STAR."

```bash
cd reference/Homo_sapiens
# Download genome FASTA
wget http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download GTF annotation
wget http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz
```

**3. Build STAR index:**
```bash
celatlas_spatial rna mkref \
    --genome_name Homo_sapiens \
    --fasta reference/Homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --genes reference/Homo_sapiens/Homo_sapiens.GRCh38.110.gtf \
    --thread 8
```

## Quick Start

### Using the Celatlas Pipeline Script

After installation, you can run the complete spatial transcriptomics analysis pipeline using the `Celatlas.sh` script:

```bash
# Basic usage (7 parameters)
Celatlas.sh <sample_name> <chip_number> <casno> <chemistry> <species> <method> <mode>

# Alternative usage (6 parameters)
Celatlas.sh <chip_number> <casno> <chemistry> <species> <method> <mode>
```

### Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `sample_name` | Sample identifier (optional in 6-parameter mode) | `Sample001` |
| `chip_number` | Chip number/identifier | `ST110001` |
| `casno` | Case number for result storage | `CAS250801` |
| `chemistry` | Chemistry protocol used | `BBV2.4` |
| `species` | Target species | `Homo_sapiens` or `Mus_musculus` |
| `method` | Analysis method | `image` or `gene_expr` |
| `mode` | Analysis mode | `strna` or `scrna` |

### Examples

```bash
# Image-based analysis
bash Celatlas.sh Sample001 ST110001 CAS250801 BBV2.4 Mus_musculus image strna

# Gene expression analysis
bash Celatlas.sh Sample001 ST110001 CAS250801 BBV2.4 Mus_musculus gene_expr strna

# Using 6-parameter format
bash Celatlas.sh ST110001 CAS250801 BBV2.4 Mus_musculus image strna
```

## Configuration

### Environment Variables

You can customize the pipeline behavior using environment variables:

| Variable | Description | Default Value |
|----------|-------------|---------------|
| `CELATLAS_WORKSPACE` | Main workspace directory | `/your/workspace/path` |
| `MAX_PARALLEL_FILES` | Maximum parallel file processing | 3 |

### Setting Custom Paths

**Option 1: Export environment variables**
```bash
export CELATLAS_WORKSPACE="/home/user/my_workspace"
bash Celatlas.sh Sample001 ST110001 CAS250801 BBV2.4 Mus_musculus image strna
```

**Option 2: One-time usage**
```bash
CELATLAS_WORKSPACE="/home/user/workspace" \
bash Celatlas.sh MySample001 Chip001 CAS250801 strnaV2 human image strna
```

**Option 3: Create a configuration script**
```bash
# create_config.sh
#!/bin/bash
export CELATLAS_WORKSPACE="/your/workspace/path"
export MAX_PARALLEL_FILES=4

# Usage: source create_config.sh && Celatlas.sh ...
```

## Directory Structure

The pipeline expects and creates the following directory structure:

```
$CELATLAS_WORKSPACE/
├── binSegment/           # Segmentation results
├── images/               # Input images
├── fastq/               # FASTQ files organized by chemistry
│   ├── BBV2.4
│   └── strnaV3/
├── reference/           # Reference genomes
│   ├── Homo_sapiens/
│   └── Mus_musculus/
│   
├── src/                 # Source files and models
├── rawdata/            # Raw data per sample
└── results/            # Analysis results
    └── {casno}/        # Organized by case number
        └── {sample}/   # Sample-specific results
```

## Pipeline Steps

The Celatlas.sh script runs the following analysis steps:

1. **Sample Processing** - Initial sample validation and setup
2. **Barcode Extraction** - Spatial barcode identification
3. **Adapter Trimming** - Quality control and adapter removal
4. **Sequence Alignment** - STAR alignment to reference genome
5. **Feature Counting** - Gene expression quantification
6. **Cell Calling** - Cell detection and filtering
7. **Spatial Binning** - Image segmentation and spatial binning
8. **Analysis & Visualization** - Comprehensive spatial analysis
9. **Report Generation** - HTML report with interactive plots

## Output Files

Key output files include:

- `{sample}_spatial_analysis_report.html` - Comprehensive analysis report
- `square_bin/` - Spatial bin matrices
- `pipeline.log` - Detailed execution log
- Various intermediate files for each processing step

## Command Line Tools

Individual pipeline components are also available as command-line tools:

```bash
# RNA analysis subcommands
celatlas_spatial rna sample --help
celatlas_spatial rna barcode --help
celatlas_spatial rna cutadapt --help
celatlas_spatial rna star --help
celatlas_spatial rna featureCounts --help
celatlas_spatial rna count --help
celatlas_spatial rna binSegment --help
celatlas_spatial rna analysis --help
```

## Additional Notes

All software dependencies and system requirements are detailed in the [Installation & Setup](#installation--setup) section above. Python packages are automatically installed with pip, but bioinformatics tools like STAR and SAMtools need to be installed separately as described in the Prerequisites section.

## Support

For issues and questions:
- Create an issue on GitHub
- Contact: `rd@celatlas.com` or `lqs60667106@gmail.com`

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use Celatlas Spatial in your research, please cite:

```
https://github.com/
```

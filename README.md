# Celatlas Spatial - Spatial Transcriptomics Analysis Pipeline

[![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.7.0-brightgreen.svg)](https://github.com/)

Celatlas Spatial is a comprehensive spatial transcriptomics analysis pipeline that provides end-to-end processing from raw sequencing data to publication-ready results. Version 1.7.0 introduces major enhancements including **H&E staining-guided analysis (HE mode)**, flexible parameter configuration, and optimized multi-threading support.

## Installation & Setup

### Installation Methods

Choose the installation method that best fits your situation:

```bash
# Option 1: Download and unzip the package, then install
# 1. Download the zip file from the release page
# 2. Navigate to the directory where you downloaded the file
# 3. Unzip the file to a folder
unzip celatlas_spatial-mian.zip
cd celatlas_spatial-mian
pip install dist/celatlas_spatial-1.7.0-py3-none-any.whl

# Option 2: git clone
git clone [repository_url]
cd Celatlas_spatial
pip install dist/celatlas_spatial-1.7.0-py3-none-any.whl
```

### Download Deep Learning Model (Required)

**IMPORTANT: Download the Swin-UNET model before running analysis:**

```bash
# 1. Download swin_tiny_model.tar.gz (97MB) from GitHub Release

# 2. Extract to workspace src directory
mkdir -p /mnt/strna/celatlas_spatial/src
tar -xzf swin_tiny_model.tar.gz -C /mnt/strna/celatlas_spatial/src/

# If using custom workspace:
# mkdir -p $CELATLAS_WORKSPACE/src
# tar -xzf swin_tiny_model.tar.gz -C $CELATLAS_WORKSPACE/src/

# 3. Verify extraction
ls -lh /mnt/strna/celatlas_spatial/src/swin_tiny.pth
# Should show ~106MB swin_tiny.pth file
```

**Model Details:**
- **Filename**: `swin_tiny.pth`
- **Size**: ~106MB
- **Location**: `$CELATLAS_WORKSPACE/src/swin_tiny.pth`
- **Purpose**: Tissue segmentation for all analysis modes

### Verification

Verify installation success:

```bash
celatlas_spatial --help
```

### Requirements

- **Python**: >= 3.9
- **External tools**: STAR, SAMtools, Subread, TRUST4 (install via conda)
- **Python packages**: Automatically installed (PyTorch, Pandas, NumPy, Plotly, Jinja2)
- **Hardware** (HE mode):
  - **GPU**: NVIDIA GPU with CUDA 12 support (recommended for deep learning tissue segmentation)
  - **CPU**: Multi-core processor (16+ cores recommended)
  - **RAM**: 64GB+ (128GB recommended for large tissue sections)
  - **Storage**: 100GB+ free space for temporary files and results

## Prerequisites

### Must install...

**Important:** You must install the STAR aligner and other necessary software separately before using Celatlas Spatial:

```bash
conda install --file conda_pkgs.txt
```

### Deep Learning Model Setup (Required)

**IMPORTANT: You must download the Swin-UNET model for tissue segmentation:**

```bash
# 1. Download swin_tiny_model.tar.gz from GitHub Release or provided link

# 2. Extract model to src directory
mkdir -p $CELATLAS_WORKSPACE/src
tar -xzf swin_tiny_model.tar.gz -C $CELATLAS_WORKSPACE/src/

# 3. Verify the model file exists
ls -lh $CELATLAS_WORKSPACE/src/swin_tiny.pth
# Should show ~106MB file
```

**Notes:**
- Model size: ~106MB
- Required for tissue segmentation in all analysis modes
- Default workspace: `/mnt/strna/celatlas_spatial` (customizable via `CELATLAS_WORKSPACE`)

### Reference Genome Setup

Before running analysis, you need to build reference genome indices:

**1. Create reference directory structure:**
```bash
mkdir -p reference/Homo_sapiens
or
mkdir -p reference/Mus_musculus
or
mkdir -p reference/others
```

**2. Download genome files (example for human):** 

Note: For detailed parameters and the latest guidelines, please refer to the official STAR documentation at https://github.com/alexdobin/STAR."

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
    --gtf reference/Homo_sapiens/Homo_sapiens.GRCh38.110.gtf \
    --thread 8
```

## Quick Start

### Using the Celatlas Pipeline Script

After installation, you can run the complete spatial transcriptomics analysis pipeline using the `Celatlas.sh` script. **Version 1.7.0 introduces flexible parameter modes** - you can use traditional positional arguments or modern named arguments, and optionally add performance tuning parameters.

```bash
# Modern usage (Recommended - with named arguments)
Celatlas.sh --sample <name> --chip <chip> --casno <case> --chemistry <chem> \
            --species <species> --method <method> --mode <mode> [OPTIONS]

# Legacy usage (Positional arguments - backward compatible)
Celatlas.sh <sample_name> <chip_number> <casno> <chemistry> <species> <method> <mode> [OPTIONS]
Celatlas.sh <chip_number> <casno> <chemistry> <species> <method> <mode> [OPTIONS]

# Hybrid mode (Recommended for flexibility)
Celatlas.sh <positional args...> --thread 32 --bin 10,50,100
```

### Required Parameters

| Parameter | Description | Valid Values | Example |
|-----------|-------------|--------------|---------|
| `sample_name` | Sample identifier (optional in 6-param mode) | Any string | `DEMO` |
| `chip_number` | Chip/slide number | Unique identifier | `ST110001_A1` |
| `casno` | Case/project number for organizing results | Project ID | `HE_TEST` |
| `chemistry` | Chemistry version of the kit | `BBV0`, `BBV2.4`, `BBV3` | `BBV2.4` |
| `species` | Target species for alignment | `Homo_sapiens`, `Mus_musculus`, `others` | `Mus_musculus` |
| `method` | Analysis method **(see detailed explanation below)** | `HE`, `ssDNA`, `gene_expr` | `HE` |
| `mode` | Pipeline mode | `strna`, `scrna` | `strna` |

### Analysis Methods - Detailed Explanation

**Understanding the `--method` parameter is crucial for successful analysis:**

| Method | Full Name | When to Use | Required Files | Key Features |
|--------|-----------|-------------|----------------|--------------|
| **`HE`** | H&E Staining Guided Analysis | **Most common** - When you have H&E stained histology images | • FASTQ files<br>• `{chip}_he.(tif\|png\|jpg)`<br>• Barcode position files | ✨ **NEW in v1.7.0**<br>• Gene expression-based tissue detection<br>• HE image registration<br>• Combined visualization<br>• Optimal for histological context |
| `ssDNA` | Single-strand DNA Probes | When using fluorescent ssDNA tissue imaging | • FASTQ files<br>• `{chip}.tif` (tissue image)<br>• Barcode position files | • Tissue boundary from imaging<br>• Traditional image segmentation<br>• Alias for `image` mode |
| `gene_expr` | Gene Expression Only | When you have **NO imaging data** | • FASTQ files<br>• Barcode position files | • Pure computational tissue detection<br>• UMI-based segmentation<br>• No image required |

> **🚨 CRITICAL WARNING: HE Image File Naming**
>
> **The `_he` suffix is MANDATORY and required!**
>
> | Your Chip Number | ✅ Correct Filename | ❌ Common Mistakes |
> |------------------|---------------------|---------------------|
> | `ST110001_A1` | `ST110001_A1_he.tif` | `ST110001_A1.tif` (missing `_he`) |
| `ST110001` | `ST110001_he.png` OR `ST110001_HE.png` | `ST110001.png` (missing `_he`/`_HE`) |
| `ST110001` | `ST110001_he.jpg` OR `ST110001_HE.jpg` | `ST110001he.tif` (missing underscore) |

**Note:** Case-insensitive - both `_he` and `_HE` (and `_He`, `_hE`) are supported!
> | `ST110001` | `ST110001_he.png` | `ST110001_HE.tif` (wrong case) |
> | `ST110001` | `ST110001_he.jpg` | `ST110001he.tif` (missing underscore) |
>
> **Pipeline behavior:**
> - With `_he` suffix → HE mode (gene expression + H&E registration)
> - Without `_he` suffix → ssDNA mode (looks for `{chip}.tif` fluorescent image)
> - Wrong case/format → **ERROR: HE mode requires H&E staining image!**


**Important Notes:**
- **HE mode is the recommended method** for most applications as it combines gene expression data with histological information
- The HE image filename **must** follow the pattern: `{chip_number}_he.(tif|png|jpg|jpeg)`
- For example, if `chip_number=ST110001_A1`, the HE image should be named `ST110001_A1_he.tif` (or `.png`, `.jpg`)

### Optional Parameters (NEW in v1.7.0)

These parameters allow fine-tuning of performance and resolution:

| Parameter | Description | Default | Valid Range | Recommendation |
|-----------|-------------|---------|-------------|----------------|
| `--thread` | Number of CPU threads | Auto-detect | 1-128 | Use 80% of available cores<br>e.g., `--thread 64` on 80-core machine |
| `--bin` | Bin sizes in micrometers (comma-separated) | `50` | `10,20,50,100` | • `10`: Highest resolution (~single cell)<br>• `20`: High resolution<br>• **`50`: Balanced (recommended)**<br>• `100`: More genes, lower resolution |
| `--insertR2` | Insert R2 fragment size (bp) | `150` | 50-300 | Match your sequencing read length |
| `--cell_num` | Expected number of cells/spots | `50000` | 1000-200000 | Adjust based on tissue size |
| `--pixelSize` | Pixel size in micrometers | `0.5` | 0.1-2.0 | Check your imaging system specs |

### Path Customization (NEW in v1.7.0)

For users with shared data environments or custom storage layouts, you can now specify custom paths for reference genomes, images, and FASTQ files:

| Parameter | Description | Default | Use Case |
|-----------|-------------|---------|----------|
| `--reference_dir` | Custom reference genome directory | `$WORKSPACE/reference` | Shared reference across multiple projects |
| `--mask_dir` | Custom mask/barcode directory | `$WORKSPACE/ST_mask` | Centralized barcode storage |
| `--image_dir` | Custom image directory | `$WORKSPACE/images` | Separate image storage location |
| `--fastq_dir` | Custom FASTQ directory | `$WORKSPACE/fastq/$chemistry` | Avoid duplicate data copies, use original sequencing location |
| `--fastq_name` | Custom FASTQ filename prefix | `chip_number` | When FASTQ files have different names than chip numbers |

**Benefits:**
- ✅ **Avoid data duplication**: Point directly to original sequencing files
- ✅ **Shared resources**: Multiple users can share reference genomes and images
- ✅ **Flexible organization**: Adapt to existing storage infrastructure
- ✅ **Storage savings**: No need to copy large files into workspace

### FASTQ File Format Support (Enhanced in v1.7.0)

The pipeline now supports **4 FASTQ naming formats** with automatic detection:

| Priority | Format | Pattern | Example |
|----------|--------|---------|---------|
| 1 | Multi-lane (Recommended) | `{chip}_S*_L*_R1_*.fastq.gz` | `ST110001_A1_S1_L001_R1_001.fastq.gz`<br>`ST110001_A1_S1_L002_R1_001.fastq.gz` |
| 2 | Multi-fold (Legacy) | `{chip}_fold{1-5}_1.fq.gz` | `ST110001_A1_fold1_1.fq.gz`<br>`ST110001_A1_fold2_1.fq.gz` |
| 3 | Simple format | `{chip}_1.fq.gz` | `ST110001_A1_1.fq.gz`<br>`ST110001_A1_2.fq.gz` |
| 4 | _R1/_R2 format | `{chip}_R1.fq.gz` | `ST110001_A1_R1.fq.gz`<br>`ST110001_A1_R2.fq.gz` |

**Note:** The pipeline automatically detects and uses the first available format.

### Usage Examples

#### Example 1: Basic HE Mode (Most Common)
```bash
# Using H&E stained image for tissue analysis
# Required files: ST110001_A1_he.tif ⚠️ MUST have "_he" suffix!, FASTQ files, barcode files
bash Celatlas.sh DEMO ST110001_A1 HE_test BBV2.4 Mus_musculus HE strna
```
> **🚨 CRITICAL: HE Mode File Naming Rule**
> 
> **The HE image MUST have `_he` suffix before the file extension!**
> 
> ✅ **Correct:** `ST110001_A1_he.tif` (or `.png`, `.jpg`, `.jpeg`)  
> ❌ **Wrong:** `ST110001_A1.tif` → This will be treated as ssDNA mode!  
> ✅ **Also Correct:** `ST110001_A1_HE.tif` → Case-insensitive (both `_he` and `_HE` work!)  
> ❌ **Wrong:** `ST110001_A1he.tif` → Missing underscore before `he`
> 
> **Why this matters:** Without the `_he` suffix, the pipeline cannot detect HE mode and will fail with an error.


#### Example 2: HE Mode with Performance Optimization
```bash
# High-performance mode with 64 threads and single bin size
bash Celatlas.sh DEMO ST110001_A1 HE_test BBV2.4 Mus_musculus HE strna \
  --thread 64 --bin 50
```

#### Example 3: Multi-Resolution Analysis
```bash
# Generate multiple bin sizes for comparison (10µm, 20µm, 50µm, 100µm)
# Useful for exploring optimal spatial resolution
bash Celatlas.sh DEMO ST110001_A1 multi_res BBV2.4 Mus_musculus HE strna \
  --thread 32 --bin 10,20,50,100
```

#### Example 4: Gene Expression Only (No Image)
```bash
# When you don't have tissue images
bash Celatlas.sh DEMO ST110001_A1 expr_only BBV2.4 Mus_musculus gene_expr strna
```

#### Example 5: ssDNA Probe Mode
```bash
# When using fluorescent ssDNA tissue imaging
# Required: ST110001_A1.tif (tissue image)
bash Celatlas.sh DEMO ST110001_A1 ssDNA_test BBV2.4 Homo_sapiens ssDNA strna \
  --thread 48
```

#### Example 6: scRNA-seq Mode
```bash
# Single-cell RNA-seq analysis (no spatial information)
bash Celatlas.sh DEMO ST110001_A1 scrna_test BBV0 Homo_sapiens gene_expr scrna
```

#### Example 7: Full Custom Configuration
```bash
# Complete parameter customization
bash Celatlas.sh DEMO ST110001_A1 custom_test BBV2.4 Homo_sapiens HE strna \
  --thread 48 \
  --bin 10,20,50,100 \
  --insertR2 120 \
  --cell_num 80000 \
  --pixelSize 0.5
```

#### Example 8: Pure Named Arguments (Good for Scripts)
```bash
# Recommended for automation and reproducibility
bash Celatlas.sh \
  --sample DEMO \
  --chip ST110001_A1 \
  --casno batch_001 \
  --chemistry BBV2.4 \
  --species Mus_musculus \
  --method HE \
  --mode strna \
  --thread 32 \
  --bin 50
```

#### Example 9: Custom Directories for Shared Data Environment ✨ NEW
```bash
# Ideal for production environments with centralized storage
bash Celatlas.sh \
  --chip ST110001_A1 \
  --casno shared_proj_001 \
  --chemistry BBV2.4 \
  --species Mus_musculus \
  --method HE \
  --mode strna \
  --reference_dir /data/shared/reference \
  --fastq_dir /data/sequencing/run001 \
  --mask_dir /data/shared/masks \
  --image_dir /data/shared/images
```

**Why use custom directories?**
- Avoid copying 100+ GB FASTQ files
- Share reference genomes across multiple projects
- Use original sequencing output location directly

#### Example 10: Custom FASTQ Name (Different from Chip Number) ✨ NEW
```bash
# When FASTQ files have different naming than chip numbers
bash Celatlas.sh \
  --chip ST110001_A1 \
  --casno test_001 \
  --chemistry BBV2.4 \
  --species Mus_musculus \
  --method HE \
  --mode strna \
  --fastq_name "Sample_ABC_XYZ"

# Pipeline will look for:
#   Sample_ABC_XYZ_S1_L001_R1_001.fastq.gz (multi-lane)
#   Sample_ABC_XYZ_fold1_1.fq.gz (multi-fold)
#   Sample_ABC_XYZ_1.fq.gz (simple)
#   Sample_ABC_XYZ_R1.fq.gz (R1/R2 format)
```

## Configuration

### Sample Naming Convention (IMPORTANT - Changed in v1.7.0)

**Understanding the sample naming system is critical for file organization:**

#### Two-Part Naming System

The pipeline uses a **two-part naming system** to distinguish between sample identity and technical replicate:

1. **`sample_name`** (Optional): Biological sample identifier
   - Example: `DEMO`, `Sample001`, `MouseBrain_A`
   - Used for: Grouping related samples, metadata tracking
   - Can be omitted (will show as "N/A" in reports)

2. **`chip_number`** (Required): Chip/slide identifier (the actual data identifier)
   - Example: `ST110001_A1`, `ST110001`, `ST110001`
   - Used for: Finding FASTQ files, matching image files, creating output directories
   - **This is the primary identifier used throughout the pipeline**

#### File Naming Requirements

**All input files must use the `chip_number` as the base name:**

| File Type | Naming Pattern | Example | Notes |
|-----------|----------------|---------|-------|
| **FASTQ (Multi-lane)** | `{chip}_S*_L*_R1_*.fastq.gz`<br>`{chip}_S*_L*_R2_*.fastq.gz` | `ST110001_A1_S1_L001_R1_001.fastq.gz`<br>`ST110001_A1_S1_L002_R1_001.fastq.gz` | ✅ **Recommended format**<br>Automatically detects all lanes |
| **FASTQ (Multi-fold)** | `{chip}_fold1_1.fq.gz`<br>`{chip}_fold2_1.fq.gz` | `ST110001_A1_fold1_1.fq.gz`<br>`ST110001_A1_fold1_2.fq.gz` | Legacy format<br>Supports fold1-fold5 |
| **FASTQ (Simple)** | `{chip}_1.fq.gz`<br>`{chip}_2.fq.gz` | `ST110001_A1_1.fq.gz`<br>`ST110001_A1_2.fq.gz` | Simple format |
| **FASTQ (_R1/_R2)** | `{chip}_R1.fq.gz`<br>`{chip}_R2.fq.gz` | `ST110001_A1_R1.fq.gz`<br>`ST110001_A1_R2.fq.gz` | Common alternative format |
| **HE Image** | `{chip}_he.(tif\|png\|jpg\|jpeg)` | `ST110001_A1_he.tif` | **Required for HE mode** |
| **Tissue Image** | `{chip}.tif` | `ST110001_A1.tif` | Required for ssDNA mode |
| **Barcode Position** | `{chip}.barcodeToPos.h5`<br>`{chip}_FilterBarcodes.csv`<br>`{chip}_tissue_bbox.csv` | `ST110001_A1.barcodeToPos.h5` | Required for spatial modes |

**⚠️ IMPORTANT: HE Mode Image Naming**

The most common error in HE mode is incorrect image file naming. Please note:

- **Mandatory suffix:** The HE image MUST have `_he` before the file extension
- **Case sensitive:** Must be lowercase `_he`, not `_HE` or `_He`
- **With underscore:** Must be `_he`, not `he` (missing underscore will fail)
- **No extra text:** Must be `{chip}_he.ext`, not `{chip}_he_anything.ext`

**Examples:**
- ✅ `ST110001_A1_he.tif` → Correct
- ✅ `ST110001_he.png` → Correct
- ❌ `ST110001_A1.tif` → Will be treated as ssDNA mode
- ❌ `ST110001_A1_HE.tif` → Case mismatch, will not be detected
- ❌ `ST110001_A1he.tif` → Missing underscore, will not be detected


#### FASTQ File Detection Priority

The pipeline searches for FASTQ files in this order:

```
Priority 1: Multi-lane format (Recommended)
  ├── Pattern: {chip}_S*_L*_R1_*.fastq.gz
  ├── Example: ST110001_A1_S1_L001_R1_001.fastq.gz
  │           ST110001_A1_S1_L002_R1_001.fastq.gz
  └── Auto-detects and merges all lanes

Priority 2: Multi-fold format (Legacy)
  ├── Pattern: {chip}_fold{1-5}_1.fq.gz
  ├── Example: ST110001_A1_fold1_1.fq.gz
  │           ST110001_A1_fold2_1.fq.gz
  └── Supports up to 5 fold files

Priority 3: Single file format (Simple)
  ├── Pattern: {chip}_1.fq.gz
  ├── Example: ST110001_A1_1.fq.gz
  └── Single file pair

Priority 4: _R1/_R2 format 
  ├── Pattern: {chip}_R1.fq.gz / {chip}_R2.fq.gz
  ├── Example: ST110001_A1_R1.fq.gz
  │           ST110001_A1_R2.fq.gz
  └── Common alternative naming convention
```

**Custom FASTQ Directory & Name:**
You can override the default FASTQ location and naming using `--fastq_dir` and `--fastq_name` parameters (see Example 9 and 10 above).

#### Comparison with Previous Versions

| Version | Sample Identifier | File Base Name | Example Command |
|---------|------------------|----------------|-----------------|
| **v1.7.0** | `chip_number` (primary)<br>`sample_name` (optional metadata) | Uses `chip_number` | `Celatlas.sh DEMO ST110001_A1 ...` |
| v1.6.x | `sample` only | Uses `sample` | `Celatlas.sh Sample001 ...` |

**Key Change:** In v1.7.0, the `chip_number` parameter is now the primary identifier for file matching, while `sample_name` is optional metadata. This allows better organization when multiple samples share the same chip.

#### Practical Examples

**Example 1: Standard HE Mode with Multi-lane FASTQ**
```bash
# File structure:
# /mnt/strna/celatlas_spatial/
# ├── fastq/BBV2.4/
# │   ├── ST110001_A1_S1_L001_R1_001.fastq.gz
# │   ├── ST110001_A1_S1_L001_R2_001.fastq.gz
# │   ├── ST110001_A1_S1_L002_R1_001.fastq.gz
# │   └── ST110001_A1_S1_L002_R2_001.fastq.gz
# ├── images/
# │   └── ST110001_A1_he.tif
# └── ST_mask/
#     ├── ST110001_A1.barcodeToPos.h5
#     ├── ST110001_A1_FilterBarcodes.csv
#     └── ST110001_A1_tissue_bbox.csv

bash Celatlas.sh DEMO ST110001_A1 HE_test BBV2.4 Mus_musculus HE strna
```
> **🚨 CRITICAL: HE Mode File Naming Rule**
> 
> **The HE image MUST have `_he` suffix before the file extension!**
> 
> ✅ **Correct:** `ST110001_A1_he.tif` (or `.png`, `.jpg`, `.jpeg`)  
> ❌ **Wrong:** `ST110001_A1.tif` → This will be treated as ssDNA mode!  
> ✅ **Also Correct:** `ST110001_A1_HE.tif` → Case-insensitive (both `_he` and `_HE` work!)  
> ❌ **Wrong:** `ST110001_A1he.tif` → Missing underscore before `he`
> 
> **Why this matters:** Without the `_he` suffix, the pipeline cannot detect HE mode and will fail with an error.


**Example 2: Gene Expression Mode with Simple FASTQ**
```bash
# File structure:
# /mnt/strna/celatlas_spatial/
# ├── fastq/BBV2.4/
# │   ├── ST110001_1.fq.gz
# │   └── ST110001_2.fq.gz
# └── ST_mask/
#     ├── ST110001.barcodeToPos.h5
#     ├── ST110001_FilterBarcodes.csv
#     └── ST110001_tissue_bbox.csv

bash Celatlas.sh DEMO ST110001 test_001 BBV2.4 Homo_sapiens gene_expr strna
```

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
bash Celatlas.sh DEMO ST110001 CAS250801 BBV2.4 Mus_musculus image strna
```

**Option 2: One-time usage**
```bash
CELATLAS_WORKSPACE="/home/user/workspace" \
bash Celatlas.sh DEMO ST110001 CAS250801 strnaV2 human image strna
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
│   └── BBV2.4/
│   
├── reference/           # Reference genomes
│   ├── Homo_sapiens/
│   ├── Mus_musculus/
|   └── others/
│   
├── src/                 # Source files and models
├── rawdata/            # Raw data per sample
└── results/            # Analysis results
    └── {casno}/        # Organized by case number
        └── {sample}/   # Sample-specific results
```

## Pipeline Steps

The Celatlas.sh script runs the following analysis steps. The workflow varies depending on the selected method (`HE`, `ssDNA`, or `gene_expr`):

### Standard Workflow (All Methods)

1. **Sample Processing** (`00.sample`)
   - Sample validation and metadata collection
   - Chemistry detection and verification
   - System resource check

2. **Barcode Extraction** (`01.barcode`)
   - Spatial barcode identification from R1 reads
   - Whitelist-based barcode correction
   - UMI extraction
   - **Performance:** Multi-file parallel processing (configurable with `MAX_PARALLEL_FILES`)

3. **Adapter Trimming** (`02.cutadapt`)
   - Quality control and adapter removal from R2 reads
   - NextSeq-specific quality trimming
   - Insert size filtering based on `--insertR2` parameter

4. **Sequence Alignment** (`03.star`)
   - STAR alignment to reference genome
   - **Threading:** Utilizes `--thread` parameter for parallel alignment
   - Adaptive parameters based on mode (stricter for spatial, relaxed for scRNA)

5. **Feature Counting** (`04.featureCounts`)
   - Gene expression quantification
   - Counts reads per gene per barcode
   - **Threading:** Multi-core counting enabled

6. **UMI Counting** (`05.count`)
   - UMI deduplication
   - Cell/spot calling using EmptyDrops algorithm
   - Expected cell number based on `--cell_num` parameter
   - Generates filtered and raw count matrices

### Method-Specific Segmentation (Step 6)

The **binSegment** step varies significantly based on the selected method:

#### 🔵 HE Mode (Recommended - NEW in v1.7.0)

**`06.binSegment` - Gene Expression + H&E Image Registration**

This is the **most advanced mode** combining computational and visual tissue detection:

```
Step 6A: Gene Expression-based Tissue Detection
  ├── Aggregate UMI counts into spatial grid (default: 20µm bins)
  ├── Apply UMI threshold filtering (default: ≥30 UMIs)
  ├── Enhance signal using percentile normalization (p_low=5, p_high=95)
  ├── Generate tissue mask from gene expression heatmap
  └── Create GEM (Gene Expression Matrix) tissue boundary

Step 6B: H&E Image Registration
  ├── Load H&E stained histology image ({chip}_he.tif/png/jpg)
  ├── Tissue segmentation using Swin-UNET deep learning model
  ├── Extract tissue region from H&E image
  ├── Register H&E tissue to GEM tissue boundary
  │   ├── Registration method: Affine transformation (default)
  │   ├── Alignment: SimpleITK with mutual information metric
  │   └── Initial alignment: Contour-based initialization
  └── Generate aligned overlay visualization

Step 6C: Spatial Binning
  ├── Generate bins at specified resolutions (--bin parameter)
  │   ├── bin10/  - 10µm bins (~single cell resolution)
  │   ├── bin20/  - 20µm bins (high resolution)
  │   ├── bin50/  - 50µm bins (balanced, recommended)
  │   └── bin100/ - 100µm bins (more genes, lower resolution)
  ├── Assign barcodes to spatial bins
  ├── Aggregate UMI counts per bin
  └── Generate count matrices for each bin size
```

**Key Features:**
- **Dual tissue detection:** Uses both gene expression AND histology
- **Robust registration:** Aligns H&E image to gene expression coordinates
- **Interactive visualization:** HTML report shows HE-GEM overlay with toggle controls
- **Quality metrics:** Registration accuracy, tissue coverage statistics

**Output Files (HE Mode):**
```
06.binSegment/
├── {sample}_tissue_HE.jpg              # Segmented H&E tissue region
├── {sample}_1_tissue.png               # Gene expression tissue mask
├── {sample}_2_he_registered.jpg        # Registered H&E image
├── {sample}_3_overlay_combined.jpg     # HE + GEM overlay
├── {sample}_3c_gem_heatmap_only.png    # Pure GEM heatmap (for download)
└── square_bin/
    ├── {sample}_bin10/
    │   ├── filtered_feature_bc_matrix/
    │   ├── stat.txt
    │   └── downsample.tsv
    ├── {sample}_bin50/
    └── {sample}_bin100/
```

#### 🟢 ssDNA Mode (Image-based)

**`06.binSegment` - Fluorescent Image Segmentation**

Traditional image-based tissue detection using fluorescent ssDNA probes:

```
Step 6: Image Segmentation + Spatial Binning
  ├── Load tissue image ({chip}.tif)
  ├── Tissue segmentation using Swin-UNET model
  ├── Extract tissue boundary
  ├── Generate spatial bins (--bin parameter)
  ├── Assign barcodes to bins
  └── Create count matrices per bin
```

**Required:** High-quality fluorescent tissue image (`{chip}.tif`)

#### 🟡 gene_expr Mode (No Image)

**`06.binSegment` - Pure Computational Segmentation**

Tissue detection based solely on gene expression data:

```
Step 6: Gene Expression-based Segmentation
  ├── Aggregate UMI counts into spatial grid (--gem-bin-size, default: 20µm)
  ├── Apply UMI threshold (--umi-min-threshold, default: ≥30)
  ├── Enhance signal (--enhance-params)
  ├── Detect tissue boundary from expression pattern
  ├── Generate spatial bins (--bin parameter)
  └── Create count matrices per bin
```

**Advantages:** No imaging required, works with any spatial platform

### Analysis and Visualization (Steps 7-8)

7. **Spatial Analysis** (`07.analysis`)
   - Dimensionality reduction (UMAP, t-SNE*)
   - Clustering (Leiden algorithm)
   - Marker gene identification
   - Spatial expression patterns
   - **Performance:**
     - Uses `--thread` for parallel computation
     - **Smart t-SNE handling:** For high-resolution data (bin10, bin20), t-SNE is automatically skipped to prevent memory issues
     - UMAP is always computed for all bin sizes

8. **Report Generation**
   - Interactive HTML report with Plotly visualizations
   - **HE Mode:** Includes Image Alignment Viewer with HE-GEM overlay
   - **All Modes:** QC metrics, spatial plots, clustering results
   - Downloadable figures (PNG format)

### Performance Optimization

The pipeline includes several performance enhancements in v1.7.0:

- **Auto-threading:** Automatically detects CPU cores if `--thread` not specified
- **BLAS thread limiting:** Prevents segmentation faults on high-core servers (auto-limits to 16 threads)
- **Parallel barcode extraction:** Multiple FASTQ files processed concurrently
- **Smart t-SNE skipping:** High-resolution bins skip t-SNE to avoid memory issues
- **Intermediate file cleanup:** Optional automatic cleanup of intermediate files (controlled by `CLEAN_INTERMEDIATE` variable)

**Thread Configuration Details:**
```bash
# Auto-detected (recommended for most users)
bash Celatlas.sh ... HE strna  # Uses all available cores

# Manual specification (for high-performance servers)
bash Celatlas.sh ... HE strna --thread 64

# Conservative mode (for shared servers)
bash Celatlas.sh ... HE strna --thread 16
```

**Internal Thread Management:**
- **Analysis threads:** Set by `--thread` parameter (e.g., 64)
- **BLAS threads:** Auto-limited to 16 (prevents OpenBLAS segfault)
- **Environment variables set:**
  - `OMP_NUM_THREADS`: Analysis threads
  - `OPENBLAS_NUM_THREADS`: Limited to 16 (max)
  - `MKL_NUM_THREADS`: Limited to 16 (max)
  - `NUMEXPR_NUM_THREADS`: Analysis threads

## Output Files

The pipeline generates comprehensive outputs organized by case number and sample:

### Output Directory Structure

```
$CELATLAS_WORKSPACE/results/{casno}/{chip_number}/
├── 00.sample/
│   └── stat.json                           # Sample metadata and chemistry info
├── 01.barcode/
│   ├── {chip}_1.fq.gz                      # Barcode-corrected R1
│   ├── {chip}_2.fq.gz                      # Barcode-corrected R2
│   └── stat.json                           # Barcode extraction statistics
├── 02.cutadapt/
│   ├── {chip}_clean_2.fq.gz                # Adapter-trimmed R2
│   └── stat.json                           # Trimming statistics
├── 03.star/
│   ├── {chip}_Aligned.sortedByCoord.out.bam
│   ├── {chip}_Log.final.out                # Alignment summary
│   └── stat.json                           # STAR statistics
├── 04.featureCounts/
│   ├── {chip}_nameSorted.bam               # Name-sorted BAM for counting
│   ├── {chip}_counts.txt                   # Gene-level counts
│   └── stat.json                           # Feature counting stats
├── 05.count/
│   ├── {chip}_count_detail.txt             # Detailed UMI counts per barcode
│   ├── {chip}_filtered_feature_bc_matrix/  # Filtered count matrix (cells)
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── {chip}_raw_feature_bc_matrix/       # Raw count matrix (all barcodes)
│   └── stat.json                           # Cell calling statistics
├── 06.binSegment/
│   ├── {chip}_tissue_HE.jpg                # [HE mode] Segmented H&E tissue
│   ├── {chip}_1_tissue.png                 # [HE mode] Gene expression mask
│   ├── {chip}_2_he_registered.jpg          # [HE mode] Registered H&E
│   ├── {chip}_3_overlay_combined.jpg       # [HE mode] HE+GEM overlay
│   ├── {chip}_3c_gem_heatmap_only.png      # [HE mode] Pure GEM heatmap
│   └── square_bin/
│       ├── {chip}_bin10/                   # 10µm resolution
│       │   ├── filtered_feature_bc_matrix/
│       │   ├── stat.txt                    # Bin-specific statistics
│       │   └── downsample.tsv              # Saturation curve data
│       ├── {chip}_bin20/                   # 20µm resolution
│       ├── {chip}_bin50/                   # 50µm resolution (recommended)
│       └── {chip}_bin100/                  # 100µm resolution
├── 07.analysis/
│   ├── {chip}_bin{size}_cluster.tsv        # Cluster assignments
│   ├── {chip}_bin{size}_umap.tsv           # UMAP coordinates
│   ├── {chip}_bin{size}_tsne.tsv           # t-SNE coordinates (if generated)
│   ├── {chip}_bin{size}_markers.csv        # Marker genes per cluster
│   └── spatial_plots/                      # Spatial visualization PNGs
├── {chip}_spatial_analysis_report.html     # ⭐ Main interactive report
└── pipeline.log                            # Complete pipeline execution log
```

### Key Output Files

| File | Description | Use Case |
|------|-------------|----------|
| **`{chip}_spatial_analysis_report.html`** | **Main deliverable** - Interactive HTML report with all QC metrics, plots, and visualizations | Share with collaborators, publication figures |
| `06.binSegment/square_bin/{chip}_bin50/filtered_feature_bc_matrix/` | Count matrix for downstream analysis | Load into Seurat, Scanpy, or other tools |
| `06.binSegment/{chip}_3_overlay_combined.jpg` | [HE mode] HE-GEM alignment visualization | Verify registration quality |
| `07.analysis/{chip}_bin{size}_cluster.tsv` | Spatial cluster assignments | Downstream spatial analysis |
| `pipeline.log` | Complete execution log with timing and errors | Troubleshooting, performance analysis |

### Report Features (HTML)

The interactive HTML report includes:

**QC Metrics Section:**
- Total reads, mapped reads, genes detected
- Sequencing saturation curve
- Barcode rank plot
- Median genes per square

**Spatial Visualization:**
- UMI distribution heatmap
- Gene expression overlay on tissue
- Cluster spatial distribution

**Clustering Analysis:**
- UMAP/t-SNE embeddings colored by cluster
- Top marker genes per cluster
- Cluster composition statistics

**[HE Mode Only] Image Alignment Viewer:**
- Interactive overlay of H&E and gene expression
- Toggle visibility of HE/GEM layers
- Opacity slider for blend control
- Download buttons for HE-only, GEM-only, or combined views

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
[https://github.com/](https://github.com/celatlas/Celatlas_spatial/)
```

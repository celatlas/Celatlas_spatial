# Celatlas Spatial v1.8 Linux Deployment Guide

This document describes how to install the `celatlas-spatial-v1.8.0-linux-x86_64`
deployment bundle on a new Linux server. The package contains the Celatlas
wheel, source archive, Conda environment definition, H&E StarDist model, and
the Swin tissue-segmentation model. Reference genomes remain separate assets.

## 1. Prerequisites

- Linux x86_64 server with sufficient CPU, memory, and disk space for STAR.
- A normal Linux account with read/write access to the chosen software and data
  directories.
- Miniforge, Mambaforge, or another Conda installation. Mamba is recommended.
- Network access to Conda/PyPI, or a configured internal mirror. An offline
  server needs a prepared Conda package cache or a prebuilt Conda environment.
- A STAR reference directory for every species/workflow to be run.

The deployment package does not contain FASTQ files, sample barcode files,
images, or reference indexes.

## 2. Transfer and Extract the Bundle

Copy the following archive to the new server:

```text
celatlas-spatial-v1.8.0-linux-x86_64.tar.gz
```

Extract it in the software directory:

```bash
mkdir -p /data/software
tar -xzf celatlas-spatial-v1.8.0-linux-x86_64.tar.gz -C /data/software
cd /data/software/celatlas-spatial-v1.8.0-linux-x86_64
sha256sum -c CHECKSUMS.sha256
```

All checksum lines must end in `OK` before installation continues.

## 3. Install Miniforge When Conda Is Unavailable

Skip this section if `conda --version` works. For an x86_64 server with network
access:

```bash
wget -O /tmp/Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash /tmp/Miniforge3.sh -b -p /data/miniforge3
source /data/miniforge3/etc/profile.d/conda.sh
conda install -n base -c conda-forge mamba -y
```

For a disconnected server, transfer an approved Miniforge installer and use
the organization's Conda mirror or offline package cache.

## 4. Create the Celatlas Environment

Make Conda available on `PATH`, then run the bundled installer:

```bash
source /data/miniforge3/etc/profile.d/conda.sh
cd /data/software/celatlas-spatial-v1.8.0-linux-x86_64
./release/install.sh
```

The installer creates `celatlas18`, installs the Celatlas wheel,
and explicitly installs `subread`/`featureCounts` and `samtools`. These two
packages are required even though they are not listed in the current YAML file.

Verify the finished environment:

```bash
./release/verify.sh
```

Expected tools include `STAR`, `featureCounts`, `samtools`, `cutadapt`,
`celatlas`, and `celatlas_spatial`.

## 5. Prepare Models and References

The installer copies the bundled H&E StarDist model into:

```text
~/.keras/models/StarDist2D/2D_versatile_he/
```

Copy the bundled Swin model to the directory that will be configured as
`paths.src_dir`:

```bash
mkdir -p /data/celatlas/src
cp models/swin_tiny.pth /data/celatlas/src/swin_tiny.pth
```

For H&E or fluorescence StarDist cell segmentation, the corresponding model
must be available to the Linux account that runs Celatlas. This bundle includes
H&E only; `2D_versatile_fluo` is needed for offline fluorescence/ssDNA cell
segmentation.

Copy STAR references from a validated server. ST and SN references use this
layout:

```text
/data/celatlas/reference/Homo_sapiens/
  celatlas_spatial_genome.config
  Genome
  SA
  SAindex
  ...
```

For SX, `sx_reference_dir` must point directly at the corresponding FFPE/panel
reference directory. Do not copy only FASTA and GTF files; the generated STAR
index and `celatlas_spatial_genome.config` are also required.

## 6. Configure Paths

Create a site-specific configuration file outside the extracted release so that
future software upgrades do not overwrite it. For example, create
`/data/celatlas/config/runner.yaml` from `config-templates/runner.yaml.example`:

```yaml
paths:
  workspace: /data/celatlas
  results_root: /data/celatlas/results
  fastq_root: /data/celatlas/fastq
  mask_dir: /data/celatlas/mask
  image_dir: /data/celatlas/images
  reference_dir: /data/celatlas/reference
  sx_reference_dir: /data/celatlas/reference/Homo_sapiens_wtpanel
  src_dir: /data/celatlas/src

defaults:
  default_thread: 32
  default_star_thread: 12
  default_bin: 50
  default_pixel_size: 0.5

sx:
  featurecounts_thread: 8

env:
  CELATLAS_CONDA_ROOT: /data/miniforge3
  CELATLAS_ENV_NAME: celatlas18
  CELATLAS_ENV_PREFIX: /data/miniforge3/envs/celatlas18
  CELATLAS_AUTO_ACTIVATE: "1"
  CELATLAS_USE_ENV_PATH: "1"
```

Create the data directories if needed:

```bash
mkdir -p /data/celatlas/{fastq,mask,images,results,reference,src,config}
```

## 7. Runtime Environment

Before every manual run, activate the Conda environment and isolate user Python
packages and cache files:

```bash
source /data/miniforge3/etc/profile.d/conda.sh
conda activate celatlas18
export PYTHONNOUSERSITE=1
export MPLCONFIGDIR=/tmp/celatlas_mplconfig
export NUMBA_CACHE_DIR=/tmp/celatlas_numba_cache
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
```

## 8. First Dry Run

Use one real sample and dry-run first. `--dry-run` checks FASTQ discovery,
barcode positions, image paths, and reference paths without starting analysis.

```bash
celatlas count \
  --id SX000293_A1 \
  --sample-name test \
  --tissue lung \
  --targetdir /data/celatlas/results/test/SX000293_A1 \
  --chemistry BBV2.4 \
  --species Homo_sapiens \
  --workflow SX \
  --image HE \
  --mode strna \
  --fastqs /data/celatlas/fastq/BBV2.4 \
  --config /data/celatlas/config/runner.yaml \
  --dry-run
```

For an HE sample, the preflight check requires:

```text
/data/celatlas/mask/SX000293_A1_FilterBarcodes.csv
/data/celatlas/images/SX000293_A1_he.tif
```

The H&E image extension may also be `tiff`, `png`, `jpg`, or `jpeg`. Remove
`--dry-run` only after the printed plan and all input paths are correct.

## 9. Upgrade and Reproducibility

- Preserve the release archive, `CHECKSUMS.sha256`, site configuration, and
  reference version information with every analysis project.
- Keep references and runtime data outside the software release directory.
- Install a new version in a separate directory and reuse the site-specific
  configuration only after reviewing new options and model requirements.
- Record `celatlas --help`, `python --version`, and the release SHA256 in the
  deployment record.

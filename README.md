# Celatlas Spatial v1.8.0

[![Version](https://img.shields.io/badge/version-1.8.0-1f883d.svg)](https://github.com/celatlas/Celatlas_spatial/releases/tag/v1.8.0)
[![License](https://img.shields.io/badge/license-MIT-2ea44f.svg)](LICENSE.txt)
[![Platform](https://img.shields.io/badge/platform-Linux%20x86__64-555.svg)](#installation-and-deployment)

Celatlas Spatial is a Linux pipeline for spatial transcriptomics and related RNA workflows. It turns FASTQ files, spatial coordinates, and optional HE or ssDNA images into reproducible count matrices, spatial bins, cell-level results, and HTML reports.

**v1.8.0 is a major standalone update.** It keeps the v1.7-compatible shell entry points while adding the unified `celatlas` CLI, ST/SX/SN workflows, preflight checks, resumable execution, reanalysis and report-only modes, StarDist cell segmentation, and cell-level analysis. The v1.8 source is on the `v1.8.0` branch, and the deployment bundle is published in [GitHub Releases](https://github.com/celatlas/Celatlas_spatial/releases/tag/v1.8.0).

## What's new in v1.8.0

- **Unified public CLI:** `celatlas count`, `celatlas reanalyze`, and `celatlas mkreport`, with `celatlas_spatial` retained as an advanced entry point.
- **Three workflows:** ST, SX, and SN library designs, with `strna` and `scrna` analysis modes.
- **Image-aware processing:** gene expression, ssDNA, and HE inputs with tissue segmentation, HE registration, spatial binning, and interactive reports.
- **Cell-level analysis:** optional StarDist HE/fluorescence segmentation with cell matrices, AnnData, spatial coordinates, QC, clustering, and marker results.
- **Reproducible execution:** FASTQ layout preflight, execution plans, step-level control, failure reports, resume support, and reanalysis from existing count outputs.
- **Deployment bundle:** Conda environment definitions, Linux install and verification scripts, bilingual manuals, and de-identified synthetic FASTQ smoke-test data.

## Quick start

### 1. Download the release bundle

Download the Linux bundle from the [v1.8.0 Release](https://github.com/celatlas/Celatlas_spatial/releases/tag/v1.8.0):

```text
celatlas-spatial-v1.8.0-linux-x86_64.tar.gz
```

```bash
tar -xzf celatlas-spatial-v1.8.0-linux-x86_64.tar.gz
cd celatlas-spatial-v1.8.0-linux-x86_64
./release/install.sh
conda activate celatlas18
./release/verify.sh
```

The installer sets up the software and Python/Conda dependencies. FASTQ files, barcode and mask inputs, reference indexes, and project data must be prepared separately as described in the deployment guide.

### 2. Plan and validate the inputs

```bash
celatlas count \
  --id SX000293_A1 \
  --sample-name Human_FFPE_001 \
  --tissue lung \
  --targetdir /data/celatlas/results/SXV1.1test/SX000293_A1 \
  --chemistry BBV2.4 \
  --species Homo_sapiens \
  --workflow SX \
  --image HE \
  --mode strna \
  --fastqs /data/celatlas/fastq/BBV2.4 \
  --reference-dir /data/celatlas/reference/Homo_sapiens_wtpanel \
  --dry-run
```

Remove `--dry-run` after the plan and input checks pass. Add `--enable-cell-segmentation` when cell segmentation is needed. See the deployment guide for image, model, and GPU requirements.

### 3. Reuse existing results

```bash
# Re-run downstream analysis from existing count outputs
celatlas reanalyze --targetdir /data/celatlas/results/SXV1.1test/SX000293_A1

# Regenerate only the HTML report
celatlas mkreport --targetdir /data/celatlas/results/SXV1.1test/SX000293_A1
```

## Entry points

| Entry point | Purpose |
|---|---|
| `celatlas count` | New sample: FASTQ to count, analysis, and report |
| `celatlas reanalyze` | Run downstream analysis from existing count outputs |
| `celatlas mkreport` | Generate or refresh an HTML report |
| `celatlas_spatial run` | Config-file/CSV batch processing and step-level control |
| `Celatlas.sh`, `Celatlas_FFPE.sh`, `Celatlas_SN.sh`, `Celatlas_reanalysis.sh` | Compatible shell entry points |

## Outputs

The standard spatial workflow produces QC metrics, tissue and image-registration outputs, square-bin matrices at selected resolutions, UMAP/clustering/marker results, and an interactive HTML report. With cell segmentation enabled, it also produces:

- `cell_matrix/`: 10X-compatible cell-level matrix
- `*_cell_adata.h5ad`: AnnData object
- `*_cell_positions.tsv` and `*_cell_metadata.tsv`: cell coordinates and QC metadata
- Cell boundaries and overlays, cell-level QC, clustering, spatial plots, and markers

Do not commit project data, reference indexes, or analysis results to the source repository.

## Installation and deployment

- [English user manual](docs/celatlas_spatial_manual_en.md)
- [Chinese user manual](docs/celatlas_spatial_manual_zh.md)
- [Linux deployment guide](docs/linux_deployment_guide.md)
- [Conda environment definition](envs/celatlas18.yml)
- [Pip dependency manifest](envs/celatlas18.requirements.txt)
- [Release notes](RELEASE_NOTES.md)

Minimum requirements are Linux x86_64, Python 3.11, and a working Conda or Miniforge installation. Complete workflows also require STAR, featureCounts, samtools, and cutadapt. HE and cell-segmentation workflows require the corresponding images, models, and computing resources.

## Release files

The source repository contains versioned code, documentation, and small examples. The Linux deployment bundle, models, and checksums are distributed as Release assets:

```text
celatlas-spatial-v1.8.0-linux-x86_64.tar.gz
```

After downloading, read `RELEASE_NOTES.md` and `release/MISSING_RUNTIME_ASSETS.md`, and keep the `CHECKSUMS.sha256` file for verification.

## License and feedback

This project is distributed under the [MIT License](LICENSE.txt). Please submit issues, feature requests, and installation feedback through [GitHub Issues](https://github.com/celatlas/Celatlas_spatial/issues).

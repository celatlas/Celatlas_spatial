# Celatlas Spatial

Spatial transcriptomics analysis from sequencing reads to expression matrices, spatial maps, and HTML reports.

**v1.8.0** is the current development and source release on this branch. It adds a unified command-line interface, ST/SX/SN workflow selection, input checks and resumable runs, reanalysis from existing count outputs, and optional StarDist cell segmentation. The [v1.7.0 branch](https://github.com/celatlas/Celatlas_spatial/tree/v1.7.0) remains available for users of the previous workflow.

> **Download status:** The [v1.8.0 GitHub Release](https://github.com/celatlas/Celatlas_spatial/releases/tag/v1.8.0) currently contains GitHub-generated source archives only. The Linux deployment bundle `celatlas-spatial-v1.8.0-linux-x86_64.tar.gz` has not yet been attached. The source archives do not include the bundled models or installer layout described below.

## What Celatlas Spatial does

The pipeline processes FASTQ files with the appropriate barcode layout, reference index, spatial coordinates, and optional tissue images. Depending on the run mode, it produces gene-count matrices, spatial bins, image-registration results, clustering and marker analyses, and an HTML report. Cell segmentation is optional and requires a compatible image and StarDist model.

| Mode | Typical input | Main result |
|---|---|---|
| `strna` | Spatial barcodes and, optionally, HE or ssDNA images | Spatial bin matrices and maps |
| `scrna` | Barcode-based RNA data | Cell expression matrix and analysis |
| `reanalyze` | Existing Celatlas count outputs | Updated downstream analysis without repeating alignment |
| `mkreport` | Existing Celatlas analysis outputs | Refreshed HTML report |

`ST`, `SX`, and `SN` select the library workflow. `SN` refers to the random-primer library design; choose it from the experimental protocol.

## What's new in v1.8.0

- A public `celatlas` command with `count`, `reanalyze`, and `mkreport` subcommands.
- A Python runner for explicit input preflight, dry-run plans, step-level execution, logs, and resume decisions.
- Optional StarDist segmentation for supported HE or fluorescence images, with cell-level matrices, AnnData outputs, QC, clustering, and marker analysis.
- Updated output organization and report generation, plus English and Chinese user manuals.
- Compatibility shell entry points for migration from v1.7.

## Getting started

The supported deployment path uses the Linux x86_64 bundle. Once it is attached to the [v1.8.0 Release](https://github.com/celatlas/Celatlas_spatial/releases/tag/v1.8.0), download it and run:

```bash
tar -xzf celatlas-spatial-v1.8.0-linux-x86_64.tar.gz
cd celatlas-spatial-v1.8.0-linux-x86_64
sha256sum -c CHECKSUMS.sha256
./release/install.sh
conda activate celatlas18
./release/verify.sh
```

Conda or Miniforge must be available before installation. The bundle provides the software, environment definition, documentation, and selected models. You must supply project FASTQs, barcode/mask files, images, and validated STAR reference indexes separately. See the [Linux deployment guide](docs/linux_deployment_guide.md) for server setup and model paths.

### Plan a new run

This example prints the planned steps and checks input paths. Replace the sample IDs and `/data/...` paths with your own values:

```bash
celatlas count \
  --id SX000293_A1 \
  --sample-name Human_FFPE_001 \
  --tissue lung \
  --targetdir /data/celatlas/results/project_001/SX000293_A1 \
  --chemistry BBV2.4 \
  --species Homo_sapiens \
  --workflow SX \
  --image HE \
  --mode strna \
  --fastqs /data/celatlas/fastq/BBV2.4 \
  --reference-dir /data/celatlas/reference/Homo_sapiens_wtpanel \
  --dry-run
```

When the plan and required inputs are correct, remove `--dry-run` to run the pipeline. Add `--enable-cell-segmentation` only when the appropriate image and model are available. For a shorter installation smoke test, use the synthetic FASTQ data included in the Linux bundle as described in the [user manual](docs/celatlas_spatial_manual_en.md).

### Reanalyze or regenerate a report

`reanalyze` and `mkreport` both require the sample metadata arguments, even when the output directory already exists. Use the same values as the original run and change only the subcommand:

```bash
celatlas reanalyze \
  --id SX000293_A1 \
  --sample-name Human_FFPE_001 \
  --tissue lung \
  --targetdir /data/celatlas/results/project_001/SX000293_A1 \
  --chemistry BBV2.4 \
  --species Homo_sapiens \
  --workflow SX \
  --image HE \
  --mode strna
```

Use `celatlas mkreport` with those same arguments to regenerate only the report. Run `celatlas <subcommand> --help` for the full option list.

## Documentation

- [English user manual](docs/celatlas_spatial_manual_en.md)
- [Chinese user manual](docs/celatlas_spatial_manual_zh.md)
- [Linux deployment guide](docs/linux_deployment_guide.md)
- [Conda environment definition](envs/celatlas18.yml)
- [Release notes](RELEASE_NOTES.md)

The source tree can also be built as a Python wheel and source archive with `python -m build`. Building from source still requires external command-line tools, models, reference indexes, and runtime data for a complete analysis.

## License and support

Celatlas Spatial is released under the [MIT License](LICENSE.txt). Please report installation problems and bugs through [GitHub Issues](https://github.com/celatlas/Celatlas_spatial/issues).

# Celatlas Spatial v1.8.0 User Manual

This is the public English manual for Celatlas Spatial v1.8.0. It is written for Linux command-line users, analysts, and deployment engineers. The interface described here is the interface shipped in the public release bundle. Paths, sample identifiers, and commands containing `/data/...` are illustrative and must be replaced with local values.

> Release: v1.8.0  
> Platform: Linux x86_64  
> Recommended Python: 3.11

## 1. What the software does

Celatlas Spatial processes spatial transcriptomics or single-cell/single-nucleus RNA data from FASTQ files to expression matrices, spatial bins, optional cell segmentation, downstream clustering, and HTML reports. The public architecture has three layers:

1. **Public CLI**: `celatlas count`, `celatlas reanalyze`, and `celatlas mkreport`.
2. **Python runner**: configuration merging, path resolution, preflight checks, step planning, resume decisions, and logs.
3. **Analysis steps**: barcode extraction, cutadapt, STAR, featureCounts, UMI counting, binSegment, StarDist, analysis, and report generation.

The normal plan is:

```text
00.sample
  -> 01.barcode
  -> 02.cutadapt
  -> 03.star
  -> 04.featureCounts
  -> 05.count
  -> 06.segment/01.binsegment
  -> 06.segment/02.cellsegment       (optional)
  -> 07.outs/binned_outputs
  -> 07.outs/cellsegmented_outputs   (optional)
  -> 08.report
```

`strna` keeps spatial barcode coordinates and creates spatial bins. `scrna` treats barcodes as independent cells and goes directly to cell expression analysis. `ST`, `SX`, and `SN` are workflow types; `denovo`, `reanalysis`, and `report` describe the starting point of a run.

## 2. Changes from v1.7

- The recommended interface is the `celatlas` command with explicit metadata and input validation.
- `reanalyze` reruns downstream work without repeating FASTQ alignment and counting. `mkreport` refreshes only the final report.
- FASTQ manifests record path, size, and modification time for every R1/R2 file. This prevents stale count results from being silently reused after a sequencing change.
- `--dry-run`, `--preflight`, `--show-steps`, and `--show-resume-plan` make hidden execution decisions inspectable.
- v1.8 adds optional cavity filtering, StarDist cell segmentation, the DAPI preset, and resumable tiled inference for very large images.
- Output directories use `06.segment`, `07.outs`, and `08.report`. Legacy shell entry points remain for migration.
- The public bundle excludes private production tables, web schedulers, logs, reference indexes, and experiment data.

## 3. Installation and verification

### 3.1 Install the release bundle

For example, extract the deployment bundle to:

```text
/data/software/celatlas-spatial-v1.8.0-linux-x86_64/
```

After extraction, verify every file:

```bash
cd /data/software/celatlas-spatial-v1.8.0-linux-x86_64
sha256sum -c CHECKSUMS.sha256
```

Install with the bundled Conda helpers:

```bash
source /data/miniforge3/etc/profile.d/conda.sh
./release/install.sh
conda activate celatlas18
./release/verify.sh
```

Verification should find `celatlas`, `celatlas_spatial`, `STAR`, `featureCounts`, `samtools`, and `cutadapt`. If Conda is not available, install Miniforge/Mambaforge or prepare the environment offline as described in `docs/linux_deployment_guide.md`.

### 3.2 Recommended runtime variables

```bash
source /data/miniforge3/etc/profile.d/conda.sh
conda activate celatlas18
export PYTHONNOUSERSITE=1
export MPLCONFIGDIR=/tmp/celatlas_mplconfig
export NUMBA_CACHE_DIR=/tmp/celatlas_numba_cache
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
celatlas --help
celatlas_spatial_runner --help
```

Large image jobs need enough temporary disk space. Put temporary and tile-cache directories on a writable filesystem with sufficient capacity.

## 4. Release contents and external assets

## 4.1 FASTQ smoke-test demo

The release bundle includes `demo_data/fastq/BBV2.4/` with 20,000 paired read pairs and `release/create_fastq_demo.py`. The subset is intended to verify installation, FASTQ discovery, and dry-run planning; it is not a biological benchmark and cannot replace a complete run with a reference and matching spatial assets. Every FASTQ header is synthetic, bases and quality values are synthetic, and the manifest does not record the source path.

Run the smoke test from the unpacked bundle:

```bash
cd /data/software/celatlas-spatial-v1.8.0-linux-x86_64
source /data/miniforge3/etc/profile.d/conda.sh
conda activate celatlas18
celatlas count \
  --id CELATLAS_DEMO \
  --sample-name Celatlas_Demo \
  --tissue demo \
  --targetdir "$PWD/demo_results/CELATLAS_DEMO" \
  --chemistry BBV2.4 \
  --species Mus_musculus \
  --workflow ST \
  --image gene_expr \
  --mode scrna \
  --fastqs "$PWD/demo_data/fastq/BBV2.4" \
  --reference-dir /data/celatlas/reference \
  --thread 2 \
  --bin 50 \
  --skip-preflight \
  --dry-run
```

For an actual count run, provide a complete STAR reference, chemistry-compatible inputs, and the barcode/mask files required by `strna`. The 20,000-read subset is normally too small for meaningful mapping or QC.

To create a new de-identified subset from an authorized paired FASTQ input:

```bash
python release/create_fastq_demo.py \
  --r1 /data/input/sample_R1.fq.gz \
  --r2 /data/input/sample_R2.fq.gz \
  --outdir demo_data/fastq/custom \
  --sample CELATLAS_DEMO \
  --reads 20000
```

The tool preserves paired record lengths, validates R1/R2 synchronization, rewrites headers as `CELATLAS_DEMO:NNNNNNNN/{1,2}`, replaces bases and qualities by default, and writes a SHA256 manifest without source paths or original instrument metadata. Use `--sequence-mode copy` only when the biological reads are authorized for redistribution.



```text
celatlas-spatial-v1.8.0-linux-x86_64/
├── dist/                 wheel and source archive
├── docs/                 Chinese and English manuals, deployment notes
├── envs/                 Conda and pip definitions
├── config-templates/     portable configuration examples
├── models/               Swin model and optional StarDist H&E model
├── release/              install.sh, verify.sh, runtime asset notes
├── CHECKSUMS.sha256
└── RELEASE_NOTES.md
```

Prepare these assets separately:

- FASTQ files;
- spatial barcode files or masks;
- H&E, ssDNA, or DAPI images;
- an ST/SN species STAR reference or an SX panel reference;
- `swin_tiny.pth` copied to the configured `src_dir`;
- StarDist models when cell segmentation is enabled. The bundle commonly includes the H&E model; `2D_versatile_fluo` must be prepared for fluorescence/DAPI runs.

A reference is more than FASTA and GTF. It needs the generated STAR index and the Celatlas configuration file, such as `Genome`, `SA`, `SAindex`, and `celatlas_spatial_genome.config`.

### 4.1 Build a reference from FASTA/GTF

For a GFF3 annotation, first convert it to a GTF with `gene_id` attributes:

```bash
python scripts/normalize_gff3_to_gtf.py \
  /data/reference/raw/genes.gff3 \
  /data/reference/raw/genes.celatlas.gtf
```

Build the RNA reference with:

```bash
celatlas_spatial rna mkref \
  --genome_name Mus_musculus_demo \
  --fasta /data/reference/raw/Mus_musculus.fa \
  --gtf /data/reference/raw/genes.celatlas.gtf \
  --outdir /data/celatlas/reference/Mus_musculus \
  --thread 16
```

Add `--mt_gene_list /data/reference/raw/mt_genes.txt` when a mitochondrial gene list is available. Inspect the generated configuration without building the index:

```bash
celatlas_spatial rna mkref ... --dry_run
```

Record the genome name, FASTA, GTF, and output directory. An SX panel reference should be built directly at the configured `sx_reference_dir`; do not use a generic species reference as a panel reference.

## 5. The three choices you make for every run

### 5.1 Workflow: ST, SX, or SN

| Workflow | Meaning | Typical use | Reference resolution |
|---|---|---|---|
| `ST` | regular spatial transcriptomics | fresh/frozen or standard libraries | `reference/<species>` |
| `SX` | FFPE or targeted panel | probe, FFPE, targeted-panel libraries | `sx_reference_dir` points directly to the panel reference |
| `SN` | Random-primer workflow | Spatial or nucleus-oriented libraries using the random-primer design | `reference/<species>` |

`FF` and `FFPE` remain accepted aliases by compatibility runners, while the public `celatlas` CLI accepts the canonical `ST`, `SX`, and `SN` values. `SN` identifies the random-primer library design; do not select it from the words “single nucleus” alone. Confirm the library chemistry and protocol before choosing `SN`.

### 5.2 Image method: gene_expr, ssDNA, or HE

| Method | Input | Tissue detection | Use when |
|---|---|---|---|
| `gene_expr` | no tissue image | inferred from GEM/UMI distribution | no image is available |
| `ssDNA` | fluorescence `<id>.tif` | traditional fluorescence mask | ssDNA or fluorescence imaging was used |
| `HE` | `<id>_he.tif` or supported image extension | expression mask plus H&E registration | pathology and morphology work |

H&E files must use the `_he` suffix. `<id>.tif` is the default ssDNA naming pattern and should not be used for an H&E image.

### 5.3 Analysis mode: strna or scrna

| Mode | Spatial coordinates | Main outputs |
|---|---|---|
| `strna` | used | bin10/20/50/100, spatial clusters, spatial plots |
| `scrna` | not used | barcode/cell matrix, PCA/UMAP/markers |

Choose `strna` when a spatial barcode coordinate file is available and tissue location matters. Choose `scrna` for cell or nucleus expression analysis without spatial binning.

## 6. Input preparation

### 6.1 Recommended layout

```text
/data/celatlas/
├── fastq/BBV2.4/
│   ├── SX000293_A1_S1_L001_R1_001.fastq.gz
│   ├── SX000293_A1_S1_L001_R2_001.fastq.gz
│   └── ...
├── mask/
│   ├── SX000293_A1.barcodeToPos.h5
│   ├── SX000293_A1_FilterBarcodes.csv
│   └── SX000293_A1_tissue_bbox.csv
├── images/
│   └── SX000293_A1_he.tif
├── reference/
│   └── Homo_sapiens/
├── src/
│   └── swin_tiny.pth
└── results/
```

The public spatial input is `<id>_FilterBarcodes.csv`. The barcode step prefers `<id>.barcodeToPos.h5` and falls back to the CSV. If present, `<id>_tissue_bbox.csv` contributes the tissue canvas extent.

### 6.2 FASTQ discovery

`--fastq-name` defaults to `--id`. The runner searches these layouts in order and does not mix layouts after one is found:

1. Multi-lane: `<prefix>_S*_L*_R1_*.fastq.gz` and matching R2; `.fq.gz` is also accepted.
2. Fold: `<prefix>_fold1_1.fq.gz` and `<prefix>_fold1_2.fq.gz`, through fold5.
3. Simple: `<prefix>_1.fq.gz` and `<prefix>_2.fq.gz`.

R1 and R2 counts must match. Every file must exist, be non-empty, and be readable. When a sequencing provider uses a different prefix, set `--fastq-name`:

```bash
--id SX000293_A1 --fastq-name Sequencer_Sample_20261001
```

### 6.3 Image names

- `HE`: `<id>_he.tif`, `.tiff`, `.png`, `.jpg`, or `.jpeg`.
- `ssDNA`: `<id>.tif` or the equivalent staged fluorescence image.
- `gene_expr`: no image is required.
- DAPI preset: provide the manual H&E ROI mask. After registration, the preset uses the blue channel of the full-resolution registered TIFF; it does not use the hires preview.

### 6.4 Coordinate consistency

Barcode coordinates, tissue bbox, images, registered images, and StarDist labels must use the same registered coordinate system. An existing file passed through `--stardist-labels` must already align with `FilterBarcodes.csv`; Celatlas does not infer a transform from the filename.

### 6.5 Runtime staging and manual ROI masks

Public `count` and `reanalyze` stage spatial inputs under:

```text
<targetdir>/06.segment/mask/
```

Accepted manual mask names include `manual_mask.png`, `<id>_manual_mask.png`, `<id>_mask_manual.png`, `mask_manual.png`, `he_manual_mask.png`, `<id>_he_manual_mask.png`, `<id>_he_mask_manual.png`, and `he_mask_manual.png`. Place a mask under `--mask-dir` for automatic staging, or directly under `<targetdir>/06.segment/mask`; existing runtime files are preserved during reanalysis. The DAPI preset requires a manual H&E ROI for registration and requires the registered full-resolution TIFF and mask to share coordinates.

## 7. First validation run

Inspect the installed interfaces:

```bash
celatlas --help
celatlas count --help
celatlas reanalyze --help
celatlas mkreport --help
celatlas_spatial_runner run --help
celatlas_spatial_runner batch --help
```

Run a real sample as a dry run first:

```bash
celatlas count \
  --id SX000293_A1 \
  --sample-name Demo_Human_Lung \
  --tissue lung \
  --targetdir /data/celatlas/results/demo/SX000293_A1 \
  --chemistry BBV2.4 \
  --species Homo_sapiens \
  --workflow SX \
  --image HE \
  --mode strna \
  --fastqs /data/celatlas/fastq/BBV2.4 \
  --mask-dir /data/celatlas/mask \
  --image-dir /data/celatlas/images \
  --reference-dir /data/celatlas/reference \
  --config /data/celatlas/config/runner.yaml \
  --dry-run
```

`--targetdir` must be absolute. A dry run prints the plan and checks inputs without starting the full analysis.

## 8. Public commands

### 8.1 `celatlas count`

Runs the complete FASTQ-to-report workflow. Required values are `--id`, `--sample-name`, `--tissue`, `--targetdir`, `--chemistry`, `--species`, `--workflow`, `--image`, and `--fastqs`.

```bash
celatlas count \
  --id DEMO_ST000001_A1 --sample-name Demo_Mouse_Brain --tissue brain \
  --targetdir /data/celatlas/results/demo/DEMO_ST000001_A1 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image gene_expr --mode strna \
  --fastqs /data/celatlas/fastq/BBV2.4 \
  --mask-dir /data/celatlas/mask --reference-dir /data/celatlas/reference \
  --thread 32 --bin 50
```

### 8.2 `celatlas reanalyze`

Starts from an existing `05.count` and reruns binSegment, optional post-binSegment modules, `07.outs`, and the report. It requires `<targetdir>/05.count/<id>_count_detail.txt` and does not repeat FASTQ alignment.

```bash
celatlas reanalyze \
  --id DEMO_ST000001_A1 --sample-name Demo_Mouse_Brain --tissue brain \
  --targetdir /data/celatlas/results/demo/DEMO_ST000001_A1 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image HE --mode strna --mask-dir /data/celatlas/mask \
  --image-dir /data/celatlas/images --reference-dir /data/celatlas/reference \
  --bin 50
```

`--skip-binsegment` keeps the existing binSegment output and starts later. `--skip-analysis` keeps `07.outs` and refreshes the report only. Enabling cavity filtering or cell segmentation makes reanalysis plan the required step range automatically.

### 8.3 `celatlas mkreport`

Generates only the final report from existing outputs:

```bash
celatlas mkreport \
  --id DEMO_ST000001_A1 --sample-name Demo_Mouse_Brain --tissue brain \
  --targetdir /data/celatlas/results/demo/DEMO_ST000001_A1 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image gene_expr --mode strna
```

## 9. Main options

### 9.1 Metadata and paths

| Option | Meaning |
|---|---|
| `--id` | chip/run ID, output stem, and default FASTQ prefix |
| `--sample-name` | biological sample name stored in metadata and reports |
| `--tissue` | tissue name stored in the report |
| `--targetdir` | final sample directory; must be absolute |
| `--project` | project/case name; defaults to the target directory parent |
| `--fastqs` | FASTQ directory; required by `count` |
| `--fastq-name` | actual FASTQ file prefix; defaults to `--id` |
| `--mask-dir` | barcode and tissue-coordinate directory |
| `--image-dir` | source H&E or fluorescence image directory |
| `--reference-dir` | ST/SN reference root or SX panel reference |
| `--config` | YAML, JSON, or simple env configuration |
| `--profile` | profile name or profile file |
| `--log-dir` | runner log and summary directory |

### 9.2 Threads, bins, and filters

| Option | Function | Guidance |
|---|---|---|
| `--thread` | general CPU thread count | limit on shared servers |
| `--star-thread` | STAR thread count | often lower than general threads |
| `--featurecounts-thread` | featureCounts threads | 8 is common for SX |
| `--bin` | spatial bin size, usually 10/20/50/100 | 50 is a balanced default |
| `--pixel-size` | microns per image pixel | match the imaging system |
| `--cluster-resolution` | Leiden resolution for spatial bins | larger values generally create more clusters |
| `--cell-cluster-resolution` | Leiden resolution for cells | used with cell analysis |
| `--gem-bin-size` | GEM heatmap aggregation | 10 is finer, 20 uses fewer resources |
| `--cell-num` | expected cells/spots for count calling | tune to tissue size |
| `--cell-min-genes` / `--cell-min-counts` | lower cell filters | inspect QC before changing |
| `--cell-max-genes` / `--cell-max-counts` / `--cell-max-mt` | upper and mitochondrial filters | record every change |

bin10 is close to single-cell resolution but is large. bin20 is useful for fine structure; bin50 is the everyday balance; bin100 is useful for global patterns. A bin is a grid edge length and is not a cell diameter.

### 9.3 ssDNA options

| Option | Function |
|---|---|
| `--ssdna-threshold-scale` | scales the Otsu threshold; below 1 retains dim edges |
| `--ssdna-mask-expand-pixels` | expands the final tissue mask in registered-image pixels |
| `--ssdna-min-hole-area` | fills enclosed mask holes below this area |
| `--fluorescence-background` | keeps fluorescence tissue backgrounds black during H&E ROI registration |

These values are image-pixel or threshold-scale parameters, not microns. Save the mask QC when changing them.

### 9.4 Cavity filter

`--cavity-mode off|qc|apply` controls tissue-mask filtering:

- `off`: do not run the filter.
- `qc`: create filtered outputs and QC without replacing the original `square_bin`.
- `apply`: back up the original and replace the bin output used downstream.

Related options are `--cavity-min-component-ratio`, `--cavity-min-hole-area`, `--cavity-close-radius`, `--preserve-gem-support`, and `--cavity-force`. Run `qc` first. `apply` creates a backup, but the sample directory should still be copied or versioned before replacement.

### 9.5 StarDist cell segmentation

Enable ordinary segmentation with:

```bash
--enable-cell-segmentation
```

For H&E plus a manual ROI and DAPI/blue fluorescence, use:

```bash
--cell-segmentation-preset dapi
```

The DAPI preset is valid only with `--image HE --mode strna`. It requires the manual H&E ROI, enables a black background, performs full-resolution tiled fluorescence segmentation on the registered TIFF, and runs cell matrix and report steps. The default channel is `blue`:

```bash
--stardist-fluorescence-channel blue
```

The following advanced controls are available through compatibility shell entry points, configuration/profiles, or the low-level `celatlas_spatial rna stardistCellSegment` command. The public `celatlas` command directly exposes `--enable-cell-segmentation`, `--cell-segmentation-preset dapi`, and `--stardist-fluorescence-channel`; the DAPI preset supplies the full-resolution tiled defaults.

Useful advanced controls include:

| Option | Function |
|---|---|
| `--stardist-model` | model name; HE defaults to `2D_versatile_he`, fluorescence/DAPI to `2D_versatile_fluo` |
| `--stardist-model-dir` | local model directory |
| `--stardist-labels` | an existing coordinate-aligned label TIFF; skips inference |
| `--stardist-image` | explicit segmentation image |
| `--stardist-tissue-bbox` | label canvas bbox |
| `--stardist-prob-thresh` | probability threshold, approximately 0.30 by default |
| `--stardist-nms-thresh` | NMS threshold |
| `--stardist-max-dim` | ordinary inference long-edge limit; DAPI uses full resolution by default |
| `--stardist-scale` | internal StarDist scale |
| `--stardist-n-tiles` | ordinary inference tile grid, such as `4,4` |
| `--stardist-expand-pixels` | label expansion before barcode assignment |
| `--stardist-min-umi` / `--stardist-min-genes` | minimum expression for retained cells |
| `--stardist-skip-counts` | labels/assignments only |
| `--stardist-no-label-output` | skip slide-sized label TIFF output |
| `--stardist-tiled-inference` | resumable overlapping tiles for large images |
| `--stardist-tile-size` | tile edge, default 4096 |
| `--stardist-tile-overlap` | tile context halo, default 256 |
| `--stardist-tile-merge-overlap` | label merge overlap threshold |
| `--stardist-tile-cache-dir` | tile cache directory |
| `--stardist-exclude-rectangle` | excluded `x0,y0,x1,y1` region; semicolon separates multiple regions |

Compatibility shell example for a large registered image:

```bash
bash Celatlas_reanalysis.sh DEMO_HE0002_A1 demo_case BBV2.4 Homo_sapiens HE strna \
  --sample Demo_HE --sampledir /data/celatlas/results/demo/DEMO_HE0002_A1 \
  --image_dir /data/celatlas/images --mask_dir /data/celatlas/mask \
  --reference_dir /data/celatlas/reference \
  --enable-stardist-cell-segment --stardist-tiled-inference \
  --stardist-tile-size 4096 --stardist-tile-overlap 256
```

### 9.6 Advanced workflow parameters

These controls are mainly used through `celatlas_spatial_runner run`, CSV/profile configuration, or compatibility shell entry points:

| Option/config key | Function |
|---|---|
| `--insert-r2` / `--insertR2` | R2 insert length; match the sequencing design |
| `--star-match-min` | STAR `outFilterMatchNmin` |
| `--star-match-ratio` | STAR minimum match-length ratio |
| `--star-score-ratio` | STAR minimum score ratio |
| `--star-multimap` | STAR `outFilterMultimapNmax` |
| `--bbv4-strna-barcode-mismatch` | allowed spatial-whitelist mismatches for BBV4/BBV4_L9 |
| `--resolve-multigene-umi` | keep the highest-read-support gene for multi-gene UMI groups and discard ties |
| `--fastq-root` | FASTQ root used to resolve the chemistry directory |
| `--workspace` / `--results-root` | default workspace and result roots |
| `--genome-dir` / `--src-dir` | explicit reference and model directories |
| `--stage-inputs` | stage spatial inputs and report the actions |

SN plant/single-nucleus tuning normally belongs in a profile or YAML `sn` section, including `sn_barcode_lownum`, `sn_cutadapt_min_length`, `sn_star_multimap`, `sn_featurecounts_param`, and `sn_resolve_multigene_umi`. The `configs/profiles/bbv4_hd_loose.env` profile is for BBV4/HD mapping-sensitivity diagnostics; keep a control run and QC before changing STAR filters.

## 10. Python runner

### 10.1 Step-level planning

`celatlas_spatial_runner run` accepts `shell` and `python` engines. `--engine python` exposes Python step plans and is useful for dry runs and explicit step selection.

```bash
celatlas_spatial_runner run \
  --engine python --pipeline reanalysis --workflow ST \
  --chip-number DEMO_ST000006_A1 --casno demo_project \
  --sample-name Demo --tissue brain --chemistry BBV2.4 \
  --species Mus_musculus --method gene_expr --mode strna \
  --sampledir /data/celatlas/results/demo/DEMO_ST000006_A1 \
  --from-step 06 --to-step 08.report --show-steps --dry-run
```

Stable selectors include `00.sample`, `01.barcode`, `02.cutadapt`, `03.star`, `04.featureCounts`, `05.count`, `06`, `06.segment/01.binsegment`, `06.segment/02.cellsegment`, `07`, and `08.report`.

Useful inspection flags are:

```text
--show-steps --show-step-status --show-resume-plan
--show-execution-preflight --preflight --dry-run
```

Use `--allow-existing-step-outputs` only after confirming that the existing outputs match the current inputs.

### 10.2 CSV batch interface

The generic batch runner requires at least:

```text
workflow,chip_number,casno,chemistry,species,method,mode
```

Common optional fields are `run_id`, `enabled`, `sample_name`, `tissue`, `thread`, `bin`, `fastq_dir`, `fastq_name`, `mask_dir`, `image_dir`, `reference_dir`, `resume_existing`, `gene_mask_filter`, `enable_cavity_filter`, `cavity_apply`, `enable_stardist_cell_segment`, and `extra_args`.

Example `samples.csv`:

```csv
run_id,enabled,workflow,sample_name,tissue,chip_number,casno,chemistry,species,method,mode,thread,bin,fastq_dir,mask_dir,image_dir,reference_dir
demo_st,1,ST,Mouse brain,brain,DEMO_ST0007_A1,demo_case,BBV2.4,Mus_musculus,gene_expr,strna,32,50,/data/celatlas/fastq/BBV2.4,/data/celatlas/mask,/data/celatlas/images,/data/celatlas/reference
demo_sx,1,SX,Human lung,lung,DEMO_SX0008_A1,demo_case,BBV2.4,Homo_sapiens,HE,strna,32,50,/data/celatlas/fastq/BBV2.4,/data/celatlas/mask,/data/celatlas/images,/data/celatlas/reference
```

Dry-run the table:

```bash
celatlas_spatial_runner batch \
  --samples /data/celatlas/config/samples.csv \
  --config /data/celatlas/config/runner.yaml \
  --preflight --dry-run --show-steps
```

Use `--only demo_sx` to select a run ID and `--stop-on-error` to stop after the first failure.

### 10.3 Low-level cavity and StarDist commands

For a standalone post-binSegment rerun, call the RNA subcommands directly:

```bash
celatlas_spatial rna cavityFilter \
  --sampledir /data/celatlas/results/demo/DEMO_HE0002_A1 \
  --bins 20,50 --mask-preset default --dry-run
```

With an existing label TIFF in the same registered coordinate system, skip inference:

```bash
celatlas_spatial rna stardistCellSegment \
  --sampledir /data/celatlas/results/demo/DEMO_HE0002_A1 \
  --labels /data/celatlas/labels/DEMO_HE0002_A1.labels.tif \
  --skip-inference --min-umi 10 --min-genes 5
```

For full-resolution large-image inference:

```bash
celatlas_spatial rna stardistCellSegment \
  --sampledir /data/celatlas/results/demo/DEMO_HE0002_A1 \
  --tiled-inference --tile-size 4096 --tile-overlap 256 \
  --fluorescence-channel blue --tile-cache-dir /data/celatlas/tile_cache
```

These low-level commands still require `FilterBarcodes.csv`, `tissue_bbox.csv`, the registered TIFF, and the label canvas to share coordinates.

## 11. Configuration and profiles

### 11.1 YAML

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
  gem_bin_size: 20

sx:
  featurecounts_thread: 8

env:
  CELATLAS_ENV_NAME: celatlas18
  CELATLAS_AUTO_ACTIVATE: "1"
  CELATLAS_USE_ENV_PATH: "1"
```

### 11.2 Environment file

```bash
export CELATLAS_WORKSPACE=/data/celatlas
export CELATLAS_REFERENCE_DIR=/data/celatlas/reference
export CELATLAS_MASK_DIR=/data/celatlas/mask
export CELATLAS_IMAGE_DIR=/data/celatlas/images
export CELATLAS_FASTQ_ROOT=/data/celatlas/fastq
export CELATLAS_RESULTS_ROOT=/data/celatlas/results
export CELATLAS_SRC_DIR=/data/celatlas/src
```

Pass it with `--config /path/to/site.env`. The parser reads `KEY=value` and `export KEY=value` records and does not execute shell control flow.

Configuration is merged from site defaults, profiles, workflow/method sections, CSV row values, and explicit command-line values. Later, more specific values override earlier ones. Keep credentials, patient data, and machine-specific paths outside the public repository.

## 12. Resume and rerun behavior

### 12.1 Safe resume

Start with:

```bash
celatlas count ... --dry-run
celatlas count ... --resume-existing
```

After a successful run, the baseline is stored at:

```text
<targetdir>/.fastq_inputs.tsv
```

The input set attempted by an in-progress run is stored as `.fastq_inputs.current.tsv`. It is promoted only after success, so a failed attempt does not become the new baseline.

### 12.2 FASTQ changes

`celatlas count` automatically creates and checks the FASTQ manifest. When files are added, reordered, or changed in size or modification time, sequencing-dependent outputs are not reused silently. Request an explicit reset:

```bash
celatlas count ... --reset-fastq-outputs
```

To keep the old output tree:

```bash
celatlas count ... --reset-fastq-outputs \
  --archive-fastq-outputs
```

Old outputs are moved under `.rerun_archive`. Avoid `--allow-existing-step-outputs` until the manifest and input compatibility are confirmed.

### 12.3 Downstream reruns

Use `reanalyze` after changing images, masks, bins, cavity settings, or cell segmentation. Use `mkreport` when only the report needs to be regenerated.

## 13. Output tree and result formats

```text
<targetdir>/
├── 00.sample/                         sample metadata and statistics
├── 01.barcode/                        barcode extraction
├── 02.cutadapt/                       adapter trimming and QC
├── 03.star/                           STAR BAM and logs
├── 04.featureCounts/                  gene-count BAM/tables
├── 05.count/                          UMI count and raw/filtered matrices
│   ├── <id>_count_detail.txt
│   └── <id>_filtered_feature_bc_matrix/
├── 06.segment/
│   ├── 01.binsegment/
│   │   ├── images/                    masks, registered image, previews
│   │   └── square_bin/                bin10/20/50/100/Raw
│   └── 02.cellsegment/                optional labels and assignments
├── 07.outs/
│   ├── binned_outputs/                spatial analysis and images
│   └── cellsegmented_outputs/         optional cell analysis
├── 08.report/                         report-step logs/intermediates
├── <id>_spatial_analysis_report.html
├── .fastq_inputs.tsv                  successful FASTQ baseline
└── pipeline.log / runner logs
```

Matrices normally use a 10X-style sparse layout:

```text
filtered_feature_bc_matrix/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

They can be read with Seurat `Read10X()` or Scanpy `sc.read_10x_mtx()`. Interpret UMI, gene, mapping, cell, cluster, and image metrics together with experimental design. Tissue-outside barcodes, low-quality calls, and registration errors directly affect downstream results.

## 14. Eight virtual demos

### Demo 1: gene expression only

```bash
celatlas count \
  --id DEMO_ST000001_A1 --sample-name Demo_Mouse_Brain --tissue brain \
  --targetdir /data/celatlas/results/demo/DEMO_ST000001_A1 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image gene_expr --mode strna \
  --fastqs /data/celatlas/fastq/BBV2.4 \
  --mask-dir /data/celatlas/mask --reference-dir /data/celatlas/reference \
  --thread 32 --bin 50 --dry-run
```

### Demo 2: SX + H&E + StarDist

```bash
celatlas count \
  --id DEMO_SX000002_A1 --sample-name Demo_Human_Lung --tissue lung \
  --targetdir /data/celatlas/results/demo/DEMO_SX000002_A1 \
  --chemistry BBV2.4 --species Homo_sapiens --workflow SX \
  --image HE --mode strna --fastqs /data/celatlas/fastq/BBV2.4 \
  --mask-dir /data/celatlas/mask --image-dir /data/celatlas/images \
  --reference-dir /data/celatlas/reference \
  --enable-cell-segmentation --dry-run
```

### Demo 3: DAPI preset

```bash
celatlas reanalyze \
  --id DEMO_ST000003_A1 --sample-name Demo_DAPI_Liver --tissue liver \
  --targetdir /data/celatlas/results/demo/DEMO_ST000003_A1 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image HE --mode strna --mask-dir /data/celatlas/mask \
  --image-dir /data/celatlas/images --reference-dir /data/celatlas/reference \
  --cell-segmentation-preset dapi --stardist-fluorescence-channel blue \
  --dry-run
```

### Demo 4: ssDNA cavity QC

```bash
celatlas reanalyze \
  --id DEMO_ST000004_A1 --sample-name Demo_Mouse_Kidney --tissue kidney \
  --targetdir /data/celatlas/results/demo/DEMO_ST000004_A1 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image ssDNA --mode strna --mask-dir /data/celatlas/mask \
  --image-dir /data/celatlas/images --reference-dir /data/celatlas/reference \
  --cavity-mode qc --ssdna-threshold-scale 0.85 \
  --ssdna-mask-expand-pixels 16 --dry-run
```

### Demo 5: report only

```bash
celatlas mkreport \
  --id DEMO_SX000005_A1 --sample-name Demo_Report --tissue colon \
  --targetdir /data/celatlas/results/demo/DEMO_SX000005_A1 \
  --chemistry BBV2.4 --species Homo_sapiens --workflow SX \
  --image HE --mode strna
```

### Demo 6: Python step control

```bash
celatlas_spatial_runner run --engine python --pipeline reanalysis \
  --workflow ST --chip-number DEMO_ST000006_A1 --casno demo_project \
  --sample-name Demo --tissue brain --chemistry BBV2.4 \
  --species Mus_musculus --method gene_expr --mode strna \
  --sampledir /data/celatlas/results/demo/DEMO_ST000006_A1 \
  --from-step 06 --to-step 08.report --show-steps --dry-run
```

### Demo 7: configuration file

```bash
celatlas count \
  --id DEMO_CONFIG_001 --sample-name Demo_Config --tissue brain \
  --targetdir /data/celatlas/results/demo/DEMO_CONFIG_001 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image gene_expr --fastqs /data/celatlas/fastq/BBV2.4 \
  --config /data/celatlas/config/runner.yaml --dry-run
```

### Demo 8: resume and manifest

```bash
celatlas count \
  --id DEMO_RESUME_001 --sample-name Demo_Resume --tissue brain \
  --targetdir /data/celatlas/results/demo/DEMO_RESUME_001 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image gene_expr --fastqs /data/celatlas/fastq/BBV2.4 \
  --resume-existing --dry-run
```

## 15. Legacy shell entry points

The release still provides:

```text
Celatlas.sh              ST/regular workflow
Celatlas_FFPE.sh         SX/FFPE workflow
Celatlas_SN.sh           SN workflow
Celatlas_reanalysis.sh   downstream reanalysis
Celatlas_run.sh          configuration/scheduler compatibility entry
```

Legacy positional example:

```bash
bash Celatlas.sh DEMO ST110001_A1 project01 BBV2.4 Mus_musculus HE strna \
  --thread 32 --bin 50
```

Use the public CLI for new projects. Shell entry points are useful when migrating an existing scheduler; they still support `--reference_dir`, `--mask_dir`, `--image_dir`, `--fastq_dir`, and `--fastq_name` for shared storage layouts.

## 16. Troubleshooting

### FASTQ not found

`--fastqs` is a directory, not a single file. Check the prefix and the pairs:

```bash
find /data/celatlas/fastq/BBV2.4 -maxdepth 1 -type f | sort | head
```

Do not mix lanes from different samples.

### H&E not found

Use `<id>_he.tif` and point `--image-dir` at its directory. `.tiff`, `.png`, `.jpg`, and `.jpeg` are supported.

### Why does `strna` need a mask?

`strna` maps barcodes onto a spatial canvas and therefore needs `FilterBarcodes.csv` or `barcodeToPos.h5`. `scrna` analyzes an expression matrix without spatial binning.

### `reanalyze` cannot find count_detail

Check `<targetdir>/05.count/<id>_count_detail.txt`. Use `count` when only FASTQ is available, or correct `--targetdir`/sampledir when count outputs are elsewhere.

### Continue after interruption

Repeat the command with `--dry-run --show-resume-plan --preflight`. Use `--resume-existing` only when the FASTQ manifest is unchanged. If inputs changed, use `--reset-fastq-outputs` and optionally `--archive-fastq-outputs`.

### Empty or shifted DAPI segmentation

Use `--cell-segmentation-preset dapi`, keep the manual H&E ROI and registered TIFF, verify that the barcode CSV and labels share the registered coordinate system, and check the blue channel. Do not use `tissue_hires_image.png` as the DAPI inference input.

### Large-image memory pressure

Lower `--thread`, enable `--stardist-tiled-inference`, choose a suitable tile size and overlap, and put the tile cache on a large writable filesystem. `--stardist-no-label-output` reduces output storage but does not remove inference memory requirements.

### Cavity apply behavior

Run `--cavity-mode qc` first. `apply` creates a backup before replacement; use `--cavity-force` only when intentionally regenerating an existing backup.

### Image and spatial plot do not overlap

Check `pixel-size`, barcode coordinates, tissue bbox, image orientation, and registration metadata. Re-run downstream steps with `reanalyze` after correcting registration.

### What to include in a bug report

Keep the command, `celatlas --help`, `python --version`, dry-run output, runner logs, `pipeline.log`, failure report, reference version, model version, and FASTQ manifest. Do not upload FASTQ files, patient data, or private paths to a public issue.

## 17. Performance and reproducibility

- STAR memory is driven by the reference and sorting stage; do not set `--star-thread` to every available core without checking memory.
- bin10, full-resolution labels, and tiled inference increase memory, temporary disk, and I/O.
- Limit threads on shared hosts and set writable `MPLCONFIGDIR`, `NUMBA_CACHE_DIR`, and tile-cache directories.
- Preserve the release archive, `CHECKSUMS.sha256`, configuration, reference version, model version, command line, and FASTQ manifest for every project.
- Keep references, images, FASTQs, and results outside the software release directory so upgrades can be installed side by side.
- Preserve dry-run output, QC images, reports, and parameters when thresholds are changed.

## 18. Known limitations

- Species and panel STAR references are not bundled; users must prepare indexes matching the chemistry and workflow.
- Cell segmentation depends on image quality, registration, model, and tissue type and requires visual QC.
- `gene_expr` infers tissue from expression signal and cannot replace a true tissue image.
- Spatial plots and UMAP/clusters are analysis outputs; cell-type annotation still requires biological interpretation.
- Internal v1.7 production commands, old directory names, and unpublished parameters are not part of the v1.8 public API.

## 19. Related files

- [中文完整手册](celatlas_spatial_manual_zh.md)
- [Linux deployment guide](linux_deployment_guide.md)
- [v1.8 release summary](celatlas-v1.8-release-summary.md)
- [Runner configuration template](../configs/runner.yaml.example)
- [Environment definition](../envs/celatlas18.yml)

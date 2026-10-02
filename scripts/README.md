# scripts/

This directory contains the shell modules used by the Celatlas Spatial v1.8
entrypoints. Keep the root of `scripts/` small and production-facing.

## Production Architecture Modules

### `create_fastq_demo.py`

Creates a small paired FASTQ smoke-test dataset from an authorized input
read pair. It copies complete R1/R2 records, checks that the mates stay in
sync, rewrites every header as `CELATLAS_DEMO:NNNNNNNN/{1,2}`, and writes a
manifest containing only output names, counts, sizes, and SHA256 values. The
original instrument, flow-cell, lane, coordinate, index, and source-path
metadata are not retained.

```bash
python scripts/create_fastq_demo.py \
  --r1 /data/input/sample_R1.fq.gz \
  --r2 /data/input/sample_R2.fq.gz \
  --outdir demo_data/fastq/BBV2.4 \
  --sample CELATLAS_DEMO \
  --reads 20000
```

### `celatlas_config.sh`

Loads portable runtime configuration from one of:

- `$CELATLAS_CONFIG`
- `<repo>/celatlas.local.env`
- `<repo>/configs/celatlas.env`
- `$HOME/.config/celatlas_spatial/celatlas.env`

It defines workspace paths, reference/image/FASTQ/results roots, the runtime
environment name, and optional StarDist model paths. Spatial barcode, image,
and manual mask runtime files are staged under each sample's
`06.segment/mask` directory.

### `celatlas_env.sh`

Shared environment setup for:

- `Celatlas.sh`
- `Celatlas_FFPE.sh`
- `Celatlas_reanalysis.sh`

It exports the runtime Python path, `PYTHONNOUSERSITE`, `MPLCONFIGDIR`,
`NUMBA_CACHE_DIR`, optional `LD_LIBRARY_PATH`, and optional environment `PATH`
hooks. It also supports optional auto-activation through
`CELATLAS_AUTO_ACTIVATE=1`.

### `celatlas_post_binsegment.sh`

Shared post-binSegment option parser and runner for the production entrypoints.
It owns the official v1.8 shell interface for:

- cavity filter options
- optional StarDist cell segmentation options

The main entrypoints call `celatlas_run_post_binsegment_optional_modules`
directly, so standalone legacy wrappers are no longer required in the public
root of `scripts/`.

### `normalize_gff3_to_gtf.py`

Converts plant-style GFF3 annotations to a Celatlas-compatible GTF. It writes
`gene_id` on `gene`, `transcript`, `exon`, `CDS`, and UTR records so STAR,
featureCounts, and Celatlas GTF parsers all see the same gene model.

Use this before `celatlas_spatial rna mkref` when the source annotation is GFF3
or when a converted GTF has exon records without `gene_id`.

```bash
python3 scripts/normalize_gff3_to_gtf.py source.gff3 genes.gtf
```

## Standalone Image Tools

[`background_swap/`](background_swap/README.md) replaces bin/cell display
backgrounds with a matching fluorescence channel after completed DAPI analysis.
It replays and verifies registration, stages previews, and supports backed-up
apply/rollback without changing DAPI inference inputs, labels, or counts.
It does not refresh report HTML or embedded h5ad images.

## Development / Legacy Tools

`scripts/dev/` contains internal migration, debugging, and historical helper
scripts. These files are not required for the production v1.8 architecture and
can be excluded from open-source release packages.

Current examples include:

- `run_cavity_filter.sh`
- `run_post_binsegment.sh`
- `run_stardist_he.py`
- `stardist_assign_counts.py`
- `repair_square_bin_spatial_10x.py`
- `00_Rawdata_celatlas_squarebin.R`
- `test_cavity_filter_SX000131_B1.sh`

Do not add new public user-facing entrypoints under `scripts/dev/`. If a dev
tool becomes part of the supported pipeline, promote the logic into
`celatlas_spatial/tools/` or one of the production shell modules above.

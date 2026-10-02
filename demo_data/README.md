# Celatlas Spatial v1.8 smoke-test data

This directory contains a small paired FASTQ subset for checking installation,
FASTQ discovery, input validation, and command planning. It is not a biological
benchmark and is not expected to produce a meaningful expression matrix from
only 20,000 read pairs.

The read headers were replaced with synthetic IDs (`CELATLAS_DEMO:NNNNNNNN/1`
and `/2`). Bases and quality values in this bundled demo are also synthetic;
only read lengths and pairing structure come from the input. The manifest contains
output names, counts, sizes, and SHA256 values; it does not record the source
server path or original sequencing metadata.

The subset was generated with:

```bash
python scripts/create_fastq_demo.py \
  --r1 /path/to/input_R1.fq.gz \
  --r2 /path/to/input_R2.fq.gz \
  --outdir demo_data/fastq/BBV2.4 \
  --sample CELATLAS_DEMO \
  --reads 20000
```

The commands below exercise the public CLI without requiring a reference,
spatial mask, image, or model. Use `--skip-preflight` because this smoke test
only checks FASTQ discovery and step planning.

```bash
source /path/to/miniforge3/etc/profile.d/conda.sh
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

A real `count` run requires a matching chemistry barcode definition, spatial
barcode files for `strna`, a complete STAR reference, and the relevant image
or tissue mask. The FASTQ subset alone is intentionally only a smoke test.

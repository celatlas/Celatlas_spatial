# Required Runtime Assets

This deployment bundle includes the H&E StarDist model when it was available
at package time. A Swin model is included only when `models/swin_tiny.pth`
exists in the bundle; confirm its inclusion in `RELEASE_NOTES.md`.

The following installation-specific data is not included:

- STAR reference indexes and `celatlas_spatial_genome.config`
- SX/FFPE panel references
- The optional `2D_versatile_fluo` StarDist model for fluorescence or ssDNA
  cell segmentation

Copy `swin_tiny.pth` to the directory configured as `paths.src_dir`. For ST
and SN workflows, configure `paths.reference_dir` so that each species uses:

```text
<reference_dir>/<species>/celatlas_spatial_genome.config
<reference_dir>/<species>/Genome
<reference_dir>/<species>/SA
...
```

For SX, `paths.sx_reference_dir` must point directly at the matching panel
reference directory. Keep reference archives separate from this software bundle
so their large files and genome versions can be tracked independently.

# Celatlas v1.8 Environment Packaging

This directory contains the deployable conda environment definition for
Celatlas Spatial v1.8.

## Primary install file

Use this file for new deployments:

```bash
mamba env create -f envs/celatlas18.yml
conda activate celatlas18
pip install -e . --no-deps
```

## Version source

The environment file is pinned from the validated production stack used for:

- Celatlas Spatial v1.8.0
- Python 3.11.15
- numpy 2.0.2
- pandas 2.3.3
- scanpy 1.10.3
- StarDist 0.9.2
- TensorFlow 2.21.0
- CSBDeep 0.8.2
- PyTorch 2.11.0
- torchvision 0.26.0
- SimpleITK 2.5.5

## Pip version manifest

The pip dependency snapshot lives in:

```text
envs/celatlas18.requirements.txt
```

The repository root `requirements.txt` remains the package install input used
by `setup.py`, while this directory-level file is the deployment-facing alias.

## Recommended validation

```bash
python --version
python -c "import tensorflow; print(tensorflow.__version__)"
python -c "import stardist, csbdeep; print(stardist.__version__, csbdeep.__version__)"
celatlas_spatial rna stardistCellSegment --help
```

## Deployment note

For production deployment, prefer the YAML file in this directory over manual
package installation. Keep the root `requirements.txt` as the pip dependency
input for editable installs and packaging workflows.

# An experimental package to coordinate Chandra x-ray data reduction

This package is not yet ready for production use, and the API is not stable.

## Important Note

I decided to publicly archive this repo and move the efforts towards https://github.com/adonath/snakemake-workflow-chandra instead.


## Installation

The package is not yet released on PyPI, so you need to clone the
repository and install it in development mode:

```bash
git clone https://github.com/adonath/gpxray.git
cd gpxray
```

Create a new conda environment

```bash
conda env create -f environment.yml
conda activate gpxray-dev
```

Install in development mode

```bash
python -m pip install -e .
```


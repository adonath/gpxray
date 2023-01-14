# An experimental package to use Gammapy for x-ray data analysis

This package is a prototype for using Gammapy for X-ray data analysis.
It is not yet ready for production use, and the API is not stable.

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

### Install SAOTrace DB

```bash
curl https://cxc.harvard.edu/cal/Downloads/Hrma/data/saotrace/chandra/orbit-200809-01f-a.tar.gz
tar -xvz -C /data/saotrace/chandra/orbit-200809-01f-a.tar.gz
```


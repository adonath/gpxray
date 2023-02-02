# Licensed under a 3-clause BSD style license - see LICENSE.rst
import sys

import pytest
from gammapy.utils.testing import run_cli

from gpxray.chandra.config import ChandraConfig
from gpxray.cli.core import cli

CONFIG_STR = """
name: my-analysis
sub_name: my-config
obs_ids:
- 62558
obs_id_ref: 62558
roi:
    center:
        frame: icrs
        lon: "06h35m46.5079301472s"
        lat: "-75d16m16.816418256s"
    width: 5.0 arcsec
    energy_range:
        min: 0.5 keV
        max: 7.0 keV
irfs:
    pks-0637:
        spectrum:
            center:
                frame: icrs
                lon: "06h35m46.5079301472s"
                lat: "-75d16m16.816418256s"
            radius: 30 arcsec
            energy_range:
                min: 0.5 keV
                max: 7.0 keV
            energy_groups: 5
            energy_step: 0.01
        psf:
            readout_streak: false
            pileup: false
            extended: true
            minsize: null
"""


@pytest.fixture(scope="session")
def path_config(tmp_path_factory):
    path = tmp_path_factory.mktemp("test-chandra-pks")

    path_config = path / "config.yaml"

    config = ChandraConfig.from_yaml(CONFIG_STR)
    config.write(path_config)

    return path_config


def test_cli_chandra_config(tmp_path):
    path_config = tmp_path / "config.yaml"

    args = ["chandra", f"--filename={path_config}", "init-config"]
    run_cli(cli, args)
    assert path_config.exists()


def test_cli_chandra_download(path_config):
    args = ["chandra", f"--filename={path_config}", "download", "--exclude", "vvref"]
    run_cli(cli, args)

    path = path_config.parent / "data/62558/"
    assert path.exists()


def test_cli_chandra_reprocess(path_config):
    args = [
        "chandra",
        f"--filename={path_config}",
        "reprocess",
    ]
    run_cli(cli, args)

    path = path_config.parent / "data/62558/repro"
    assert path.exists()


def test_cli_chandra_reproject_events(path_config):
    args = [
        "chandra",
        f"--filename={path_config}",
        "reproject-events",
    ]
    run_cli(cli, args)

    path = (
        path_config.parent / "data/62558/repro/acisf62558_repro_evt2_reprojected.fits"
    )
    assert path.exists()


def test_cli_chandra_bin_events(path_config):
    args = [
        "chandra",
        f"--filename={path_config}",
        "bin-events",
    ]
    run_cli(cli, args)

    path = path_config.parent / "my-config/62558/counts.fits.gz"
    assert path.exists()


def test_cli_chandra_extract_spectra(path_config):
    args = [
        "chandra",
        f"--filename={path_config}",
        "extract-spectra",
    ]
    run_cli(cli, args)

    path = path_config.parent / "my-config/62558/spectrum/pks-0637/pks-0637.pi"
    assert path.exists()


def test_cli_chandra_fit_spectra(path_config):
    args = [
        "chandra",
        f"--filename={path_config}",
        "fit-spectra",
    ]
    run_cli(cli, args)

    path = (
        path_config.parent
        / "my-config/62558/spectrum/pks-0637/source-flux-chart-pks-0637.dat"
    )
    assert path.exists()

    path = (
        path_config.parent
        / "my-config/62558/spectrum/pks-0637/spectral-fit-model-pks-0637.yaml"
    )
    assert path.exists()


def test_cli_chandra_compute_exposure(path_config):
    args = [
        "chandra",
        f"--filename={path_config}",
        "compute-exposure",
    ]
    run_cli(cli, args)

    path = path_config.parent / "my-config/62558/exposure.fits.gz"
    assert path.exists()


@pytest.mark.skipif(sys.platform.startswith("darwin"), reason="Skip on MacOS")
def test_cli_chandra_simulate_psf(path_config):
    args = [
        "chandra",
        f"--filename={path_config}",
        "simulate-psf",
    ]
    run_cli(cli, args)

    path = path_config.parent / "my-config/62558/psf-marx-pks-0637.fits.gz"
    assert path.exists()

# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest
from gammapy.utils.testing import run_cli

from gpxray.chandra.config import ChandraConfig
from gpxray.cli.core import cli

CONFIG_STR = """
name: my-analysis
sub_name: my-config
obs_ids:
- 1093
obs_id_ref: 1093
roi:
    center:
        frame: icrs
        lon: "06h35m46.5079301472s"
        lat: "-75d16m16.816418256s"
    width: 5.0 arcsec
energy_range:
    min: 0.5 keV
    max: 7.0 keV
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

    args = [
        "chandra",
        f"--filename={path_config}",
        "download",
    ]
    run_cli(cli, args)

    path = path_config.parent / "data/1093/"
    assert path.exists()


def test_cli_chandra_reprocess(path_config):
    args = [
        "chandra",
        f"--filename={path_config}",
        "reprocess",
    ]
    run_cli(cli, args)

    path = path_config.parent / "data/1093/repro"
    assert path.exists()

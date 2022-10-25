# Licensed under a 3-clause BSD style license - see LICENSE.rst
import logging
import warnings
from pathlib import Path

import click

from gpxray import __version__
from gpxray.chandra.config import ChandraConfig
from gpxray.chandra.io import ChandraFileIndex

from .chandra import (
    cli_chandra_download,
    cli_chandra_init_config,
    cli_chandra_reprocess,
    cli_chandra_reproject_events,
)


class ContextObject(object):
    """Context"""

    def __init__(self, config, filename, obs_id=-1, overwrite=False):
        self.config = config

        if obs_id == -1:
            obs_ids = config.obs_ids
        else:
            obs_ids = [obs_id]

        self.obs_ids = obs_ids
        self.overwrite = overwrite
        self.filename = Path(filename)

    @property
    def file_indices(self):
        """File indices"""
        return [
            ChandraFileIndex(obs_id=obs_id, path=self.filename.parent)
            for obs_id in self.obs_ids
        ]

    @property
    def file_index_ref(self):
        """Reference file index"""
        return ChandraFileIndex(
            obs_id=self.config.obs_id_ref, path=self.filename.parent
        )


# We implement the --version following the example from here:
# http://click.pocoo.org/5/options/#callbacks-and-eager-options
def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print(f"gpxray version {__version__}")
    ctx.exit()


# http://click.pocoo.org/5/documentation/#help-parameter-customization
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

# https://click.palletsprojects.com/en/5.x/python3/#unicode-literals
click.disable_unicode_literals_warning = True


@click.group("gpxray", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--log-level",
    default="info",
    help="Logging verbosity level.",
    type=click.Choice(["debug", "info", "warning", "error"]),
)
@click.option("--ignore-warnings", is_flag=True, help="Ignore warnings?")
@click.option(
    "--version",
    is_flag=True,
    callback=print_version,
    expose_value=False,
    is_eager=True,
    help="Print version and exit.",
)
def cli(log_level, ignore_warnings):  # noqa: D301
    """gpxray command line interface (CLI).

    \b
    Examples
    --------

    \b
    $ gpxray --help
    $ gpxray --version
    """
    logging.basicConfig(level=log_level.upper())

    if ignore_warnings:
        warnings.simplefilter("ignore")


@cli.group("chandra")
@click.option("--filename", help="Config file name", default="config.yaml")
@click.option("--obs-id", help="Obs ID", default=-1)
@click.option(
    "--overwrite", default=False, is_flag=True, help="Overwrite existing files."
)
@click.pass_context
def cli_chandra(ctx, filename, obs_id, overwrite):
    """Chandra sub-commmands"""
    filename = Path(filename)

    if filename.exists():
        config = ChandraConfig.read(path=filename)
    else:
        config = ChandraConfig()

    ctx.obj = ContextObject(
        config=config, obs_id=obs_id, overwrite=overwrite, filename=filename
    )


cli_chandra.add_command(cli_chandra_init_config)
cli_chandra.add_command(cli_chandra_download)
cli_chandra.add_command(cli_chandra_reprocess)
cli_chandra.add_command(cli_chandra_reproject_events)

"""Run Chandra data preparation for jolideco"""
import logging
import subprocess

import click

from gpxray.chandra.config import ChandraConfig

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)


def execute_command(command, cwd="."):
    log.info(f"Executing: {' '.join(command)}")
    subprocess.run(command, cwd=cwd)


@click.command(name="init-config")
@click.option(
    "--filename",
    default="config.yaml",
    help="Filename to store the default configuration values.",
    show_default=True,
)
@click.option(
    "--overwrite", default=False, is_flag=True, help="Overwrite existing file."
)
def cli_chandra_init_config(filename, overwrite):
    """Writes default configuration file."""
    config = ChandraConfig()
    config.write(filename, overwrite=overwrite)
    log.info(f"Configuration file produced: {filename}")

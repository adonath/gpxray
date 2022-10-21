"""Run Chandra data preparation for jolideco"""
import logging
import subprocess

import click

from gpxray.chandra.config import ChandraConfig

log = logging.getLogger(__name__)


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
    log.info(f"Writing: {filename}")


@click.command(name="download")
@click.option(
    "--filename",
    default="config.yaml",
    help="Filename to read configuration from.",
    show_default=True,
)
@click.option(
    "--overwrite", default=False, is_flag=True, help="Overwrite existing file."
)
def cli_chandra_download(filename, overwrite):
    """Download data"""
    config = ChandraConfig.read(filename)

    for index in config.file_indices:
        index.path_data.mkdir(exist_ok=True)

        if not index.path_obs_id.exists() or overwrite:
            command = ["download_chandra_obsid", f"{index.obs_id}"]
            execute_command(command=command, cwd=f"{index.path_data}")
        else:
            log.info(f"Skipping download, {index.path_obs_id} already exists.")

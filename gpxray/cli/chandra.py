"""Run Chandra data preparation for jolideco"""
import logging
import subprocess

import click

from gpxray.chandra.utils import run_ciao_tool

log = logging.getLogger(__name__)


def execute_command(command, cwd="."):
    log.info(f"Executing: {' '.join(command)}")
    subprocess.run(command, cwd=cwd)


@click.command(name="init-config", short_help="Init config file")
@click.pass_obj
def cli_chandra_init_config(obj):
    """Writes default configuration file."""
    obj.config.write(obj.filename, overwrite=obj.overwrite)
    log.info(f"Writing: {obj.filename}")


@click.command(name="download", short_help="Download chandra observation data")
@click.pass_obj
def cli_chandra_download(obj):
    """Download data"""
    for index in obj.file_indices:
        if index.path_obs_id.exists() and not obj.overwrite:
            log.info(f"Skipping download, {index.path_obs_id} already exists.")
            continue

        index.path_data.mkdir(exist_ok=True)

        command = ["download_chandra_obsid", f"{index.obs_id}"]
        execute_command(command=command, cwd=f"{index.path_data}")


@click.command("reprocess", short_help="Reprocess chandra observation data")
@click.pass_obj
def cli_chandra_reprocess(obj):
    """Reprocess data"""
    for index in obj.file_indices:
        if index.path_repro.exists() and not obj.overwrite:
            log.info(f"Skipping reprocessing, {index.path_repro} already exists.")
            continue

        run_ciao_tool("chandra_repro", config=obj.config, file_index=index)

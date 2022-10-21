"""Run Chandra data preparation for jolideco"""
import logging
import subprocess

import click

log = logging.getLogger(__name__)


def execute_command(command, cwd="."):
    log.info(f"Executing: {' '.join(command)}")
    subprocess.run(command, cwd=cwd)


@click.command(name="init-config")
@click.pass_obj
def cli_chandra_init_config(obj):
    """Writes default configuration file."""
    obj.config.write(obj.filename, overwrite=obj.overwrite)
    log.info(f"Writing: {obj.filename}")


@click.command(name="download")
@click.pass_obj
def cli_chandra_download(obj):
    """Download data"""
    for index in obj.file_indices:
        index.path_data.mkdir(exist_ok=True)

        if not index.path_obs_id.exists() or obj.overwrite:
            command = ["download_chandra_obsid", f"{index.obs_id}"]
            execute_command(command=command, cwd=f"{index.path_data}")
        else:
            log.info(f"Skipping download, {index.path_obs_id} already exists.")

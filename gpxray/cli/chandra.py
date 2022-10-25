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
@click.option(
    "-e",
    "--exclude",
    type=click.STRING,
    help="Sub selection of data to download",
    default="",
)
@click.pass_obj
def cli_chandra_download(obj, exclude):
    """Download data"""
    for index in obj.file_indices:
        if index.path_obs_id.exists() and not obj.overwrite:
            log.info(f"Skipping download, {index.path_obs_id} already exists.")
            continue

        index.path_data.mkdir(exist_ok=True)

        command = ["download_chandra_obsid", f"{index.obs_id}", "--exclude", exclude]
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


@click.command("reproject-events", short_help="Reproject events to common WCS")
@click.pass_obj
def cli_chandra_reproject_events(obj):
    """Reproject events"""
    index_ref = obj.file_index_ref

    for index in obj.file_indices:
        if index.filename_repro_evt2_reprojected.exists() and not obj.overwrite:
            log.info(
                f"Skipping reproject events, {index.filename_repro_evt2_reprojected} "
                "already exists."
            )
            continue

        run_ciao_tool(
            "reproject_events",
            config=obj.config,
            file_index=index,
            file_index_ref=index_ref,
        )

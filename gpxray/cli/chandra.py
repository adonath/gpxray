"""Run Chandra data preparation for jolideco"""
import logging
import subprocess

import click
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from gpxray.chandra.utils import run_ciao_tool, run_sherpa_spectral_fit

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
    default="vvref",
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

        run_ciao_tool(
            config=obj.config.ciao.chandra_repro,
            file_index=index,
            clobber=obj.overwrite,
        )


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
            config=obj.config.ciao.reproject_events,
            file_index=index,
            file_index_ref=index_ref,
            clobber=obj.overwrite,
        )


@click.command("bin-events", short_help="Bin events into FITS image")
@click.pass_obj
def cli_chandra_bin_events(obj):
    """Bin events"""
    for index in obj.file_indices:
        if index.filename_counts.exists() and not obj.overwrite:
            log.info(f"Skipping bin events, {index.filename_counts} " "already exists.")
            continue

        run_ciao_tool(
            config=obj.config.roi,
            file_index=index,
            clobber=obj.overwrite,
        )


@click.command("compute-exposure", short_help="Compute exposure FITS image")
@click.pass_obj
def cli_chandra_compute_exposure(obj):
    """Compute exposure image"""
    # TODO: take into account spatial dependence and maybe compute absolute exposure...
    exposure_ref = obj.file_index_ref.index_table.meta["EXPOSURE"]

    for index in obj.file_indices:
        if index.filename_exposure.exists() and not obj.overwrite:
            log.info(
                f"Skipping compute exposure, {index.filename_exposure} "
                "already exists."
            )
            continue

        value = index.index_table.meta["EXPOSURE"]
        header = fits.getheader(index.filename_counts)
        shape = header["NAXIS2"], header["NAXIS1"]
        data = value * np.ones(shape) / exposure_ref
        hdu = fits.PrimaryHDU(data=data, header=WCS(header).to_header())

        hdulist = fits.HDUList([hdu])

        filename = index.filename_exposure
        log.info(f"Writing {filename}")

        hdulist.writeto(filename, overwrite=obj.overwrite)


@click.command("extract-spectra", short_help="Extract spectra, arfs and rmfs")
@click.pass_obj
def cli_chandra_extract_spectra(obj):
    """Extract spectra"""
    for index in obj.file_indices:
        for name, irf_config in obj.config.irfs.items():
            filename_spectrum = index.filenames_spectra[name]

            if filename_spectrum.exists() and not obj.overwrite:
                log.info(
                    f"Skipping extact spectrum, {filename_spectrum} " "already exists."
                )
                continue

            run_ciao_tool(
                config=irf_config.spectrum,
                file_index=index,
                irf_label=name,
                clobber=obj.overwrite,
            )


@click.command("fit-spectra", short_help="Fit spectra")
@click.pass_obj
def cli_chandra_fit_spectra(obj):
    """Fit spectra"""
    for index in obj.file_indices:
        for name, config_irf in obj.config.irfs.items():
            filename_spectrum = index.filenames_spectra[name]

            if filename_spectrum.exists() and not obj.overwrite:
                log.info(f"Skipping fit spectrum, {filename_spectrum} already exists.")
                continue

            run_sherpa_spectral_fit(
                config_irf=config_irf.spectrum, file_index=index, irf_label=name
            )


def copy_file(path_input, path_output):
    """Copy file from path input to output"""
    command = ["cp", f"{path_input}", f"{path_output}"]
    execute_command(command=command)


@click.command("simulate-psf", short_help="Simulate PSF FITS image")
@click.pass_obj
def cli_chandra_simulate_psf(obj):
    """Simulate psf"""
    for index in obj.file_indices:
        for name, irf_config in obj.config.irfs.items():
            filename_psf = index.filenames_psf[name]

            if filename_psf.exists() and not obj.overwrite:
                log.info(f"Skipping simulate-psf, {filename_psf} " "already exists.")
                continue

            run_ciao_tool(
                config=irf_config.psf,
                file_index=index,
                irf_label=name,
                clobber=obj.overwrite,
            )
            path_input = index.paths_psf[name] / "psf"
            copy_file(path_input=path_input, path_output=filename_psf)


@click.command("all", short_help="Run all commands")
@click.pass_context
def cli_chandra_all(ctx):
    """Run all commands"""
    ctx.forward(cli_chandra_download)
    ctx.forward(cli_chandra_reprocess)
    ctx.forward(cli_chandra_reproject_events)
    ctx.forward(cli_chandra_bin_events)
    ctx.forward(cli_chandra_compute_exposure)
    ctx.forward(cli_chandra_extract_spectra)
    ctx.forward(cli_chandra_fit_spectra)
    ctx.forward(cli_chandra_simulate_psf)

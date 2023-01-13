"""Run Chandra data preparation for jolideco"""
import logging
import subprocess

import click
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from gpxray.chandra.config import PSFSimulatorEnum
from gpxray.chandra.utils import run_ciao_tool, run_sao_trace, run_sherpa_spectral_fit

log = logging.getLogger(__name__)


def execute_command(command, cwd="."):
    log.info(f"Executing: {' '.join(command)}")
    subprocess.run(command, cwd=cwd)


def gzip_fits_file(filename):
    """Gzip file using Astropy"""
    hdulist = fits.open(filename)

    log.info(f"GZipping {filename}")

    hdulist.writeto(filename.parent / (filename.name + ".gz"), overwrite=True)

    filename.unlink()


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
    for file_index in obj.file_indices:
        if file_index.path_obs_id.exists() and not obj.overwrite:
            log.info(f"Skipping download, {file_index.path_obs_id} already exists.")
            continue

        file_index.path_data.mkdir(exist_ok=True)

        command = [
            "download_chandra_obsid",
            f"{file_index.obs_id}",
            "--exclude",
            exclude,
        ]
        execute_command(command=command, cwd=f"{file_index.path_data}")


@click.command("reprocess", short_help="Reprocess chandra observation data")
@click.pass_obj
def cli_chandra_reprocess(obj):
    """Reprocess data"""
    for file_index in obj.file_indices:
        if file_index.path_repro.exists() and not obj.overwrite:
            log.info(f"Skipping reprocessing, {file_index.path_repro} already exists.")
            continue

        run_ciao_tool(
            config=obj.config.ciao.chandra_repro,
            file_index=file_index,
            overwrite=obj.overwrite,
        )


@click.command("reproject-events", short_help="Reproject events to common WCS")
@click.pass_obj
def cli_chandra_reproject_events(obj):
    """Reproject events"""
    file_index_ref = obj.file_index_ref

    for file_index in obj.file_indices:
        if file_index.filename_repro_evt2_reprojected.exists() and not obj.overwrite:
            log.info(
                f"Skipping reproject events, {file_index.filename_repro_evt2_reprojected} "
                "already exists."
            )
            continue

        run_ciao_tool(
            config=obj.config.ciao.reproject_events,
            file_index=file_index,
            file_index_ref=file_index_ref,
            overwrite=obj.overwrite,
        )


@click.command("bin-events", short_help="Bin events into FITS image")
@click.pass_obj
def cli_chandra_bin_events(obj):
    """Bin events"""
    for file_index in obj.file_indices:
        if file_index.filename_counts.exists() and not obj.overwrite:
            log.info(
                f"Skipping bin events, {file_index.filename_counts} " "already exists."
            )
            continue

        run_ciao_tool(
            config=obj.config.roi,
            file_index=file_index,
            overwrite=obj.overwrite,
        )
        gzip_fits_file(filename=file_index.filename_counts)


def compute_exposure_simple(file_index, obj):
    """Compute exposure simple"""
    exposure_ref = obj.file_index_ref.index_table.meta["EXPOSURE"]
    value = file_index.index_table.meta["EXPOSURE"]

    filename = file_index.filename_counts
    header = fits.getheader(filename.parent / (filename.name + ".gz"))

    shape = header["NAXIS2"], header["NAXIS1"]
    data = value * np.ones(shape) / exposure_ref
    hdu = fits.PrimaryHDU(data=data, header=WCS(header).to_header())

    hdulist = fits.HDUList([hdu])

    filename = file_index.filename_exposure
    log.info(f"Writing {filename}")

    hdulist.writeto(filename, overwrite=obj.overwrite)


def compute_exposure_ciao(file_index, obj):
    """Compute exposure ciao"""
    irf_label = list(obj.config.irfs)[0]

    if file_index.filename_repro_asp_hist.exists() and not obj.overwrite:
        log.info(
            f"Skipping compute asp hist, {file_index.filename_repro_asp_hist} "
            "already exists."
        )
        pass
    else:
        run_ciao_tool(
            config=obj.config.ciao.asphist,
            file_index=file_index,
            overwrite=obj.overwrite,
        )

    if file_index.filename_repro_inst_map.exists() and not obj.overwrite:
        log.info(
            f"Skipping compute inst map, {file_index.filename_repro_inst_map} "
            "already exists."
        )
        pass
    else:
        run_ciao_tool(
            config=obj.config.irfs[irf_label].aeff,
            file_index=file_index,
            overwrite=obj.overwrite,
            irf_label=irf_label,
        )

    run_ciao_tool(
        config=obj.config.irfs[irf_label].exposure,
        file_index=file_index,
        overwrite=obj.overwrite,
    )


@click.command("compute-exposure", short_help="Compute exposure FITS image")
@click.pass_obj
def cli_chandra_compute_exposure(obj):
    """Compute exposure image"""
    # TODO: take into account spatial dependence and maybe compute absolute exposure...

    for file_index in obj.file_indices:

        if file_index.filename_exposure.exists() and not obj.overwrite:
            log.info(
                f"Skipping compute exposure, {file_index.filename_exposure} "
                "already exists."
            )
            continue

        compute_exposure_simple(file_index=file_index, obj=obj)
        gzip_fits_file(file_index.filename_exposure)


@click.command("extract-spectra", short_help="Extract spectra, arfs and rmfs")
@click.pass_obj
def cli_chandra_extract_spectra(obj):
    """Extract spectra"""
    for file_index in obj.file_indices:
        for irf_label, irf_config in obj.config.irfs.items():
            filename_spectrum = file_index.filenames_spectra[irf_label]

            if filename_spectrum.exists() and not obj.overwrite:
                log.info(
                    f"Skipping extact spectrum, {filename_spectrum} " "already exists."
                )
                continue

            run_ciao_tool(
                config=irf_config.spectrum,
                file_index=file_index,
                irf_label=irf_label,
                overwrite=obj.overwrite,
            )


@click.command("fit-spectra", short_help="Fit spectra")
@click.pass_obj
def cli_chandra_fit_spectra(obj):
    """Fit spectra"""
    for file_index in obj.file_indices:
        for irf_label, config_irf in obj.config.irfs.items():
            filename_spectrum = file_index.filenames_spectra[irf_label]

            if filename_spectrum.exists() and not obj.overwrite:
                log.info(f"Skipping fit spectrum, {filename_spectrum} already exists.")
                continue

            run_sherpa_spectral_fit(
                config=config_irf.spectrum,
                file_index=file_index,
                irf_label=irf_label,
                overwrite=obj.overwrite,
            )


def copy_file(path_input, path_output):
    """Copy file from path input to output"""
    command = ["cp", f"{path_input}", f"{path_output}"]
    execute_command(command=command)


@click.command("simulate-psf", short_help="Simulate PSF FITS image")
@click.pass_obj
def cli_chandra_simulate_psf(obj):
    """Simulate psf"""
    for file_index in obj.file_indices:
        for irf_label, irf_config in obj.config.irfs.items():
            filename_psf = file_index.filenames_psf[irf_label]

            if filename_psf.exists() and not obj.overwrite:
                log.info(f"Skipping simulate-psf, {filename_psf} " "already exists.")
                continue

            if obj.config.psf_simulator == PSFSimulatorEnum.saotrace:
                for idx in range(irf_config.psf.numiter):
                    run_sao_trace(
                        saotrace_config=obj.config.saotrace,
                        file_index=file_index,
                        irf_config=irf_config,
                        irf_label=irf_label,
                        idx=idx,
                    )

            run_ciao_tool(
                config=irf_config.psf,
                file_index=file_index,
                irf_label=irf_label,
                overwrite=obj.overwrite,
            )
            path_input = file_index.paths_psf_marx[irf_label] / "psf"
            copy_file(path_input=path_input, path_output=filename_psf)
            gzip_fits_file(filename=filename_psf)


@click.command("all", short_help="Run all commands")
@click.pass_context
def cli_chandra_all(ctx):
    """Run all commands"""
    ctx.forward(cli_chandra_download)
    ctx.forward(cli_chandra_reprocess)
    ctx.forward(cli_chandra_reproject_events)
    ctx.forward(cli_chandra_bin_events)
    ctx.forward(cli_chandra_extract_spectra)
    ctx.forward(cli_chandra_fit_spectra)
    ctx.forward(cli_chandra_compute_exposure)
    ctx.forward(cli_chandra_simulate_psf)

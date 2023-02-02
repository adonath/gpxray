import logging
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import sherpa.astro.ui as sau
from ciao_contrib import runtool
from sherpa_contrib.chart import save_chart_spectrum

from gpxray.chandra.io import convert_spectrum_chart_to_rdb, write_sherpa_model_to_yaml

log = logging.getLogger(__name__)

RUN_SAO_SCRIPT_PATH = Path(__file__).parent.parent.parent / "scripts/run-saotrace.sh"


def run_ciao_tool(
    config, file_index, file_index_ref=None, irf_label=None, overwrite=False
):
    """Run ciao tool

    Parameters
    ----------
    config : `~gpxray.chandra.config.CiaoBaseConfig`
        Tools config
    file_index : `ChandraFileIndex`
        Chandra file index
    file_index_ref :  `ChandraFileIndex`
        Reference file index
    irf_label : str
        IRF label
    overwrite ; bool
        Overwrite output files
    """
    with runtool.new_pfiles_environment(ardlib=False):
        tool = getattr(runtool, config._tool_name)

        kwargs = config.to_ciao(
            file_index=file_index, file_index_ref=file_index_ref, irf_label=irf_label
        )

        if "clobber" in kwargs:
            kwargs["clobber"] = overwrite

        tool.punlearn()
        tool(**kwargs)


def run_sherpa_spectral_fit(config, file_index, irf_label, overwrite, pileup=True):
    """Run sherpa spectral fit

    Parameters
    ----------
    config : `~gpxray.chandra.config.PerSourceSpecExtractConfig`
        Tools config
    file_index : `ChandraFileIndex`
        Chandra file index
    file_index_ref :  `ChandraFileIndex`
        Reference file index
    irf_label : str
        IRF label
    overwrite ; bool
        Overwrite output files
    """
    filename_pha = file_index.paths_spectra_pha[irf_label] / f"{irf_label}.pi"
    sau.load_data(str(filename_pha))

    sau.group_counts(config.energy_groups)

    e_min = config.energy_range.min.to_value("keV")
    e_max = config.energy_range.max.to_value("keV")
    sau.notice(e_min, e_max)

    sau.set_stat("cash")

    sau.set_source(sau.xsphabs.absorption * sau.powlaw1d.pwl)
    sau.xsphabs.absorption.nh.val = 1
    sau.powerlaw.pwl.ampl.val = 1e-5
    sau.powerlaw.pwl.index.val = 1.5

    if pileup:
        sau.set_pileup_model(sau.jdpileup.jdp)
        sau.jdpileup.jdp.f.min = 0.85
        sau.jdpileup.jdp.ftime = file_index.exptime
        sau.jdpileup.jdp.fracexp = 0.987

    sau.guess(sau.powlaw1d.pwl)
    sau.fit()

    # sau.plot_fit_delchi()
    filename = file_index.filenames_spectral_fit_png[irf_label]
    log.info(f"Writing {filename}")
    plt.savefig(filename, dpi=300)

    sau.set_analysis(1, "energy", "rate", factor=1)

    filename = file_index.filenames_spectra[irf_label]

    save_chart_spectrum(str(filename), elow=e_min, ehigh=e_max, clobber=overwrite)
    convert_spectrum_chart_to_rdb(filename, overwrite=overwrite)

    filename = file_index.filenames_spectral_fit_model[irf_label]
    write_sherpa_model_to_yaml(filename, overwrite=overwrite)


def run_sao_trace(
    saotrace_config, irf_config, file_index, irf_label, idx, use_docker=True
):
    """Run sao trace

    Parameters
    ----------
    config : `~gpxray.chandra.config.SAOTraceConfig`
        Tools config
    file_index : `ChandraFileIndex`
        Chandra file index
    irf_label : str
        IRF label
    idx : int
        Iteration index
    use_docker : bool
        Use docker command
    """
    path = file_index.paths_psf_saotrace[irf_label]
    filename = path / "src.lua"

    with filename.open("w") as fh:
        text = saotrace_config.to_src_pars(
            file_index=file_index, irf_config=irf_config, irf_label=irf_label
        )
        log.info(f"Writing {filename}")
        fh.write(text)

    command = ["trace-nest"]
    command += saotrace_config.to_trace_nest_config(
        file_index=file_index, irf_label=irf_label, idx=idx
    )

    if use_docker:
        path = RUN_SAO_SCRIPT_PATH.absolute()
        command = ["sh", f"{path}", "-l", "docker"] + command

    log.info(f"Running command: {' '.join(command)}")

    subprocess.call(command)

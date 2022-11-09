import logging

import matplotlib.pyplot as plt
import sherpa.astro.ui as sau
from ciao_contrib import runtool
from sherpa_contrib.chart import save_chart_spectrum

log = logging.getLogger(__name__)


def run_ciao_tool(tool_name, config, file_index, file_index_ref=None, irf_label=None):
    """Run ciao tool

    Parameters
    ----------
    tool_name : str
        Tool name to run
    config : `~gpxtay.chandra.config.CiaoToolsConfig`
        Tools config
    file_index : `ChandraFileIndex`
        Chandra file index
    file_index_ref :  `ChandraFileIndex`
        Reference file index
    """
    with runtool.new_pfiles_environment(ardlib=True):
        tool = getattr(runtool, tool_name)
        tool_config = getattr(config.ciao, tool_name)

        kwargs = tool_config.to_ciao(
            file_index=file_index, file_index_ref=file_index_ref, irf_label=irf_label
        )

        if tool_name == "dmcopy":
            selection = f"[EVENTS][{config.roi.to_ciao(wcs=file_index.wcs)}]"
            selection += f"[{config.energy_range.to_ciao()}]"
            kwargs["infile"] += selection

        tool.punlearn()
        tool(**kwargs)


def run_sherpa_spectral_fit(config_irf, file_index, irf_label):
    """Run sherpa spectral fit"""

    filename_pha = file_index.paths_spectra_pha[irf_label] / f"{irf_label}.pi"
    sau.load_data(str(filename_pha))

    sau.group_counts(config_irf.energy_groups)

    e_min = config_irf.energy_range.min.to_value("keV")
    e_max = config_irf.energy_range.max.to_value("keV")
    sau.notice(e_min, e_max)

    sau.set_source(sau.xsphabs.absorption * sau.powerlaw.pwl)
    sau.xsphabs.absorption.nh = 0.5
    sau.powerlaw.pwl.ampl.val = 1e-10
    sau.powerlaw.pwl.index.val = -1.5

    sau.fit()
    sau.set_analysis(1, "energy", "rate", factor=1)
    sau.plot_data()
    sau.plot_source(overplot=True)
    sau.plot_model(overplot=True)
    plt.xlim(0.3 * e_min, 3 * e_max)

    plt.xscale("log")
    plt.yscale("log")

    filename = file_index.filenames_spectra_png[irf_label]
    log.info(f"Writing {filename}")
    plt.savefig(filename)

    filename = file_index.filenames_spectra[irf_label]
    save_chart_spectrum(str(filename), elow=e_min, ehigh=e_max)

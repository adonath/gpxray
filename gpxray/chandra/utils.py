import logging

import sherpa.astro.ui as sau
from ciao_contrib import runtool
from sherpa_contrib.chart import save_chart_spectrum

log = logging.getLogger(__name__)


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


def run_sherpa_spectral_fit(config, file_index, irf_label, overwrite):
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

    sau.set_source(sau.xsphabs.absorption * sau.powerlaw.pwl)
    sau.xsphabs.absorption.nh = 0.5
    sau.powerlaw.pwl.ampl.val = 1e-10
    sau.powerlaw.pwl.index.val = -1.5

    sau.fit()
    sau.set_analysis(1, "energy", "rate", factor=1)

    filename = file_index.filenames_spectra[irf_label]

    save_chart_spectrum(str(filename), elow=e_min, ehigh=e_max, clobber=overwrite)

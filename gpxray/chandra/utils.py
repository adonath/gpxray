from ciao_contrib import runtool


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

        if tool_name == "specextract":
            selection = config.to_ciao_spec_extract(wcs=file_index.wcs)
            kwargs["infile"] += selection

        tool.punlearn()
        tool(**kwargs)

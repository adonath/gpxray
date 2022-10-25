from ciao_contrib import runtool


def to_ciao_name(name):
    """Convert parameter name to ciao name"""
    return name.replace("_", "-")


def run_ciao_tool(tool_name, config, file_index, file_index_ref=None):
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

        kwargs = tool_config.dict()

        for name in tool_config.required_names:
            ciao_name = to_ciao_name(name)

            kwargs[ciao_name] = kwargs[name].format(
                file_index=file_index, file_index_ref=file_index_ref
            )

        tool.punlearn()
        tool(**kwargs)

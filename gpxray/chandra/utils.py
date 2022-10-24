from ciao_contrib import runtool


def run_ciao_tool(tool_name, config, file_index):
    """Run ciao tool

    Parameters
    ----------
    tool_name : str
        Tool name to run
    config : `~gpxtay.chandra.config.CiaoToolsConfig`
        Tools config
    """
    tool = getattr(runtool, tool_name)
    tool_config = getattr(config.ciao, tool_name)

    kwargs = tool_config.dict()

    for name in tool_config.required_names:
        kwargs[name] = kwargs[name].format(file_index=file_index)

    tool(**kwargs)

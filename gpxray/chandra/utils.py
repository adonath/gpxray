from ciao_contrib import runtool


def run_ciao_tool(tool_name, config):
    """Run ciao tool

    Parameters
    ----------
    tool_name : str
        Tool name to run
    config : `~gpxtay.chandra.config.CiaoToolsConfig`
        Tools config
    """
    tool = getattr(runtool, tool_name)
    tool_config = getattr(config, tool_name).to_dict()
    tool(**tool_config)

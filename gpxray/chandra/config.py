import json
import logging
from typing import List

import numpy as np
import yaml
from astropy import units as u
from astropy.coordinates import Angle
from astropy.time import Time
from astropy.units import Quantity
from ciao_contrib import runtool
from gammapy.analysis.config import AngleType, EnergyType, FrameEnum, GammapyBaseConfig
from gammapy.utils.scripts import make_path, read_yaml
from pydantic import BaseModel, create_model

log = logging.getLogger(__name__)


CIAO_TOOLS_TYPES = {"f": str, "i": int, "s": str, "b": bool, "r": float}


GammapyBaseConfig.Config.json_encoders[Angle] = lambda v: "None"


class BaseConfig(BaseModel):
    """Base config"""

    class Config:
        validate_all = True
        validate_assignment = True
        extra = "forbid"
        json_encoders = {
            Angle: lambda v: v.to_string(),
            Quantity: lambda v: f"{v.value} {v.unit}",
            Time: lambda v: f"{v.value}",
        }

    @property
    def required_names(self):
        """Get required field names"""
        return [field.name for field in self.__fields__.values() if field.required]


def create_ciao_config(toolname, model_name):
    """Create config class"""
    par_info = runtool.parinfo[toolname]

    parameters = {}

    for par in par_info["req"]:
        parameters[par.name] = (CIAO_TOOLS_TYPES[par.type], Ellipsis)

    for par in par_info["opt"]:
        parameters[par.name] = (CIAO_TOOLS_TYPES[par.type], par.default)

    model = create_model(model_name, __base__=BaseConfig, **parameters)

    return model


DMCopyConfig = create_ciao_config("dmcopy", "DMCopyConfig")
ChandraReproConfig = create_ciao_config("chandra_repro", "ChandraReproConfig")
ReprojectEventsConfig = create_ciao_config("reproject_events", "ReprojectEventsConfig")
SimulatePSFConfig = create_ciao_config("simulate_psf", "SimulatePSFConfig")


class CiaoToolsConfig(BaseConfig):
    dmcopy: DMCopyConfig = DMCopyConfig(infile="{file_index.}", outfile="{file_index.}")
    chandra_repro: ChandraReproConfig = ChandraReproConfig(
        indir="{file_index.path_obs_id}", outdir="{file_index.path_repro}"
    )
    reproject_events: ReprojectEventsConfig = ReprojectEventsConfig(
        infile="{file_index.filename_repro_evt2}",
        outfile="{file_index.filename_repro_evt2_reprojected}",
        match="{file_index_ref.filename_repro_evt2}}",
    )
    simulate_psf: SimulatePSFConfig = SimulatePSFConfig(
        infile="{file_index.path_obs_id}",
        outroot="{file_index.path_repro}",
        ra=np.nan,
        dec=np.nan,
        spectrumfile="",
    )


class SkyCoordConfig(BaseConfig):
    frame: FrameEnum = FrameEnum.icrs
    lon: AngleType = Angle("06h35m46.5079301472s")
    lat: AngleType = Angle("-75d16m16.816418256s")


class ROIConfig(BaseConfig):
    center: SkyCoordConfig = SkyCoordConfig()
    width: AngleType = Angle("5 arcsec")


class EnergyRangeConfig(BaseConfig):
    min: EnergyType = 0.5 * u.keV
    max: EnergyType = 7 * u.keV


class ChandraConfig(BaseConfig):
    name: str = "my-analysis"
    sub_name: str = "my-config"
    obs_ids: List[int] = [1093]
    obs_id_ref: int = 1093
    roi: ROIConfig = ROIConfig()
    energy_range: EnergyRangeConfig = EnergyRangeConfig()
    ciao: CiaoToolsConfig = CiaoToolsConfig()

    def __str__(self):
        """Display settings in pretty YAML format."""
        info = self.__class__.__name__ + "\n\n\t"
        data = self.to_yaml()
        data = data.replace("\n", "\n\t")
        info += data
        return info.expandtabs(tabsize=2)

    @classmethod
    def read(cls, path):
        """Reads from YAML file."""
        log.info(f"Reading {path}")
        config = read_yaml(path)
        return ChandraConfig(**config)

    @classmethod
    def from_yaml(cls, config_str):
        """Create from YAML string."""
        settings = yaml.safe_load(config_str)
        return ChandraConfig(**settings)

    def write(self, path, overwrite=False):
        """Write to YAML file."""
        path = make_path(path)

        if path.exists() and not overwrite:
            raise IOError(f"File exists already: {path}")

        path.write_text(self.to_yaml())

    def to_yaml(self):
        """Convert to YAML string."""
        # Here using `dict()` instead of `json()` would be more natural.
        # We should change this once pydantic adds support for custom encoders
        # to `dict()`. See https://github.com/samuelcolvin/pydantic/issues/1043
        config = json.loads(self.json())
        return yaml.dump(
            config, sort_keys=False, indent=4, width=80, default_flow_style=False
        )

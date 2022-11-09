import json
import logging
from typing import Dict, List

import numpy as np
import yaml
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time
from astropy.units import Quantity
from ciao_contrib import runtool
from gammapy.analysis.config import AngleType, EnergyType, FrameEnum, GammapyBaseConfig
from gammapy.utils.scripts import make_path, read_yaml
from pydantic import BaseModel, create_model
from regions import CircleSkyRegion, RectangleSkyRegion

log = logging.getLogger(__name__)


CIAO_TOOLS_TYPES = {"f": str, "i": int, "s": str, "b": bool, "r": float}

CIAO_TOOLS_REQUIRED = {
    "dmcopy": {
        "infile": "{file_index.filename_repro_evt2_reprojected}",
        "outfile": "{file_index.filename_counts}",
    },
    "simulate_psf": {
        "infile": "{file_index.filename_repro_evt2_reprojected}",
        "outroot": "{{file_index.paths_psf_marx[{irf_label}]}}",
        "ra": np.nan,
        "dec": np.nan,
        "spectrumfile": "{{file_index.filenames_spectra[{irf_label}]}}",
    },
    "chandra_repro": {
        "indir": "{file_index.path_obs_id}",
        "outdir": "{file_index.path_repro}",
    },
    "specextract": {
        "infile": "{file_index.filename_repro_evt2_reprojected}",
        "outroot": "{{file_index.paths_spectra_pha[{irf_label}]}}/{irf_label}",
    },
    "reproject_events": {
        "infile": "{file_index.filename_repro_evt2}",
        "outfile": "{file_index.filename_repro_evt2_reprojected}",
        "match": "{file_index_ref.filename_repro_evt2}",
    },
}

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


def to_ciao_name(name):
    """Convert parameter name to ciao name"""
    return name.replace("_", "-")


class CiaoBaseConfig(BaseConfig):
    """Ciao tools base config"""

    def to_ciao(self, file_index, file_index_ref=None, irf_label=None):
        """Convert to ciao config dict"""
        kwargs = self.dict()

        for name, value in CIAO_TOOLS_REQUIRED[self._tool_name].items():
            ciao_name = to_ciao_name(name)

            if not isinstance(value, str):
                continue

            if irf_label and "irf_label" in value:
                value = value.format(irf_label=irf_label)

            kwargs[ciao_name] = value.format(
                file_index=file_index,
                file_index_ref=file_index_ref,
            )

        return kwargs


def create_ciao_config(toolname, model_name):
    """Create config class"""
    par_info = runtool.parinfo[toolname]

    parameters = {"_tool_name": toolname}

    for par in par_info["opt"]:
        parameters[par.name] = (CIAO_TOOLS_TYPES[par.type], par.default)

    model = create_model(model_name, __base__=CiaoBaseConfig, **parameters)

    return model


DMCopyConfig = create_ciao_config("dmcopy", "DMCopyConfig")
ChandraReproConfig = create_ciao_config("chandra_repro", "ChandraReproConfig")
ReprojectEventsConfig = create_ciao_config("reproject_events", "ReprojectEventsConfig")
SimulatePSFConfig = create_ciao_config("simulate_psf", "SimulatePSFConfig")
SpecExtractConfig = create_ciao_config("specextract", "SpecExtractConfig")


class CiaoToolsConfig(BaseConfig):
    dmcopy: DMCopyConfig = DMCopyConfig()
    chandra_repro: ChandraReproConfig = ChandraReproConfig()
    reproject_events: ReprojectEventsConfig = ReprojectEventsConfig()
    simulate_psf: SimulatePSFConfig = SimulatePSFConfig()
    specextract: SpecExtractConfig = SpecExtractConfig()


class EnergyRangeConfig(BaseConfig):
    min: EnergyType = 0.5 * u.keV
    max: EnergyType = 7 * u.keV

    def to_ciao(self):
        """dmcopy argument"""
        return f"energy={self.min.to_value('eV')}:{self.max.to_value('eV')}"


class SkyCoordConfig(BaseConfig):
    frame: FrameEnum = FrameEnum.icrs
    lon: AngleType = Angle("06h35m46.5079301472s")
    lat: AngleType = Angle("-75d16m16.816418256s")

    @property
    def sky_coord(self):
        """SkyCoord"""
        return SkyCoord(self.lon, self.lat, frame=self.frame)


class PerSourceSimulatePSFConfig(SimulatePSFConfig):
    class Config:
        fields = {
            "pileup": {"include": True},
            "readout_streak": {"include": True},
            "extended": {"include": True},
            "minsize": {"include": True},
        }


class PerSourceSpectrumConfig(SpecExtractConfig):
    center: SkyCoordConfig = SkyCoordConfig()
    radius: AngleType = Angle(3 * u.arcsec)
    energy_range: EnergyRangeConfig = EnergyRangeConfig()
    energy_groups: int = 5
    energy_step: float = 0.01

    @property
    def region(self):
        """Spectral extraction region"""
        return CircleSkyRegion(center=self.center.sky_coord, radius=self.radius)

    def to_ciao(self, file_index, file_index_ref=None, irf_label=None):
        """Spectrum extract region to ciao string"""
        kwargs = super().to_ciao(
            file_index=file_index, file_index_ref=file_index_ref, irf_label=irf_label
        )
        region_pix = self.region.to_pixel(wcs=file_index.wcs)
        x, y = region_pix.center.x, region_pix.center.y
        kwargs["infile"] += f"[sky=circle({x},{y},{region_pix.radius})]"
        return kwargs


class IRFConfig(BaseConfig):
    spectrum: PerSourceSpectrumConfig = PerSourceSpectrumConfig()
    psf: PerSourceSimulatePSFConfig = PerSourceSimulatePSFConfig()

    @property
    def psf_config_update(self):
        """Simulate psf config update"""
        config_psf = self.psf.dict()
        center = self.center.sky_coord
        config_psf["ra"] = center.icrs.ra.deg
        config_psf["dec"] = center.icrs.dec.deg
        return config_psf

    @property
    def spec_extract_config_update(self):
        """Specextract config update"""
        data = {}
        energy = self.energy_range
        data[
            "energy"
        ] = f"{energy.min.to_value('keV')}:{energy.max.to_value('keV')}:{self.energy_step}"
        return data

    @property
    def ciao(self):
        """Simulate PSF config"""
        config = CiaoToolsConfig()
        config.simulate_psf = config.simulate_psf.copy(update=self.psf_config_update)
        config.specextract = config.specextract.copy(
            update=self.spec_extract_config_update
        )
        return config


class ROIConfig(BaseConfig):
    center: SkyCoordConfig = SkyCoordConfig()
    width: AngleType = Angle("5 arcsec")
    bin_size: float = 1.0

    @property
    def region(self):
        """ROI region"""
        region = RectangleSkyRegion(
            center=self.center.sky_coord, width=self.width, height=self.width
        )
        return region

    def to_ciao(self, wcs):
        """dmcopy argument"""
        bbox = self.region.to_pixel(wcs).bounding_box
        return f"bin x={bbox.ixmin}:{bbox.ixmax}:0.5, y={bbox.iymin}:{bbox.iymax}:{self.bin_size}"


class ChandraConfig(BaseConfig):
    name: str = "my-analysis"
    sub_name: str = "my-config"
    obs_ids: List[int] = [62558]
    obs_id_ref: int = 62558
    roi: ROIConfig = ROIConfig()
    energy_range: EnergyRangeConfig = EnergyRangeConfig()
    irfs: Dict[str, IRFConfig] = {"pks-0637": IRFConfig()}
    ciao: CiaoToolsConfig = CiaoToolsConfig()

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

    def __str__(self):
        """Display settings in pretty YAML format."""
        info = self.__class__.__name__ + "\n\n\t"
        data = self.to_yaml()
        data = data.replace("\n", "\n\t")
        info += data
        return info.expandtabs(tabsize=2)

    def to_yaml(self):
        """Convert to YAML string."""
        # Here using `dict()` instead of `json()` would be more natural.
        # We should change this once pydantic adds support for custom encoders
        # to `dict()`. See https://github.com/samuelcolvin/pydantic/issues/1043

        data = self.json()

        config = json.loads(data)
        return yaml.dump(
            config, sort_keys=False, indent=4, width=80, default_flow_style=False
        )

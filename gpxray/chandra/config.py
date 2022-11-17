import json
import logging
import os
from pathlib import Path
from typing import ClassVar, Dict, List

import yaml
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time
from astropy.units import Quantity
from ciao_contrib import runtool
from gammapy.analysis.config import AngleType, EnergyType, FrameEnum
from gammapy.utils.scripts import make_path, read_yaml
from pydantic import BaseModel, create_model
from regions import CircleSkyRegion, RectangleSkyRegion

log = logging.getLogger(__name__)


MARX_ROOT = os.environ.get("CONDA_PREFIX", "${MARX_ROOT}")


CIAO_TOOLS_TYPES = {"f": str, "i": int, "s": str, "b": bool, "r": float}

CIAO_TOOLS_REQUIRED = {
    "dmcopy": {
        "infile": "{file_index.filename_repro_evt2_reprojected}",
        "outfile": "{file_index.filename_counts}",
    },
    "simulate_psf": {
        "infile": "{file_index.filename_repro_evt2_reprojected}",
        "outroot": "{{file_index.paths_psf_marx[{irf_label}]}}",
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
    "asphist": {
        "infile": "{file_index.filenames_repro_asol}",
        "outfile": "{file_index.filename_repro_asp_hist}",
        "evtfile": "{file_index.filename_repro_evt2_reprojected}",
    },
    "mkinstmap": {
        "outfile": "{file_index.filename_repro_inst_map}",
        "spectrumfile": "{{file_index.filenames_spectra[{irf_label}]}}",
        "obsfile": "{file_index.filename_repro_evt2_reprojected}[EVENTS]",
        "detsubsys": "ACIS-1",  # TODO: this should include all CCDs?
    },
    "mkexpmap": {
        "asphistfile": "{file_index.filename_repro_asp_hist}",
        "outfile": "{file_index.filename_exposure}",
        "instmapfile": "{file_index.filename_repro_inst_map}",
    },
}


class PathType(Path):
    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, v):
        return Path(v)


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

    _tool_name: ClassVar

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


def create_ciao_config(tool_name, model_name):
    """Create config class"""
    par_opts = runtool.parinfo[tool_name]["opt"]

    if tool_name == "simulate_psf":
        pars = [
            par
            for par in runtool.parinfo[tool_name]["req"]
            if par.name in ["ra", "dec"]
        ]
        par_opts.extend(pars)

    parameters = {}

    for par in par_opts:
        parameters[par.name] = (CIAO_TOOLS_TYPES[par.type], par.default)

    model = create_model(model_name, __base__=CiaoBaseConfig, **parameters)
    model._tool_name = tool_name

    return model


DMCopyConfig = create_ciao_config("dmcopy", "DMCopyConfig")
ChandraReproConfig = create_ciao_config("chandra_repro", "ChandraReproConfig")
ReprojectEventsConfig = create_ciao_config("reproject_events", "ReprojectEventsConfig")
SimulatePSFConfig = create_ciao_config("simulate_psf", "SimulatePSFConfig")
SpecExtractConfig = create_ciao_config("specextract", "SpecExtractConfig")
AspHistConfig = create_ciao_config("asphist", "AspHistConfig")
MkInstMapConfig = create_ciao_config("mkinstmap", "MkInstMapConfig")
MkExpMapConfig = create_ciao_config("mkexpmap", "MkExpMapConfig")


class CiaoToolsConfig(BaseConfig):
    dmcopy: DMCopyConfig = DMCopyConfig()
    chandra_repro: ChandraReproConfig = ChandraReproConfig()
    reproject_events: ReprojectEventsConfig = ReprojectEventsConfig()
    simulate_psf: SimulatePSFConfig = SimulatePSFConfig(marx_root=MARX_ROOT)
    specextract: SpecExtractConfig = SpecExtractConfig()
    asphist: AspHistConfig = AspHistConfig()
    mkinstmap: MkInstMapConfig = MkInstMapConfig()
    mkexpmap: MkExpMapConfig = MkExpMapConfig(normalize=False)


class EnergyRangeConfig(BaseConfig):
    min: EnergyType = 0.5 * u.keV
    max: EnergyType = 7 * u.keV


class SkyCoordConfig(BaseConfig):
    frame: FrameEnum = FrameEnum.icrs
    lon: AngleType = Angle("06h35m46.5079301472s")
    lat: AngleType = Angle("-75d16m16.816418256s")

    @property
    def sky_coord(self):
        """SkyCoord"""
        return SkyCoord(self.lon, self.lat, frame=self.frame)


class ROIConfig(DMCopyConfig):
    center: SkyCoordConfig = SkyCoordConfig()
    width: AngleType = Angle("5 arcsec")
    bin_size: float = 1.0
    energy_range: EnergyRangeConfig = EnergyRangeConfig()

    class Config:
        fields = {
            "center": {"include": True},
            "width": {"include": True},
            "bin_size": {"include": True},
            "energy_range": {"include": True},
        }

    @property
    def region(self):
        """ROI region"""
        region = RectangleSkyRegion(
            center=self.center.sky_coord, width=self.width, height=self.width
        )
        return region

    def region_shapely(self, wcs):
        """ROI shape"""
        from shapely.geometry import Polygon

        region_pix = self.region.to_pixel(wcs=wcs)
        return Polygon(region_pix.corners)

    def intersects_fov(self, file_index):
        """ROI intersects FoV"""
        intersects = {}

        roi_shape = self.region_shapely(wcs=file_index.wcs)

        for name, shape in file_index.fov_regions_shapely.items():

            intersects[name] = shape.intersects(roi_shape)

        return intersects

    def to_ciao(self, file_index, file_index_ref=None, irf_label=None):
        """dmcopy argument"""
        config = CiaoToolsConfig().dmcopy.copy()
        kwargs = config.to_ciao(
            file_index=file_index, file_index_ref=file_index_ref, irf_label=irf_label
        )

        bbox = self.region.to_pixel(wcs=file_index.wcs).bounding_box
        spatial = (
            f"bin x={bbox.ixmin}:{bbox.ixmax}:0.5, "
            f"y={bbox.iymin}:{bbox.iymax}:{self.bin_size}"
        )

        energy = self.energy_range
        spectral = f"energy={energy.min.to_value('eV')}:{energy.max.to_value('eV')}"

        selection = f"[EVENTS][{spatial}]"
        selection += f"[{spectral}]"
        kwargs["infile"] += selection
        return kwargs


class PerSourceSimulatePSFConfig(SimulatePSFConfig):
    class Config:
        fields = {
            "pileup": {"include": True},
            "readout_streak": {"include": True},
            "extended": {"include": True},
            "minsize": {"include": True},
        }

    def to_ciao(self, file_index, file_index_ref=None, irf_label=None):
        """Spectrum extract region to ciao config"""
        config = CiaoToolsConfig().simulate_psf.copy()
        kwargs = config.to_ciao(
            file_index=file_index, file_index_ref=file_index_ref, irf_label=irf_label
        )

        # Those are non visible linked parameters
        kwargs["ra"] = self.ra
        kwargs["dec"] = self.dec
        kwargs["bin_size"] = self.bins_size
        kwargs.update(self.dict())
        return kwargs


class PerSourceSpecExtractConfig(SpecExtractConfig):
    center: SkyCoordConfig = SkyCoordConfig()
    radius: AngleType = Angle(3 * u.arcsec)
    energy_range: EnergyRangeConfig = EnergyRangeConfig()
    energy_groups: int = 5
    energy_step: float = 0.01

    class Config:
        fields = {
            "center": {"include": True},
            "radius": {"include": True},
            "energy_range": {"include": True},
            "energy_groups": {"include": True},
            "energy_step": {"include": True},
        }

    @property
    def region(self):
        """Spectral extraction region"""
        return CircleSkyRegion(center=self.center.sky_coord, radius=self.radius)

    def to_ciao(self, file_index, file_index_ref=None, irf_label=None):
        """Spectrum extract region to ciao config"""
        config = CiaoToolsConfig().specextract.copy()
        kwargs = config.to_ciao(
            file_index=file_index, file_index_ref=file_index_ref, irf_label=irf_label
        )
        region_pix = self.region.to_pixel(wcs=file_index.wcs)
        x, y = region_pix.center.x, region_pix.center.y
        kwargs["infile"] += f"[sky=circle({x},{y},{region_pix.radius})]"

        energy = self.energy_range
        selection = f"{energy.min.to_value('keV')}:{energy.max.to_value('keV')}:{self.energy_step}"
        kwargs["energy"] = selection
        return kwargs


def region_to_ciao_str(region, wcs, bin_size):
    """Convert Astropy region to ciao string"""
    bbox = region.to_pixel(wcs=wcs).bounding_box

    nx = int(bbox.shape[1] / bin_size)
    ny = int(bbox.shape[0] / bin_size)

    return f"{bbox.ixmin}:{bbox.ixmax}:#{nx},{bbox.iymin}:{bbox.iymax}:#{ny}"


class PerSourceMkInstMapConfig(MkInstMapConfig):
    roi = ROIConfig()

    def to_ciao(self, file_index, file_index_ref=None, irf_label=None):
        """To ciao config"""
        config = CiaoToolsConfig().mkinstmap

        kwargs = config.to_ciao(
            file_index=file_index, file_index_ref=file_index_ref, irf_label=irf_label
        )
        kwargs["pixelgrid"] = region_to_ciao_str(
            region=self.roi.region, wcs=file_index.wcs, bin_size=self.roi.bin_size
        )

        intersections = self.roi.intersects_fov(file_index=file_index)

        for key, value in intersections.items():
            if value:
                idx = key
                break
        else:
            raise ValueError("No overlap betwen ROI and CCD FoV")

        kwargs["detsubsys"] = f"ACIS-{idx}"
        return kwargs


class PerSourceMkExpMapConfig(MkExpMapConfig):
    roi = ROIConfig()

    def to_ciao(self, file_index, file_index_ref=None, irf_label=None):
        """To ciao config"""
        config = CiaoToolsConfig().mkexpmap

        kwargs = config.to_ciao(
            file_index=file_index, file_index_ref=file_index_ref, irf_label=irf_label
        )

        kwargs["xygrid"] = region_to_ciao_str(
            region=self.roi.region, wcs=file_index.wcs, bin_size=self.roi.bin_size
        )
        return kwargs


class IRFConfig(BaseConfig):
    spectrum: PerSourceSpecExtractConfig = PerSourceSpecExtractConfig()
    psf: PerSourceSimulatePSFConfig = PerSourceSimulatePSFConfig()
    aeff: PerSourceMkInstMapConfig = PerSourceMkInstMapConfig()
    exposure: PerSourceMkExpMapConfig = PerSourceMkExpMapConfig()

    class Config:
        fields = {
            "spectrum": {"include": True},
            "psf": {"include": True},
        }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        center = self.spectrum.center.sky_coord
        self.psf.ra = center.icrs.ra.deg
        self.psf.dec = center.icrs.dec.deg


class ChandraConfig(BaseConfig):
    name: str = "my-analysis"
    sub_name: str = "my-config"
    path_data: PathType = Path("./data")
    obs_ids: List[int] = [62558]
    obs_id_ref: int = 62558
    roi: ROIConfig = ROIConfig()
    irfs: Dict[str, IRFConfig] = {"pks-0637": IRFConfig()}
    ciao: CiaoToolsConfig = CiaoToolsConfig()

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        for config in self.irfs.values():
            config.psf = self.ciao.simulate_psf.copy(update=config.psf.dict())
            config.psf.binsize = self.roi.bin_size
            config.aeff.roi = self.roi
            config.exposure.roi = self.roi

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

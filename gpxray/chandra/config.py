import json
from typing import List

import yaml
from astropy import units as u
from astropy.coordinates import Angle
from gammapy.analysis.config import (
    AngleType,
    EnergyType,
    GammapyBaseConfig,
    SkyCoordConfig,
)
from gammapy.utils.scripts import make_path, read_yaml


class ROIConfig(GammapyBaseConfig):
    center: SkyCoordConfig = SkyCoordConfig()
    width: AngleType = Angle("3 arcsec")


class EnergyRangeConfig(GammapyBaseConfig):
    min: EnergyType = 0.5 * u.keV
    max: EnergyType = 7 * u.keV


class ChandraConfig(GammapyBaseConfig):
    name: str = "my-analysis"
    sub_name: str = "my-config"
    obs_ids: List[int] = [1, 2, 3]
    obs_id_ref: int = 1
    roi: ROIConfig = ROIConfig()
    energy_range: EnergyRangeConfig = EnergyRangeConfig()

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
            config, sort_keys=False, indent=4, width=80, default_flow_style=None
        )

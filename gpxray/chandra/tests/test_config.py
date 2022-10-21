from gpxray.chandra.config import ChandraConfig

CONFIG_STR = """
path: .
name: my-analysis
sub_name: my-config
obs_ids:
- 1
- 2
- 3
obs_id_ref: 1
roi:
    center:
        frame: icrs
        lon: 0.0 deg
        lat: 0.0 deg
    width: 3.0 arcsec
energy_range:
    min: 0.5 keV
    max: 7.0 keV
ciao:
    dmcopy:
        infile: ""
        outfile: ""
"""


def test_chandra_config():
    config = ChandraConfig.from_yaml(CONFIG_STR)

    assert str(config.path) == "."
    assert config.roi.center.frame == "icrs"

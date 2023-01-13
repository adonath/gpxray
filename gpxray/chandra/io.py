from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.utils import lazyproperty
from astropy.wcs import WCS
from gammapy.data import EventList
from gammapy.maps import MapAxis, RegionGeom, RegionNDMap
from gammapy.maps.utils import edges_from_lo_hi
from regions import PixCoord, PolygonPixelRegion, RegionVisual

__all__ = [
    "read_spectrum_chart",
    "read_event_list_chandra",
    "ChandraFileIndex",
]


def wcs_from_header_chandra(header, x_col=11):
    """Create WCS from event file header

    Parameters
    ----------
    header : `~astropy.io.fits.Header`
        FITS header

    Returns
    -------
    wcs : `~astropy.wcs.WCS`
        WCS object
    """
    y_col = x_col + 1
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [header[f"TCRPX{x_col}"], header[f"TCRPX{y_col}"]]
    wcs.wcs.cdelt = [header[f"TCDLT{x_col}"], header[f"TCDLT{y_col}"]]
    wcs.wcs.crval = [header[f"TCRVL{x_col}"], header[f"TCRVL{y_col}"]]
    wcs.wcs.ctype = [header[f"TCTYP{x_col}"], header[f"TCTYP{y_col}"]]
    return wcs


def read_spectrum_chart(filename):
    """Read CHART spectrum"""
    data = Table.read(filename, format="ascii")

    edges = edges_from_lo_hi(data["col1"], data["col2"])
    axis = MapAxis.from_edges(edges=edges, name="energy", unit="keV", interp="log")

    geom = RegionGeom.create(region=None, axes=[axis])
    spectrum = RegionNDMap.from_geom(geom=geom)
    spectrum.data = data["col3"]
    return spectrum


def read_event_list_chandra(filename, hdu="EVENTS"):
    """Read event list"""
    hdu = fits.open(filename)[hdu]
    table = Table.read(hdu)

    wcs = wcs_from_header_chandra(header=hdu.header)

    for colname in table.colnames:
        table.rename_column(colname, colname.upper())

    ra, dec = wcs.wcs_pix2world(table["X"], table["Y"], 1)
    table["RA"] = ra * u.deg
    table["DEC"] = dec * u.deg

    mjd = table.meta["MJDREF"]
    mjd_int = np.floor(mjd).astype(np.int64)
    table.meta["MJDREFI"] = mjd_int
    table.meta["MJDREFF"] = mjd - mjd_int
    table.meta["TIMESYS"] = "tt"  # TODO: not sure tt is correct here...
    return EventList(table=table)


class ChandraFileIndex:
    """File index class"""

    def __init__(self, obs_id, path=".", path_output="my-config", irf_names=None):
        self.obs_id = obs_id
        self._path = Path(path)
        self._path_output = Path(path_output)
        self._irf_names = irf_names

    @property
    def t_start(self):
        """Start time"""
        header = fits.getheader(self.filename_repro_asol1)
        t_start = header["TSTART"]
        return float(t_start)

    @property
    def t_stop(self):
        """Start time"""
        header = fits.getheader(self.filename_repro_asol1)
        t_stop = header["TSTOP"]
        return float(t_stop)

    @property
    def limit(self):
        """Limit"""
        return self.t_stop - self.t_start

    @property
    def ra_pnt(self):
        """RA pointing"""
        header = fits.getheader(self.filename_repro_evt2_reprojected, "EVENTS")
        return float(header["RA_PNT"])

    @property
    def dec_pnt(self):
        """DEC pointing"""
        header = fits.getheader(self.filename_repro_evt2_reprojected, "EVENTS")
        return float(header["DEC_PNT"])

    @property
    def roll_pnt(self):
        """ROL pointing"""
        header = fits.getheader(self.filename_repro_evt2_reprojected, "EVENTS")
        return float(header["ROLL_PNT"])

    @lazyproperty
    def pointing(self):
        """Pointing position"""
        return SkyCoord(self.ra_pnt, self.dec_pnt, unit="deg", frame="icrs")

    @property
    def irf_names(self):
        """IRF names"""
        return self._irf_names

    @lazyproperty
    def wcs(self):
        """Get WCS from reprojected events file"""
        header = fits.getheader(self.filename_repro_evt2_reprojected, "EVENTS")
        return wcs_from_header_chandra(header=header)

    @property
    def ccds(self):
        """CCDs which were part of the obs"""
        detnam = self.index_table.meta["DETNAM"]

        prefix, indices = detnam.split("-")

        ccds = [f"{prefix}-{idx}" for idx in list(indices)]
        return ccds

    @property
    def fov_regions(self):
        """FoV regions"""
        filename = self.index_filenames["FOV"]
        table = Table.read(filename, hdu="FOV")

        fovs = {}

        for row, color in zip(table, plt.color_sequences["tab10"]):
            coords = PixCoord(x=row["X"], y=row["Y"])
            region = PolygonPixelRegion(
                vertices=coords, visual=RegionVisual({"color": color})
            )
            fovs[row["CCD_ID"]] = region

        return fovs

    @property
    def fov_regions_shapely(self):
        """FoV shapes"""
        from shapely.geometry import Polygon

        shapes = {}

        for name, region in self.fov_regions.items():
            x, y = region.vertices.xy
            shape = Polygon(shell=zip(x, y))
            shapes[name] = shape

        return shapes

    @property
    def path(self):
        """Data location path"""
        return self._path

    @property
    def path_data(self):
        """Data location path"""
        return self._path

    @property
    def path_obs_id(self):
        """Data location path"""
        return self.path_data / f"{self.obs_id}"

    @property
    def source_name(self):
        """Source name (`str`)"""
        return self.index_table.meta["OBJECT"].strip()

    @property
    def path_repro(self):
        """Reprocessed data path"""
        return self.path_obs_id / "repro"

    @property
    def paths_psf_marx(self):
        """PSF data path marx"""
        paths = {}

        for name in self.irf_names:
            path = self.path_output / f"psf-marx-{name}"
            path.mkdir(parents=True, exist_ok=True)
            paths[name] = path

        return paths

    @property
    def paths_psf_saotrace(self):
        """PSF data path saotrace"""
        paths = {}

        for name in self.irf_names:
            path = self.path_output / f"psf-sao-trace-{name}"
            path.mkdir(parents=True, exist_ok=True)
            paths[name] = path

        return paths

    @property
    def paths_spectra_pha(self):
        """Spectrum data path"""
        paths = {}

        for name in self.irf_names:
            path = self.path_output / f"spectrum-{name}"
            path.mkdir(parents=True, exist_ok=True)
            paths[name] = path

        return paths

    @property
    def filename_repro_evt2(self):
        """Reprocessed data path"""
        return self.path_repro / f"acisf{self.obs_id:05d}_repro_evt2.fits"

    @property
    def filename_repro_evt2_reprojected(self):
        """Reprocessed data path"""
        return self.path_repro / f"acisf{self.obs_id:05d}_repro_evt2_reprojected.fits"

    @property
    def filename_repro_asol1(self):
        """Aspect solution file"""
        filenames = sorted(self.path_repro.glob("*asol1.fits"))
        return filenames[0]

    @property
    def filenames_repro_asol(self):
        """Aspect solution file"""
        filenames = self.path_repro.glob("*asol1.fits")
        return ",".join([str(_) for _ in filenames])

    @property
    def filename_repro_asp_hist(self):
        """Aspect histogram"""
        return self.path_repro / f"acisf{self.obs_id:05d}_asp_hist.fits"

    @property
    def filename_repro_inst_map(self):
        """Instrument map"""
        return self.path_repro / f"acisf{self.obs_id:05d}_inst_map.fits"

    @property
    def index_table(self):
        """Index table (`astropy.table.Table`)"""
        index_table = Table.read(self.path_obs_id / "oif.fits")
        index_table.add_index("MEMBER_CONTENT")
        return index_table

    @lazyproperty
    def index_filenames(self):
        """Index filenames (`dict`)"""
        filenames = {}

        keys = [_.strip() for _ in self.index_table["MEMBER_NAME"]]

        for key, filename in zip(keys, self.index_table["MEMBER_LOCATION"]):
            filename = filename.strip()

            if key not in ["FOV"]:
                filename = filename.replace(".fits", ".fits.gz")

            filenames[key] = self.path_obs_id / filename

        return filenames

    @property
    def path_output(self):
        """Output path"""
        path = self._path_output / f"{self.obs_id}"
        path.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def filename_counts(self):
        """Filename counts"""
        return self.path_output / "counts.fits"

    @property
    def filenames_psf(self):
        """Filenames psf"""
        filenames = {}

        for name in self.irf_names:
            filename = self.path_output / f"psf-{name}.fits"
            filenames[name] = filename

        return filenames

    @property
    def filenames_spectra(self):
        """Filename spectra"""
        filenames = {}

        for name in self.irf_names:
            filename = (
                self.path_output / f"spectrum-{name}" / f"source-flux-chart-{name}.dat"
            )
            filenames[name] = filename

        return filenames

    @property
    def filenames_spectra_png(self):
        """Filename spectra"""
        filenames = {}

        for name in self.irf_names:
            filename = self.path_output / f"spectrum-{name}" / f"spectrum-{name}.png"
            filenames[name] = filename

        return filenames

    @property
    def filename_exposure(self):
        """Filename counts"""
        return self.path_output / "exposure.fits"

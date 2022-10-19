from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from gammapy.data import EventList


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


def read_event_list_chandra(filename):
    """Read event list"""
    table = Table.read(filename)
    header = fits.getheader(filename, "EVENTS")
    wcs = wcs_from_header_chandra(header=header)

    ra, dec = wcs.wcs_pix2world(table["x"], table["y"], 1)
    table["RA"] = ra * u.deg
    table["DEC"] = dec * u.deg
    table.rename_column("energy", "ENERGY")
    table.rename_column("time", "TIME")

    mjd = table.meta["MJDREF"]
    mjd_int = np.floor(mjd).astype(np.int64)
    table.meta["MJDREFI"] = mjd_int
    table.meta["MJDREFF"] = mjd - mjd_int
    table.meta["TIMESYS"] = "tt"  # TODO: not sure tt is correct here...
    return EventList(table=table)


class ChandraFileIndex:
    """File index class"""

    def __init__(self, obs_id, path="."):
        self.obs_id = obs_id
        self._path = Path(path)

    @property
    def path(self):
        """Data location path"""
        return self._path

    @property
    def path_data(self):
        """Data location path"""
        return self._path / "data"

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
        return self.path_obsid / "repro"

    @property
    def path_psf(self):
        """PSF data path"""
        path = self.path_obs_id / "psf"
        path.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def filename_repro_evt2(self):
        """Reprocessed data path"""
        return self.path_repro / f"acisf{self.obs_id:05d}_repro_evt2.fits"

    @property
    def filename_repro_evt2_reprojected(self):
        """Reprocessed data path"""
        return self.path_repro / f"acisf{self.obs_id:05d}_repro_evt2_reprojected.fits"

    @property
    def index_table(self):
        """Index table (`astropy.table.Table`)"""
        index_table = Table.read(self.path_obs_id / "oif.fits")
        index_table.add_index("MEMBER_CONTENT")
        return index_table

    def get_filename(self, member_content):
        """Get file name"""
        filename = self.index_table.loc[f"{member_content:32s}"]["MEMBER_LOCATION"]
        filename = filename.strip()
        filename = filename.replace(".fits", ".fits.gz")
        return self.path_base / filename

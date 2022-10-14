import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from gammapy.data import EventList


def wcs_from_header_chandra(header, x_col=11):
    """Create WCS from event file header"""
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

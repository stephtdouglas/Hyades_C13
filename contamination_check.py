from __future__ import division, print_function
import sys, os

import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
# from astropy.table import Table, hstack, join
from astropy.time import Time
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord, Distance
import astropy.units as u
import astroquery
# import matplotlib.pyplot as plt
# import matplotlib.ticker
# import matplotlib.cm as cm
# import matplotlib.colors as colors
# from matplotlib.font_manager import FontProperties
# import palettable

from hypra.utils import cat_io

def contamination_check():

    hdat, _, _, _ = cat_io.get_data("H")
    nh = len(hdat)
    hpos = SkyCoord(hdat["RA"],hdat["Dec"],unit=u.degree)
    hpos_offset = SkyCoord(hdat["RA"]+50,hdat["Dec"],unit=u.degree)

    match10 = np.zeros(nh,bool)
    match20 = np.zeros(nh,bool)

    offset_match10 = np.zeros(nh,bool)
    offset_match20 = np.zeros(nh,bool)

    for i in range(nh):
        print(i)
        pos = hpos[i]
        offpos = hpos_offset[i]

        # For comparison, look at the basic positions
        # (even though we know they all have matches)
        gquery = Gaia.cone_search(pos,radius=30*u.arcsec)
        gout = gquery.get_data()
        gdist = gout["dist"]*u.degree

        this_match10 = len(np.where((gdist<(10*u.arcsec)) &
                        (gout["phot_g_mean_mag"]<22))[0])

        this_match20 = len(np.where((gdist<(20*u.arcsec)) &
                        (gout["phot_g_mean_mag"]<22))[0])

        match10[i] = True if (this_match10>0) else False
        match20[i] = True if (this_match20>0) else False

        # Now the offset positions
        gquery = Gaia.cone_search(offpos,radius=30*u.arcsec)
        gout = gquery.get_data()
        gdist = gout["dist"]*u.degree

        this_match10 = len(np.where((gdist<(10*u.arcsec)) &
                        (gout["phot_g_mean_mag"]<22))[0])

        this_match20 = len(np.where((gdist<(20*u.arcsec)) &
                        (gout["phot_g_mean_mag"]<22))[0])

        offset_match10[i] = True if (this_match10>0) else False
        offset_match20[i] = True if (this_match20>0) else False


        # if i>10:
        #     break

    print(len(np.where(match10)[0])," with a Gaia match within 10\" ",
          "({0:.1f}\%)".format(100*len(np.where(match10)[0])/nh))
    print(len(np.where(match20)[0])," with a Gaia match within 20\" ",
          "({0:.1f}\%)".format(100*len(np.where(match20)[0])/nh))

    print(len(np.where(offset_match10)[0]),
          " offset positions with a Gaia match within 10\" ",
          "({0:.1f}\%)".format(100*len(np.where(offset_match10)[0])/nh))
    print(len(np.where(offset_match20)[0]),
          " offset positions with a Gaia match within 20\" ",
          "({0:.1f}\%)".format(100*len(np.where(offset_match20)[0])/nh))


if __name__=="__main__":
    contamination_check()

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
import matplotlib.pyplot as plt
# import matplotlib.ticker
# import matplotlib.cm as cm
# import matplotlib.colors as colors
# from matplotlib.font_manager import FontProperties
# import palettable
import astropy

from hypra.utils import cat_io
import plot_allsky

def contamination_check(#ra_off=0*u.degree,sep_off=0*u.degree,
                        pa_off=0*u.degree,sep_off=20*u.degree,
                        to_plot=False,to_output=True):

    if to_output:
        filebase = "Gaia_contam_pa{0:3.0f}_sep{1:3.0f}".format(
                    pa_off.value,sep_off.value)
        print(filebase)

    hdat, _, _, _ = cat_io.get_data("H")
    nh = len(hdat)
    hpos = SkyCoord(hdat["RA"],hdat["Dec"],unit=u.degree)

    # Below only works in Astropy 3.1
    hpos_offset = hpos.directional_offset_by(position_angle=pa_off,
                                             separation=sep_off)

    # match10 = np.zeros(nh,bool)
    # match20 = np.zeros(nh,bool)

    offset_match10 = np.zeros(nh,bool)
    offset_match20 = np.zeros(nh,bool)

    for i in range(nh):
        print(i)
        pos = hpos[i]
        offpos = hpos_offset[i]

        # # For comparison, look at the basic positions
        # # (even though we know they all have matches)
        # gquery = Gaia.cone_search(pos,radius=30*u.arcsec)
        # gout = gquery.get_data()
        # gdist = gout["dist"]*u.degree
        #
        # this_match10 = len(np.where((gdist<(10*u.arcsec)) &
        #                 (gout["phot_g_mean_mag"]<22))[0])
        #
        # this_match20 = len(np.where((gdist<(20*u.arcsec)) &
        #                 (gout["phot_g_mean_mag"]<22))[0])
        #
        # match10[i] = True if (this_match10>0) else False
        # match20[i] = True if (this_match20>0) else False

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


        if i>50:
            break

    # print(len(np.where(match10)[0])," with a Gaia match within 10\" ",
    #       "({0:.1f}\%)".format(100*len(np.where(match10)[0])/nh))
    # print(len(np.where(match20)[0])," with a Gaia match within 20\" ",
    #       "({0:.1f}\%)".format(100*len(np.where(match20)[0])/nh))

    print(len(np.where(offset_match10)[0]),
          " offset positions with a Gaia match within 10\" ",
          "({0:.1f}\%)".format(100*len(np.where(offset_match10)[0])/nh))
    print(len(np.where(offset_match20)[0]),
          " offset positions with a Gaia match within 20\" ",
          "({0:.1f}\%)".format(100*len(np.where(offset_match20)[0])/nh))

    if to_plot:
        ax = plot_allsky.plot_equatorial_base()
        plot_allsky.plot_hyades(ax)
        plot_offset_ra = hpos_offset.ra.wrap_at(180 * u.deg).radian
        plot_offset_dec = hpos_offset.dec.radian
        ax.plot(plot_offset_ra,plot_offset_dec,'o',
                mfc="Grey",mec="Grey",ms=3,label="Offset positions",
                zorder=20)
        ax.plot(plot_offset_ra[offset_match10],
                plot_offset_dec[offset_match10],'o',
                mfc="C4",mec="C4",ms=3,label="With Gaia Matches",
                zorder=21)
        plt.show()


if __name__=="__main__":
    print(astropy.__version__)
    contamination_check(to_plot=True)

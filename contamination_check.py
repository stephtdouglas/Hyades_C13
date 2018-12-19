from __future__ import division, print_function
import sys, os
import glob

import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
from astropy.table import Table #, hstack, join
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
                        pa_off=0*u.degree,sep_off=30*u.degree,
                        to_plot=False,to_output=True):

    if to_output:
        filebase = "Gaia_contam_pa{0:03.0f}_sep{1:03.0f}".format(
                    pa_off.value,sep_off.value)
        print(filebase)

    hdat, _, _, _ = cat_io.get_data("H")
    hpos = SkyCoord(hdat["RA"],hdat["Dec"],unit=u.degree)

    hcenter = SkyCoord("4:28:25.7 +16:42:45",unit=(u.hourangle,u.degree))
    core_radius = 18*u.degree
    sep = hpos.separation(hcenter)

    hidx0 = np.arange(len(hdat))
    hpos = hpos[sep<core_radius]
    nh = len(hpos)
    hidx = hidx0[sep<core_radius]

    # Below only works in Astropy 3.1
    hpos_offset = hpos.directional_offset_by(position_angle=pa_off,
                                             separation=sep_off)

    # match10 = np.zeros(nh,bool)
    # match20 = np.zeros(nh,bool)

    offset_match10 = np.zeros(nh,bool)
    nmatch10 = np.zeros(nh,int)
    offset_match20 = np.zeros(nh,bool)
    nmatch20 = np.zeros(nh,int)

    for i in range(nh):
        # print(i)
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

        nmatch10[i] = len(np.where((gdist<(10*u.arcsec)) &
                        (gout["phot_g_mean_mag"]<22))[0])

        nmatch20[i] = len(np.where((gdist<(20*u.arcsec)) &
                        (gout["phot_g_mean_mag"]<22))[0])

        offset_match10[i] = True if (nmatch10[i]>0) else False
        offset_match20[i] = True if (nmatch20[i]>0) else False


        # if i>1:
        #     break

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
        # ax.plot(hpos.ra.wrap_at(180 * u.deg).radian,
        #         hpos.dec.radian,'o',color='C1',ms=3)
        plot_offset_ra = hpos_offset.ra.wrap_at(180 * u.deg).radian
        plot_offset_dec = hpos_offset.dec.radian
        ax.plot(plot_offset_ra,plot_offset_dec,'o',
                mfc="none",mec="Grey",ms=3,label="Offset positions",
                zorder=20)
        ax.plot(plot_offset_ra[offset_match10],
                plot_offset_dec[offset_match10],'o',
                mfc="C3",mec="C3",ms=3,label="With Gaia Matches",
                zorder=21)
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
        ax.set_title("PA={0:.0f}, Sep={1:.0f}".format(
                     pa_off.value,sep_off.value))
        ax.legend(bbox_to_anchor=(0.25, 0), framealpha=1)
        if to_output:
            plt.savefig(filebase+".png",dpi=600  )
        else:
            plt.show()

    if to_output:
        out_tab = Table([hidx,hpos_offset.ra,hpos_offset.dec,
                         offset_match10,nmatch10,
                         offset_match20,nmatch20],
                         names=["IDX","RA","Dec","Match10","NMatch10",
                                "Match20","NMatch20"])

        at.write(out_tab,filebase+".csv",delimiter=",",
                 formats={"RA":"%.6f","Dec":"%.6f"})

def check_circle(sep_off=30*u.degree,to_plot=True,to_output=True):

    for pa in np.linspace(0,360,13):
        contamination_check(pa_off=pa*u.degree,sep_off=sep_off,
                            to_plot=to_plot,to_output=to_output)

        # if pa>=30:
        #     break

def summarize_contamination(base_path="."):

    files = glob.glob(os.path.join(base_path,"Gaia_contam_*.csv"))

    if len(files)<=1:
        print("No files found! Please check directory or run check_circle")
    else:
        fraction10 = np.zeros(len(files))
        fraction20 = np.zeros(len(files))
        for i,filename in enumerate(files):
            contam = at.read(filename)
            print(contam.dtype)
            # IDX,RA,Dec,Match10,NMatch10,Match20,NMatch20
            fraction10[i] = len(np.where(contam["Match10"]=="True")[0])/len(contam)
            fraction20[i] = len(np.where(contam["Match20"]=="True")[0])/len(contam)

        print("{0:.1f}\% - {1:.1f}\% with a match in 10 arcsec".format(
              min(fraction10)*100,max(fraction10)*100))

        print("{0:.1f}\% - {1:.1f}\% with a match in 20 arcsec".format(
              min(fraction20)*100,max(fraction20)*100))

if __name__=="__main__":
    # print(astropy.__version__)
    # contamination_check(to_plot=True)
    # check_circle()
    summarize_contamination()

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.font_manager import FontProperties
import palettable
import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
import astropy.table as table
import astropy.units as u
from astropy.coordinates import SkyCoord

from hypra.utils import cat_io

mpl.rcParams['lines.markeredgewidth'] = 1.5
mpl.rcParams['lines.markersize'] = 8
mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.titlesize'] = 24
mpl.rcParams['axes.labelsize'] = 22
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.width'] = 1
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['legend.fontsize'] = 16


pdat,_,_,_ = cat_io.get_data("P")
hdat,_,_,_ = cat_io.get_data("H")

p_pos = SkyCoord(pdat["RA"],pdat["DEC"],unit=u.deg,frame="icrs")
h_pos = SkyCoord(hdat["RA"],hdat["DEC"],unit=u.deg,frame="icrs")

ecliptic = SkyCoord(np.linspace(0,360,1000),np.zeros(1000),
                    unit=u.deg,frame="barycentrictrueecliptic")

galactic_plane = SkyCoord(np.linspace(0,360,1000),np.zeros(1000),
                          unit=u.deg,frame="galactic")

def plot_equatorial_base():

    plt.figure(figsize=(10,8))
    ax2 = plt.subplot(111, projection="aitoff")

    ax2.set_xticklabels(['14h','16h','18h','20h','22h','0h',
                        '2h','4h','6h','8h','10h'])
    ax2.grid(True)

    ecliptic_ra = ecliptic.icrs.ra.wrap_at(180 * u.deg).radian
    ecliptic_dec = ecliptic.icrs.dec.radian
    new_ecl = np.argsort(ecliptic_ra)
    ecliptic_ra = ecliptic_ra[new_ecl]
    ecliptic_dec = ecliptic_dec[new_ecl]

    ax2.plot(ecliptic_ra,ecliptic_dec,'-',color="C2",lw=2,zorder=-1)

    galactic_ra = galactic_plane.icrs.ra.wrap_at(180 * u.deg).radian
    galactic_dec = galactic_plane.icrs.dec.radian
    new_gal = np.argsort(galactic_ra)
    galactic_ra = galactic_ra[new_gal]
    galactic_dec = galactic_dec[new_gal]

    ax2.plot(galactic_ra,galactic_dec,'-',color="C4",lw=2,zorder=-1)

    return ax2

def plot_hyades(ax):

    ax.plot(h_pos.ra.wrap_at(180 * u.deg).radian,h_pos.dec.radian,
         'o',mfc="C1",mec="C1",ms=3,label="Hyades")


if __name__=="__main__":

    ax = plot_equatorial_base()
    plot_hyades(ax)
    plt.show()

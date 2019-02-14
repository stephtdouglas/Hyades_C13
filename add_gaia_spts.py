import os
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import astropy.io.ascii as at



# Read in Kraus & Hillenbrand for later use
kh = at.read(os.path.expanduser('~/HyPra/models/kraushillenbrand5.dat'))
ksub = np.array([-23,-19,-15,-12,-11,-9,-7,-6,-5,-4,-3,-2,-1])
# kh_rpmK0 = (kh['Mr'] - 0.035*(kh['Mr']-kh['Mi']) - 0.007 - kh['MK'])
# kh_rpmK = kh_rpmK0[ksub]
# kh_rpmK[0] += 0.1
kh_teff = kh["Teff"][ksub]
kh_spt = kh['SpT'][ksub]

def calc_logTeff(bprp):
    # From https://arxiv.org/pdf/1008.0815.pdf
    # ARGH it only applies for bp-rp<1.5
    return 3.999  - 0.654*bprp + 0.709*(bprp**2) - 0.316*(bprp**3)

jprae = at.read(os.path.expanduser("~/my_papers/hyadesk22/NotesJC/Table-Praesepe.txt"),delimiter=" ")
def calc_bprp_jason():
    bprp_order = np.argsort(jprae["bprp"])
    jbprp = jprae["bprp"][bprp_order]
    jteff = jprae["Teff"][bprp_order]
    calc_bprp = interp1d(jteff,jbprp,bounds_error=False,fill_value=np.nan)
    return calc_bprp

def add_gaia_spts(ax,xmin=None,xmax=None,texty=None):

    if (xmin is None) or (xmax is None):
        xlims = ax.get_xlim()
        xmin, xmax = xlims

    # Determine the proper text position
    if texty is None:
        ylims = ax.get_ylim()
        texty = ylims[1]*1.05

    calc_bprp = calc_bprp_jason()

    for this_teff, this_spt in zip(kh_teff, kh_spt):

        this_bprp = calc_bprp(this_teff)
        print(this_bprp)
        if np.isnan(this_bprp) or (this_bprp<xmin) or (this_bprp>xmax):
            continue
        else:
            ax.text(this_bprp,texty,this_spt,fontsize="large")

    # Remove the tick marks to reduce confusion
    ax.tick_params(which="both",top=False)

if __name__=="__main__":

    # bprp = np.linspace(0.75,4,100)
    # logTeff = calc_logTeff(bprp)
    # plt.plot(bprp,logTeff)
    # plt.show()

    plt.figure()
    ax = plt.subplot(111)
    ax.set_xlim(0.75,4)
    add_gaia_spts(ax)
    plt.show()

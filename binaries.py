
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import palettable
import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits

# from astropy.coordinates import SkyCoord  # High-level coordinates
# import astropy.units as u

from hypra.utils import cat_io, cat_match
from hypra.plot import color_mag

hdat,hobs,hobsnr,hobsr = cat_io.get_data("H")

# A few definitions
single_figure = (8,8)
double_figure = (13,6)
double_stacked = (6,13)
triple_stacked = (6,14)
quad_square = (13,13)

std_ms = 9
std_mew = 1.5

# Define colors, same as Praesepe paper
color_cycle = palettable.colorbrewer.qualitative.Dark2_8.mpl_colors
ms_color = color_cycle[3] #"#5d58a7" #
cand_color = color_cycle[2]
lit_color = color_cycle[1]
conf_color = "k"
background_color = plt.cm.Greys(0.4)


def plot_binaries(ax, x, y, plot_all=False, plot_confirmed=True, plot_candidate=True, plot_single=False):

    # All cluster members
    if plot_all:
        good = np.where((x>0) & (y>-9))[0]
        ax.plot(x[good], y[good], '.', ms=5, color=background_color,
                label="Cluster Member")

    # Confirmed binaries
    if plot_confirmed:
        binary = np.where((x>0) & (y>-9) & (hdat["BINARY"]==2))[0]
        ax.plot(x[binary], y[binary], '*', mfc=conf_color, ms=10,
                mec=conf_color, label="Confirmed Multiple")


    our_cand = np.zeros(len(hdat),int)
    for i,cite in enumerate(hdat["BINARY_CITE"]):
        if "douglas" in cite:
            our_cand[i] = 1
        else:
            continue
    # Candidate binaries
    if plot_candidate:
        binary = (x>0) & (y>-9) & (hdat["BINARY"]==1)
        ours = our_cand==1
        ax.plot(x[binary & ours], y[binary & ours], 'o', mfc="none", mew=1.5,
                mec=cand_color, label="Candidate Multiple (Douglas+2014)")
        other = our_cand==0
        nother = len(np.where(other & binary)[0])
        print(nother)
        if nother>0:
            ax.plot(x[binary & other], y[binary & other], 'D', mfc="none", mew=1.5,
                    mec=lit_color, label="Candidate Multiple (Literature)")

    ax.tick_params(width=1.5, length=6)
    ax.tick_params(which="minor",width=1.5, length=3)

def calc_residual(abs_r, dmod):
    fig = plt.figure()
    ax = plt.subplot(111)

    ms,color = color_mag.add_ms(ax,dmod,return_ms=True,line_color=ms_color)
    plt.close("all")

    model_r = np.ones_like(abs_r)*-9999.
    good_interp = np.where((hdat["RPRIME_K"]>min(color)) & (hdat["RPRIME_K"]<max(color)))[0]
    model_r[good_interp] = ms(hdat["RPRIME_K"][good_interp])
    residual = abs_r - model_r
    residual[model_r<-9998] = -9999

    return residual

def paper_plot(abs_r,dmod):

    residual = calc_residual(abs_r,dmod)
    fig = plt.figure(figsize=triple_stacked)

    #### TOP - CMD with binary sequence overlaid
    ax = plt.subplot(311)
    #ax = plt.subplot2grid((3,1),(0,0),rowspan=2)
    color_mag.setup_axes(ax,"(r'-K)",r"M$_{r'}$",[0,6.75,17,0.5])
    color_mag.add_spt(ax, 0.2)

    plot_binaries(ax, hdat["RPRIME_K"], abs_r, True,False,False)
    ax.set_xlabel("")
    ax.tick_params(labelbottom=False)
    ax.axvline(4,color="Grey",zorder=-111)
    ms,color = color_mag.add_ms(ax,dmod,return_ms=True,line_color=ms_color)

    l1 = ax.axhline(-99, color=ms_color, ls="-", lw=2, label="Main Sequence")
    l2 = ax.axhline(-99, color=ms_color, ls="-.", lw=2, label="Candidate Binary Threshold")
    l3 = ax.axhline(-99, color=ms_color, ls=":", lw=2, label="Binary Main Sequence")
    l4 = ax.axvline(-99,color="Grey",zorder=-111,label="Candidate Binary Red Limit")
    leg = ax.legend(handles=[l1,l3,l2,l4],loc=1,borderaxespad=0)
    leg.get_frame().set_alpha(1)



    #### MIDDLE - residuals with candidates identified
    ax = plt.subplot(312)
    #ax = plt.subplot2grid((3,1),(2,0),rowspan=1)
    color_mag.setup_axes(ax,"(r'-K)",r"M$_{r'}$ - M$_{r'}$(Main Sequence)",[0,6.75,0.9,-2.49])
    plt.subplots_adjust(hspace=0)
    plot_binaries(ax, hdat["RPRIME_K"], residual, True,False,False)
    cand = np.where((hdat["RPRIME_K"]>0) & (residual>-9) & (hdat["RPRIME_K"]<4) & (residual<(-0.75/2)))[0]
    ax.plot(hdat["RPRIME_K"][cand], residual[cand], 'o', mfc="none", mew=1.5,
            mec=cand_color, label="Candidate Multiple")
    ax.axvline(4,color="Grey",zorder=-111)
    ax.set_xlabel("")
    ax.tick_params(labelbottom=False)

    l1 = ax.axhline(0, color=ms_color, ls="-", lw=2, label="Main Sequence")
    l2 = ax.axhline(-0.75/2, color=ms_color, ls="-.", lw=2, label="Candidate Binary Threshold")
    l3 = ax.axhline(-0.75, color=ms_color, ls=":", lw=2, label="Binary Main Sequence")
    l4 = ax.axvline(4,color="Grey",zorder=-111,label="Candidate Binary Red Limit")
    #leg = ax.legend(handles=[l1,l3,l2,l4],loc=2,borderaxespad=0)
    #leg.get_frame().set_alpha(1)
    ax.set_ylim(1.1,-2.8)


    #### BOTTOM- residuals with confirmed identified
    ax = plt.subplot(313)
    #ax = plt.subplot2grid((3,1),(2,0),rowspan=1)
    color_mag.setup_axes(ax,"(r'-K)",r"M$_{r'}$ - M$_{r'}$(Main Sequence)",[0,6.75,0.9,-2.49])
    plt.subplots_adjust(hspace=0)
    plot_binaries(ax, hdat["RPRIME_K"], residual, True)
    print(len(np.where(residual>-1)[0]))
    leg = ax.legend(loc=1, numpoints=1,borderaxespad=0)
    leg.get_frame().set_alpha(1)

    l1 = ax.axhline(0, color=ms_color, ls="-", lw=2, label="Main Sequence")
    l2 = ax.axhline(-0.75/2, color=ms_color, ls="-.", lw=2, label="Candidate Binary Threshold")
    l3 = ax.axhline(-0.75, color=ms_color, ls=":", lw=2, label="Binary Main Sequence")
    l4 = ax.axvline(4,color="Grey",zorder=-111,label="Candidate Binary Red Limit")

    ax.set_ylim(1.1,-2.8)

    # plt.show()

    # plt.savefig("/home/stephanie/my_papers/praeK2/binary_cmd.eps",
    #             bbox_inches="tight")
    # plt.savefig("/home/stephanie/my_papers/praeK2/fig6.eps",
    #             bbox_inches="tight")
    # plt.savefig("/home/stephanie/my_papers/praeK2/fig6.pdf",
    #             bbox_inches="tight")
    # plt.savefig("/home/stephanie/Dropbox/plots_for_sharing/binary_cmd_praesepe.png",
    #             bbox_inches="tight"))

def gaia_binaries_plot():

    gaia = at.read("Gaia_Comb_Table.csv")

    gaia_plx = np.zeros(len(hdat))*np.nan
    good_plx = gaia["Plx"].mask==False
    gaia_plx[gaia["HYADES_IDX"][good_plx]] = gaia["Plx"][good_plx]

    gaia_dist = 1000/gaia_plx

    dmod = 5 * np.log10(gaia_dist) - 5
    abs_r = hdat["RPRIME"] - dmod
    abs_r[hdat["RPRIME"]<-9998] = -9999

    paper_plot(abs_r,dmod)

    plt.suptitle("Gaia parallax for all")
    plt.savefig("binary_cmd_gaiadist_oldbin.png",bbox_inches="tight")
    plt.close("all")

def old_binaries_plot():
    dmod = 5 * np.log10(hdat["DISTANCE"]) - 5
    abs_r = hdat["RPRIME"] - dmod
    abs_r[hdat["RPRIME"]<-9998] = -9999

    paper_plot(abs_r,dmod)
    plt.savefig("binary_cmd_olddist_oldbin.png",bbox_inches="tight")
    plt.close("all")

if __name__=="__main__":

    gaia_binaries_plot()

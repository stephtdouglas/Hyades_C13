import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import palettable
import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
from scipy.interpolate import interp1d

# from astropy.coordinates import SkyCoord  # High-level coordinates
# import astropy.units as u

from hypra.utils import cat_io, cat_match
from hypra.plot import color_mag

hdat,hobs,hobsnr,hobsr = cat_io.get_data("H")

# A few definitions
single_figure = (8,8)
double_figure = (13,6)
double_stacked = (6,13)
triple_stacked = (6,12)
quad_square = (13,13)

std_ms = 9
std_mew = 1.5

# Define colors, same as Praesepe paper
color_cycle = palettable.colorbrewer.qualitative.Dark2_8.mpl_colors
ms_color = color_cycle[3] #"#5d58a7" #
cand_color = color_cycle[2]
lit_color = color_cycle[1]
old_color = color_cycle[4]
conf_color = "k"
background_color = plt.cm.Greys(0.4)

gaia_tab = at.read("Gaia_Comb_Table.csv")

def calc_residual(abs_r, dmod):
    fig = plt.figure()
    ax = plt.subplot(111)

    ms,color = color_mag.add_ms(ax,dmod,return_ms=True,line_color=ms_color)
    plt.close()

    model_r = np.ones_like(abs_r)*-9999.
    good_interp = np.where((hdat["RPRIME_K"]>min(color)) & (hdat["RPRIME_K"]<max(color)))[0]
    model_r[good_interp] = ms(hdat["RPRIME_K"][good_interp])
    residual = abs_r - model_r
    residual[model_r<-9998] = -9999

    return residual

gfit = at.read(os.path.expanduser("~/my_papers/hyadesk22/NotesJC/Hyades-CMD-fit.txt"))
def calc_gaia_ms():
    fit_gmag, color = gfit["col1"], gfit["col2"]
    ms = interp1d(fit_gmag,color,fill_value=np.nan,bounds_error=False)
    return ms, color

def calc_gaia_residual(abs_g,bp_rp,dmod):

    ms,color = calc_gaia_ms()
    model_g = np.ones_like(abs_g)*-9999.
    good_interp = np.where((bp_rp>min(color)) & (bp_rp<max(color)) &
                           (bp_rp.mask==False))[0]
    bp_rp_for_interp = np.asarray(bp_rp[good_interp])
    model_g[good_interp] = ms(bp_rp_for_interp)
    residual = abs_g - model_g
    residual[model_g<-9998] = np.nan

    return residual

def plot_gaia_ms(ax,ms_color):

    fit_gmag, color = gfit["col1"], gfit["col2"]
    ax.plot(fit_gmag,color,"-",color=ms_color)
    ax.plot(fit_gmag,color-0.75/2,"-.",color=ms_color)
    ax.plot(fit_gmag,color-0.75,":",color=ms_color)




def calc_old_binaries():
    dmod = 5 * np.log10(hdat["DISTANCE"]) - 5
    abs_r = hdat["RPRIME"] - dmod
    abs_r[hdat["RPRIME"]<-9998] = -9999
    residual = calc_residual(abs_r,dmod)
    cand = ((hdat["RPRIME_K"]>0) & (residual>-9) &
            (hdat["RPRIME_K"]<4) & (residual<(-0.75/2)))

    return cand


gaia_jason = at.read(os.path.expanduser("~/my_papers/hyadesk22/NotesJC/Table-Hyades.txt"))
def gaia_candidates():
    jason_cand = gaia_jason["IDX"][gaia_jason["dCMD"]<(-0.375)]
    print(jason_cand)

    return jason_cand


def plot_binaries(ax, x, y, plot_all=False, plot_confirmed=True,
                  plot_candidate=True, plot_single=False,
                  plot_old_candidate=False,plot_gaia_candidate=False):

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



    our_cand0 = calc_old_binaries()
    # print(our_cand0)
    our_cand = np.asarray(our_cand0,int)
    # print(our_cand)

    # our_cand = np.zeros(len(hdat),int)
    # for i,cite in enumerate(hdat["BINARY_CITE"]):
    #     if "douglas" in cite:
    #         our_cand[i] = 1
    #     else:
    #         continue
    # Candidate binaries
    if plot_candidate:
        binary = (x>0) & (y>-9) & (hdat["BINARY"]==1)
        ours = our_cand==1
        ax.plot(x[binary & ours], y[binary & ours], 'o', mfc="none", mew=1.5,
                mec=cand_color, label="Candidate D16")
        other = our_cand==0
        nother = len(np.where(other & binary)[0])
        print(nother)
        if nother>0:
            ax.plot(x[binary & other], y[binary & other], 'D', mfc="none", mew=1.5,
                    mec=lit_color, label="Candidate (Literature)")

    if plot_old_candidate:
        ax.plot(x[our_cand0], y[our_cand0], 'o', mfc="none", mew=1.5,
                mec=old_color, label="Candidate (D16)")

    if plot_gaia_candidate:
        gaia_cand = gaia_candidates()
        ax.plot(x[gaia_cand], y[gaia_cand], 'o', mfc="none", mew=1.5,
                mec='blue', label="Candidate (Gaia CMD)",zorder=12)


    ax.tick_params(width=1.5, length=6)
    ax.tick_params(which="minor",width=1.5, length=3)


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
    plot_binaries(ax, hdat["RPRIME_K"], residual, True, False, False,
                  plot_old_candidate=True)
    # cand = np.where((hdat["RPRIME_K"]>0) & (residual>-9) & (hdat["RPRIME_K"]<4) & (residual<(-0.75/2)))[0]
    # ax.plot(hdat["RPRIME_K"][cand], residual[cand], 'o', mfc="none", mew=1.5,
    #         mec=cand_color, label="Candidate Multiple")
    ax.axvline(4,color="Grey",zorder=-111)
    ax.set_xlabel("")
    ax.tick_params(labelbottom=False)
    plt.legend(loc=1)

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



def paper_plot_gaia(abs_g,color,dmod,color_name=r"(G$_{BP}$ - G$_{RP}$)",extents=[0.75,4,14.5,-1]):

    residual = calc_gaia_residual(abs_g,color,dmod)
    fig = plt.figure(figsize=triple_stacked)

    #### TOP - CMD with binary sequence overlaid
    ax = plt.subplot(311)
    #ax = plt.subplot2grid((3,1),(0,0),rowspan=2)
    color_mag.setup_axes(ax,color_name,r"M$_{G}$",
                        extents)

    # TODO: add spectral types in Gaia colors
    # color_mag.add_spt(ax, 0.2)

    plot_binaries(ax, color, abs_g, True,False,False)
    ax.set_xlabel("")
    ax.tick_params(labelbottom=False)
    # Plot main sequence for Gaia
    plot_gaia_ms(ax,ms_color)


    l1 = ax.axhline(-99, color=ms_color, ls="-", lw=2, label="Main Sequence")
    l2 = ax.axhline(-99, color=ms_color, ls="-.", lw=2, label="Candidate Binary Threshold")
    l3 = ax.axhline(-99, color=ms_color, ls=":", lw=2, label="Binary Main Sequence")
    # l4 = ax.axvline(-99,color="Grey",zorder=-111,label="Candidate Binary Red Limit")
    leg = ax.legend(handles=[l1,l3,l2],loc=1,borderaxespad=0)
    leg.get_frame().set_alpha(1)



    #### MIDDLE - residuals with candidates identified
    ax = plt.subplot(312)
    #ax = plt.subplot2grid((3,1),(2,0),rowspan=1)
    color_mag.setup_axes(ax,color_name,
                    r"M$_{G}$ - M$_{G}$(Main Sequence)",
                    extents)
    plt.subplots_adjust(hspace=0)
    plot_binaries(ax, color, residual, True, False, False,
                  plot_old_candidate=True)
    # cand = np.where((hdat["RPRIME_K"]>0) & (residual>-9) & (hdat["RPRIME_K"]<4) & (residual<(-0.75/2)))[0]
    # ax.plot(hdat["RPRIME_K"][cand], residual[cand], 'o', mfc="none", mew=1.5,
    #         mec=cand_color, label="Candidate Multiple")
    ax.set_xlabel("")
    ax.tick_params(labelbottom=False)
    plt.legend(loc=1)

    l1 = ax.axhline(0, color=ms_color, ls="-", lw=2, label="Main Sequence")
    l2 = ax.axhline(-0.75/2, color=ms_color, ls="-.", lw=2, label="Candidate Binary Threshold")
    l3 = ax.axhline(-0.75, color=ms_color, ls=":", lw=2, label="Binary Main Sequence")
    ax.set_ylim(1.5,-4)


    #### BOTTOM- residuals with confirmed identified
    ax = plt.subplot(313)
    #ax = plt.subplot2grid((3,1),(2,0),rowspan=1)
    color_mag.setup_axes(ax,color_name,
                    r"M$_{G}$ - M$_{G}$(Main Sequence)",
                    extents)
    plt.subplots_adjust(hspace=0)
    plot_binaries(ax, color, residual, True)
    print(len(np.where(residual>-1)[0]))
    leg = ax.legend(loc=1, numpoints=1,borderaxespad=0)
    leg.get_frame().set_alpha(1)

    l1 = ax.axhline(0, color=ms_color, ls="-", lw=2, label="Main Sequence")
    l2 = ax.axhline(-0.75/2, color=ms_color, ls="-.", lw=2, label="Candidate Binary Threshold")
    l3 = ax.axhline(-0.75, color=ms_color, ls=":", lw=2, label="Binary Main Sequence")
    # l4 = ax.axvline(4,color="Grey",zorder=-111,label="Candidate Binary Red Limit")

    ax.set_ylim(1.5,-4)

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

    plt.suptitle("Gaia parallax for all",y=0.92)
    plt.savefig("binary_cmd_gaiadist_oldbin.png",bbox_inches="tight")
    plt.close("all")
    # plt.show()

    abs_g = gaia["Gmag"] - dmod
    bprp = gaia["BPmag"] - gaia["RPmag"]
    # bprp_name = r"(G$_{BP}$ - G$_{RP}$)"
    paper_plot_gaia(abs_g,bprp,dmod)
    plt.suptitle("Gaia parallax & Gaia photometry",y=0.92)
    plt.savefig("binary_cmd_gaiadist_bprp.png",bbox_inches="tight")
    # plt.show()
    plt.close("all")



def jason_binaries_plot():
    dmod = 5 * np.log10(hdat["DISTANCE"]) - 5
    abs_r = hdat["RPRIME"] - dmod
    abs_r[hdat["RPRIME"]<-9998] = -9999

    paper_plot(abs_r,dmod)

    plt.suptitle("Jason's gaia matches, old parallaxes",y=0.9)
    # plt.savefig("binary_cmd_gaiadist_oldbin.png",bbox_inches="tight")
    plt.show()
    # plt.close("all")


def old_binaries_plot():
    dmod = 5 * np.log10(hdat["DISTANCE"]) - 5
    abs_r = hdat["RPRIME"] - dmod
    abs_r[hdat["RPRIME"]<-9998] = -9999

    paper_plot(abs_r,dmod)
    plt.suptitle("D16 parallax for all",y=0.92)
    plt.savefig("binary_cmd_olddist_oldbin.png",bbox_inches="tight")
    plt.close("all")



if __name__=="__main__":
    old_binaries_plot()
    gaia_binaries_plot()
    # jason_binaries_plot()

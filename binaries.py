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

import add_gaia_spts

hdat,hobs,hobsnr,hobsr = cat_io.get_data("H")

# A few definitions
single_figure = (8,8)
double_figure = (13,6)
double_stacked = (6,13)
triple_stacked = (6,13)
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

def calc_new_binaries():
    gaia = at.read("Gaia_Comb_Table.csv")
    # New binaries with Gaia parallaxes where possible
    gaia_plx = np.zeros(len(hdat))*np.nan
    # good_plx = gaia["Plx"].mask==False
    # Only want to use the Gaia parallax when it meets the quality criteria!
    good_plx = ((gaia_tab["GAIA_QUAL"]=="True")
                & (gaia_tab["GAIA_QUAL"].mask==False))

    gaia_plx[gaia["HYADES_IDX"][good_plx]] = gaia["Plx"][good_plx]

    gaia_dist = 1000/gaia_plx

    dmod = 5 * np.log10(hdat["DISTANCE"]) - 5
    dmod[np.isfinite(gaia_plx)] = 5 * np.log10(gaia_dist[np.isfinite(gaia_plx)]) - 5
    abs_r = hdat["RPRIME"] - dmod
    abs_r[hdat["RPRIME"]<-9998] = -9999

    new_residual = calc_residual(abs_r,dmod)
    new_cand = ((hdat["RPRIME_K"]>0) & (new_residual>-9) &
                (hdat["RPRIME_K"]<4) & (new_residual<(-0.75/2)))
    # print("New",np.where(new_cand))

    return new_cand

gaia_jason = at.read(os.path.expanduser("~/my_papers/hyadesk22/NotesJC/Table-Hyades.txt"))
def gaia_candidates():

    jason_cand = gaia_jason["IDX"][gaia_jason["dCMD"]<(-0.375)]
    print(jason_cand)
    print(gaia_jason["dCMD"][gaia_jason["dCMD"]<(-0.375)])

    return jason_cand

# def compare_gaia_residuals()
#
#     gaia_cand = gaia_candidates()
#
#     gaia = at.read("Gaia_Comb_Table.csv")
#     # print(gaia["GAIA_QUAL"])
#
#     gaia_plx = np.zeros(len(hdat))*np.nan
#     # good_plx = gaia["Plx"].mask==False
#     # Only want to use the Gaia parallax when it meets the quality criteria!
#     good_plx = ((gaia_tab["GAIA_QUAL"]=="True")
#                 & (gaia_tab["GAIA_QUAL"].mask==False))
#
#     gaia_plx[gaia["HYADES_IDX"][good_plx]] = gaia["Plx"][good_plx]
#
#     gaia_dist = 1000/gaia_plx
#
#     dmod = 5 * np.log10(hdat["DISTANCE"]) - 5
#     dmod[np.isfinite(gaia_plx)] = 5 * np.log10(gaia_dist[np.isfinite(gaia_plx)]) - 5
#     abs_g = gaia["Gmag"] - dmod
#     bprp = gaia["BPmag"] - gaia["RPmag"]


def compare_candidates():
    d16_cand = np.where(calc_old_binaries())[0]
    jason_cand = gaia_candidates()

    removed_cand = np.setdiff1d(d16_cand,jason_cand)
    print("Removed: ",removed_cand)
    removed_jidx = np.ones(len(removed_cand),int)*-99
    for i,ridx in enumerate(removed_cand):
        jloc = np.where(gaia_jason["IDX"]==ridx)[0]
        if len(jloc)==1:
            # print(jloc)
            removed_jidx[i] = jloc[0]
    removed_jidx = removed_jidx[removed_jidx>=0]
    # print(removed_jidx)
    # print(gaia_jason["FinalFlag"][removed_jidx])

    # Note: I was afraid that some of the old photometric candidates are
    # no longer candidates after Gaia. That's true! But Jason wasn't removing
    # any candidates from his analysis, so this only affects my plots
    j1 = np.where(gaia_jason["FinalFlag"][removed_jidx]==1)[0]
    print("\nNow not a candidate, Jason is plotting:",len(j1))
    j1_jidx = np.intersect1d(np.where(gaia_jason["FinalFlag"]==1)[0],
                            removed_jidx)
    j1_idx = gaia_jason["IDX"][j1_jidx]
    print(j1_idx)
    print(hdat["EPIC_ID"][j1_idx])
    print(hdat["PERIOD"][j1_idx])
    print(hdat["RPRIME_K"][j1_idx])
    print(hdat["KH_MASS"][j1_idx])
    print(hdat["BINARY"][j1_idx])
    # j0 = np.where(gaia_jason["FinalFlag"][removed_jidx]==0)[0]
    # print("\nNow not a candidate, Jason was not plotting:",len(j0))
    # print(hdat["BINARY"][j0])

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

    new_cand0 = calc_new_binaries()
    new_cand = np.asarray(new_cand0,int)

    # our_cand = np.zeros(len(hdat),int)
    # for i,cite in enumerate(hdat["BINARY_CITE"]):
    #     if "douglas" in cite:
    #         our_cand[i] = 1
    #     else:
    #         continue
    # Candidate binaries
    if plot_candidate:
        binary = (x>0) & (y>-9) & (new_cand==1)
        if plot_confirmed:
            binary = binary & (hdat["BINARY"]<2)
        ax.plot(x[binary], y[binary], 'o', mfc="none", mew=1.5,
                mec=cand_color, label="Candidate binary",zorder=12)

        if plot_confirmed:
            lit_cand = (x>0) & (y>-9) & (hdat["BINARY"]==1) & (our_cand==0)
            nother = len(np.where(lit_cand)[0])
            print(nother)
            if nother>0:
                ax.plot(x[lit_cand], y[lit_cand], 'D', mfc="none", mew=1.5,mec=lit_color, label="Candidate (Literature)")


    if plot_old_candidate:
        ax.plot(x[our_cand0], y[our_cand0], 'o', mfc="none", mew=1.5,
                mec=old_color, label="Candidate (D16) - removed",
                zorder=-11)

    if plot_gaia_candidate:
        gaia_cand = gaia_candidates()
        print(y[gaia_cand])
        print(x[gaia_cand])
        print(hdat["KH_MASS"][gaia_cand])
        if plot_confirmed:
            gaia_cand = np.intersect1d(gaia_cand, np.where(hdat["BINARY"]<2)[0])
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
    plot_binaries(ax, hdat["RPRIME_K"], residual, True, False, True,
                  plot_old_candidate=True)
    # cand = np.where((hdat["RPRIME_K"]>0) & (residual>-9) & (hdat["RPRIME_K"]<4) & (residual<(-0.75/2)))[0]
    # print("Plot",cand)
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
    plot_binaries(ax, hdat["RPRIME_K"], residual, True,
                  plot_candidate=True,
                  plot_old_candidate=False)
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



def paper_plot_gaia(abs_g, color, dmod,
                    color_name=r"(G$_{BP}$ - G$_{RP}$)",
                    extents=[0.75,4,14.5,-1], plot_single=False,
                    plot_old_candidate=False, plot_gaia_candidate=False):

    residual = calc_gaia_residual(abs_g,color,dmod)
    fig = plt.figure(figsize=triple_stacked)

    #### TOP - CMD with binary sequence overlaid
    ax = plt.subplot(311)
    #ax = plt.subplot2grid((3,1),(0,0),rowspan=2)
    color_mag.setup_axes(ax,color_name,r"M$_{G}$",
                        extents)

    # add spectral types in Gaia colors
    add_gaia_spts.add_gaia_spts(ax,texty=-1.4)

    plot_binaries(ax, color, abs_g, True,False,False, plot_single=False,
                  plot_old_candidate=False,
                  plot_gaia_candidate=False)
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
    plot_binaries(ax, color, residual, True, False, True,
                  plot_single=plot_single,
                  plot_old_candidate=True,
                  plot_gaia_candidate=plot_gaia_candidate)
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
    plot_binaries(ax, color, residual, True,plot_single=plot_single,
                  plot_old_candidate=False,
                  plot_gaia_candidate=plot_gaia_candidate)
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

def compute_new_rpmK_binaries():
    gaia = at.read("Gaia_Comb_Table.csv")

    # Old binaries
    old_cand = calc_old_binaries()

    # New old_binaries
    new_cand = calc_new_binaries()


    remove_candidates = np.where((old_cand==True) & (new_cand==False))[0]
    remove_final = np.copy(remove_candidates)
    for i in remove_candidates:
        print(i,hdat["BINARY"][i],hdat["BINARY_CITE"][i])
        if hdat["BINARY"][i]!=1:
            # it's confirmed or wasn't actually a candidate
            loc = np.where(remove_final==i)[0]
            remove_final = np.delete(remove_final,loc)
            print("Keep ^^, confirmed or single\n")
        elif ((hdat["BINARY"][i]==1) and (old_cand[i]==True) and
            (("douglas2014,"==hdat["BINARY_CITE"][i])==False)):
            # if we were the only ones
            loc = np.where(remove_final==i)[0]
            remove_final = np.delete(remove_final,loc)
            print("Keep ^^, there's another ref\n")
        else:
            continue

    print(remove_final)
    print(hdat["BINARY"][remove_final])

    add_candidates = np.where((old_cand==False) & (new_cand==True)
                              & (hdat["BINARY"]<1))[0]
    print(add_candidates)
    print(hdat["BINARY"][add_candidates])

    return remove_final, add_candidates, np.where(new_cand)[0]

def gaia_binaries_plot():

    gaia = at.read("Gaia_Comb_Table.csv")
    # print(gaia["GAIA_QUAL"])

    gaia_plx = np.zeros(len(hdat))*np.nan
    # good_plx = gaia["Plx"].mask==False
    # Only want to use the Gaia parallax when it meets the quality criteria!
    good_plx = ((gaia_tab["GAIA_QUAL"]=="True")
                & (gaia_tab["GAIA_QUAL"].mask==False))

    gaia_plx[gaia["HYADES_IDX"][good_plx]] = gaia["Plx"][good_plx]

    gaia_dist = 1000/gaia_plx

    dmod = 5 * np.log10(hdat["DISTANCE"]) - 5
    dmod[np.isfinite(gaia_plx)] = 5 * np.log10(gaia_dist[np.isfinite(gaia_plx)]) - 5
    abs_r = hdat["RPRIME"] - dmod
    abs_r[hdat["RPRIME"]<-9998] = -9999

    paper_plot(abs_r,dmod)

    plt.suptitle("Gaia parallax where high-quality",y=0.92)
    plt.savefig("binary_cmd_gaiadist_oldbin.png",bbox_inches="tight")
    plt.savefig(os.path.expanduser("~/my_papers/hyadesk22/binary_cmd_rK.eps"),
                bbox_inches="tight")
    # plt.savefig("/home/stephanie/my_papers/praeK2/fig6.eps",
    #             bbox_inches="tight")
    # plt.savefig("/home/stephanie/my_papers/praeK2/fig6.pdf",
    #             bbox_inches="tight")
    # plt.savefig("/home/stephanie/Dropbox/plots_for_sharing/binary_cmd_praesepe.png",
    #             bbox_inches="tight"))    plt.close("all")
    # plt.show()


    gaia_plx = np.zeros(len(hdat))*np.nan
    good_plx = gaia["Plx"].mask==False
    # Jason used all Gaia parallaxes

    gaia_plx[gaia["HYADES_IDX"][good_plx]] = gaia["Plx"][good_plx]

    gaia_dist = 1000/gaia_plx

    dmod = 5 * np.log10(hdat["DISTANCE"]) - 5
    dmod[np.isfinite(gaia_plx)] = 5 * np.log10(gaia_dist[np.isfinite(gaia_plx)]) - 5

    abs_g = gaia["Gmag"] - dmod
    bprp = gaia["BPmag"] - gaia["RPmag"]
    # bprp_name = r"(G$_{BP}$ - G$_{RP}$)"
    paper_plot_gaia(abs_g,bprp,dmod,plot_single=False,
                    plot_old_candidate=False, plot_gaia_candidate=False)
    plt.suptitle("Gaia parallax & Gaia photometry",y=0.92)
    plt.savefig("binary_cmd_gaiadist_bprp.png",bbox_inches="tight")
    # plt.show()
    plt.close("all")


    abs_g = gaia["Gmag"] - dmod
    bprp = gaia["BPmag"] - gaia["RPmag"]
    # bprp_name = r"(G$_{BP}$ - G$_{RP}$)"
    paper_plot_gaia(abs_g,bprp,dmod,plot_single=False,
                    plot_old_candidate=False, plot_gaia_candidate=True)
    plt.suptitle("Gaia parallax & Gaia photometry",y=0.92)
    plt.savefig("binary_cmd_gaiadist_bprp_gaiacand.png",bbox_inches="tight")
    plt.savefig(os.path.expanduser("~/my_papers/hyadesk22/binary_cmd_bprp.eps"),
                bbox_inches="tight")
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
    # old_binaries_plot()
    gaia_binaries_plot()
    # compare_candidates()
    # compute_new_rpmK_binaries()

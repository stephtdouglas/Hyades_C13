from __future__ import division, print_function
import sys, os

import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
import astropy.table as table
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

from hypra.utils import cat_match, cat_io, k2utils
# from hypra.plot import color_mag
import convertmass

def get_c13_periods():

    hdat,hobs,hobsnr,hobsr = cat_io.get_data("H")
    hpos = SkyCoord(hdat["RA"],hdat["DEC"],unit=u.degree)

    # Match EPIC IDs - this needs to be moved to the catalog script
    epic_list = at.read("mast_search_13064.csv",data_start=2)
    # print(epic_list.dtype)
    mpos = SkyCoord(epic_list["RA (J2000)"], epic_list["Dec (J2000)"],
                   unit=(u.hourangle,u.degree))

    idx, sep, _ = hpos.match_to_catalog_sky(mpos)
    # print(len(idx),len(hpos),len(mpos))

    good_match = np.where(sep<(5*u.arcsec))[0]
    good_idx = idx[good_match]
    # print(len(good_match),len(good_idx),len(np.unique(good_idx)))
    hdat["EPIC_ID"][good_match] = epic_list["K2 ID"][good_idx]

    # Basic code output
    period_file_base = "c13_tables/c13_k2sc_output_2018-05-29_"

    for i in range(4):
        #print(i+1)
        period_file = "{0}{1}.csv".format(period_file_base,i+1)
        new_res = at.read(period_file)

        peak_file = "{0}{1}_allpeaks.csv".format(period_file_base,i+1)
        new_peaks = at.read(peak_file)
        if i==0:
            res1 = new_res
            peaks1 = new_peaks
        else:
            res1 = table.vstack([res1,new_res])
            peaks1 = table.vstack([peaks1,new_peaks])

    res = res1 #table.vstack([res1,res2])
    peaks = peaks1 #table.vstack([peaks1,peaks2])


    # print(len(res))
    # print(len(peaks),len(np.unique(peaks["EPIC"])))
    # print(peaks.dtype)

    # Read in comments file
    # Ensure that this is still sorted by EPIC!
    comments = at.read("c13_tables/c13_k2sc_output_2018-05-29_comments_final.csv",
                       delimiter=",")
    # print(comments.dtype)
    print(len(comments),"K2 targets")
    #
    print(np.where(comments["Q"]==-1)[0])

    res.sort("EPIC")
    comments.sort("EPIC")
    epic_diff = res["EPIC"] - comments["EPIC"]
    print(np.where(epic_diff<0)[0])

    # print(res.dtype)
    k2_periods = res["sig_period"]
    k2_power = res["sig_power"]
    k2_epic = res["EPIC"]
    k2_flag = res["num_sig"]
    #k2_kpmag = res["magnitude"]
    k2_harm = res["harm_type"]
    k2_threshold = res["threshold"]
    k2_quality = comments["Q"]
    k2_spot = comments["SpotEvol"]
    k2_multi = comments["MultiProt"]
    k2_comment = comments["Notes"]
    k2_blend = comments["Blended"]

    pperiods = np.copy(hdat['PERIOD'])
    pperiods_secondary = np.zeros(len(hdat))
    powers_secondary = np.zeros(len(hdat))
    quality_secondary = np.zeros(len(hdat),int)
    pflag = np.copy(hdat['PERIOD_FLAG'])
    pmass = hdat['KH_MASS']
    pqual = np.zeros(len(hdat),int)
    pbad = np.zeros(len(hdat),int)
    ppower = np.zeros(len(hdat))
    pthreshold = np.zeros(len(hdat))
    pharm = np.empty(len(hdat),"S12")
    pharm[:] = ""
    pblend = np.zeros(len(hdat),int)

    bad_count = 0
    replace_count = 0

    # Set up for cross-matching with new periods
    potential_new = np.where((hdat["EPIC_ID"]>0) &
                             (hdat["KH_MASS"]<1.5) & (hdat["KH_MASS"]>0)
                               )[0]

    pperiods_allk2 = np.copy(pperiods)
    pperiods_onlyk2 = np.ones_like(pperiods)*-9999.0

    # peak separation threshold to be considered a binary
    min_dpp = 0.2
    close_peak_count = 0
    blend_count = 0

    # Match targets from our catalog to the EPIC results
    for i in potential_new:

        cat_epic = hdat["EPIC_ID"][i]
        k2_loc = np.where(k2_epic==cat_epic)[0]
        if len(k2_loc)>=1:
    #         print(i, k2_loc, k2_epic[k2_loc], hdat["EPIC_ID"][i])

            k2_loc = k2_loc[0]
            if ((k2_quality[k2_loc]==2) and (k2_periods[k2_loc]<0)) | (k2_quality[k2_loc]==3):
                # Indicates that no period was detected by the lomb-scargle periodogram
                bad_count += 1
                pqual[i] = 3
                continue

            # By default, use the highest significant peak in the light curve
            use_period = k2_periods[k2_loc]
            use_power = k2_power[k2_loc]
            pqual[i] = comments["Q"][k2_loc]
    #         print(c5_loc, c5_comment[c5_loc])

            # If I've flagged that peak as bad, then check my flag for the second highest peak
            # If the second period is good or OK, use that one
            if comments["Q"][k2_loc]==2:
                if comments["Q2"][k2_loc] in ["0","1"]:
                    use_period = res["sec_period"][k2_loc]
                    use_power = res["sec_power"][k2_loc]
                    pqual[i] = comments["Q2"][k2_loc]
                    replace_count += 1

            # If the first peak is good or OK, then check to see if the second or third peak
            #  is really a second period
            # AND (new criterion) the peaks must be separated by a certain fraction of the dominant peak
            else:
                if comments["Q2"][k2_loc] in ["0","1"]:

                    delta_pp = abs(res["sec_period"][k2_loc] - use_period)/use_period

                    if delta_pp > min_dpp:
                        pperiods_secondary[i] = res["sec_period"][k2_loc]
                        powers_secondary[i] = res["sec_power"][k2_loc]
                        quality_secondary[i] = comments["Q2"][k2_loc]
                    else:
                        close_peak_count += 1

                # I only saved the first and second periods/powers in my main table,
                # so have to dig into the table of all significant peaks to pull out
                # the third period/power
                elif comments["Q3"][k2_loc] in ["0","1"]:
                    epic_peaks = np.where(cat_epic==peaks["EPIC"])[0]
                    sorted_heights = np.argsort(peaks["power"][epic_peaks])
                    sorted_periods = peaks["period"][epic_peaks][sorted_heights]
                    sorted_powers = peaks["power"][epic_peaks][sorted_heights]

                    delta_pp = abs(sorted_periods[2] - use_period)/use_period

                    if delta_pp > min_dpp:
                        pperiods_secondary[i] = sorted_periods[2]
                        powers_secondary[i] = sorted_powers[2]
                        quality_secondary[i] = comments["Q3"][k2_loc]
                    else:
                        close_peak_count += 1

            # Check whether there's already a literature period in my catalog
            # If not, K2 is the dominant period
            if (pflag[i]=="-"):
                pperiods[i] = use_period
                pflag[i] = "3"

            # Save the periods to the appropriate arrays and set various flags
            pperiods_allk2[i] = use_period
            pperiods_onlyk2[i] = use_period
            ppower[i] = use_power
            pthreshold[i] = k2_threshold[k2_loc]
    #         pbad[i] = c5_flag[c5_loc]

            # If there is definitely a second period or blended neighbor,
            # flag as a candidate binary
            if ((k2_multi[k2_loc]=="y") or (k2_multi[k2_loc]=="Y")
               or (k2_blend[k2_loc]=="Y") or (k2_blend[k2_loc]=="y")
               ):
                pblend[i] = 1
                if (k2_blend[k2_loc]=="Y") or (k2_blend[k2_loc]=="y"):
                    blend_count += 1


            # Reset quality flag based on number of peaks in the light curve
            epic_peaks = np.where(cat_epic==peaks["EPIC"])[0]
            if len(epic_peaks)>0:
                epic_heights = peaks["power"][epic_peaks]
                pbad[i] = len(np.where(epic_heights>(res["fund_power"][k2_loc]*0.6))[0])
                if pbad[i]>0:
                    pbad[i] -= 1
    #             print(c5_flag[c5_loc],pbad[i])

            #print(c5_harm[c5_loc],pharm[i])
            pharm[i] = k2_harm[k2_loc][0]


    ntargets = len(comments)
    print(len(np.where((pperiods_onlyk2>0) & (pqual<2) & (pqual>=0))[0]),"K2 stars with periods")


    print(len(np.where(pperiods_secondary>0)[0]),"with secondary periods",
          "at peak separations >",min_dpp)
    print("plus",close_peak_count,"with closer secondary peaks")
    print(blend_count,"with blended neighbors, or {0:.1f}\%".format(100*blend_count/len(comments)))
    print(len(np.where(pblend>0)[0]),"blend/multi")

    # Counting - directly adapted from C5, possibly not the right numbers
    # New Periods
    new_k2 = np.where((pperiods_onlyk2>0) & (hdat["PERIOD"]<0) & (pqual<2) & (pqual>=0))[0]
    print(len(new_k2), " new K2 periods")
    repeat_k2 = np.where((pperiods_onlyk2>0) & (hdat["PERIOD"]>0) & (pqual<2) & (pqual>=0))[0]
    print(len(repeat_k2), " with literature periods")
    quality_k2 = np.where((pperiods_onlyk2>0) & (hdat["PERIOD"]<0) & (pqual==0) & (pbad==0))[0]
    print(len(quality_k2), " high quality new K2 periods")
    low_quality_k2 = np.where((pperiods_onlyk2>0) & (hdat["PERIOD"]<0) & (pqual==1))[0]
    print(len(low_quality_k2), " low quality new K2 periods as flagged by me")
    low_quality_k2 = np.where((pperiods_onlyk2>0) & (hdat["PERIOD"]<0) & (pqual==0) & (pbad>0))[0]
    print(len(low_quality_k2), " low quality new K2 periods from not-clean periodograms (but I flagged them good)")
    low_quality_k2 = np.where((pperiods_onlyk2>0) & (hdat["PERIOD"]<0) & (pqual==1) & (pbad>0))[0]
    print(len(low_quality_k2), " low quality new K2 periods as flagged by me AND the periodograms\n")

    # Non-detections
    print(bad_count," with no K2 detection",bad_count*100/ntargets,"%")
    print(len(np.where(pqual==3)[0]),"alternate count of no detections\n")

    # Ones I flagged as bad
    bad_flag_count = len(np.where((pqual==2))[0])
    print("I removed",bad_flag_count,bad_flag_count*100/ntargets,"% leaving",ntargets-bad_count-bad_flag_count)


    # 1st peak isn't real period
    print(replace_count," stars where the highest peak doesn't look real")

    # Totals
    print("Originally",ntargets,"remove", bad_count, "and", bad_flag_count, "leave", ntargets-bad_count-bad_flag_count,"\n")

    # Total rotators
    print("Pre-K2:",len(np.where(hdat["PERIOD"]>0)[0]))
    print("Total rotators:",len(np.where((hdat["PERIOD"]>0) |
          ((pperiods_onlyk2>0) & (pqual<2) & (pqual>=0)))[0]))

    # print(np.intersect1d(c13,np.where((pperiods_onlyk2<0) & ((pqual==1) | pqual==0))[0]))
    # print(np.intersect1d(c13,np.where((pperiods_allk2<0) & ((pqual==1) | pqual==0))[0]))

    output = (pperiods,pperiods_secondary,powers_secondary,
              quality_secondary,pflag,pmass,pqual,pbad,ppower,pthreshold,
              pharm,pblend,pperiods_allk2,pperiods_onlyk2)

    return output

if __name__=="__main__":
    output = get_c13_periods()

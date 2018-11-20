import os

import matplotlib
matplotlib.use("agg")
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.io.ascii as at
from astropy.table import Table

from hypra.utils import cat_match, cat_io, k2utils

if __name__=="__main__":
    hdat,hobs,hobsnr,hobsr = cat_io.get_data("H")
    hpos = SkyCoord(hdat["RA"],hdat["DEC"],unit=u.degree)
    # print(hdat["MDM_SPEC_FILE"][hdat["MDM_SPEC_FILE"]!="no string"])
    # print(hdat["TWOMASSNAME"])

    cat_files = hdat["MDM_SPEC_FILE"]
    for i, fname in enumerate(cat_files):
        cat_files[i] = cat_files[i].split("/")[-1]
    cat_files2 = hdat["MDM_SPEC_FILE2"]
    for i, fname in enumerate(cat_files2):
        cat_files2[i] = cat_files2[i].split("/")[-1]

    # Match EPIC IDs
    epic_list = at.read("mast_search_13064.csv",data_start=2)
    # print(epic_list.dtype)
    mpos = SkyCoord(epic_list["RA (J2000)"], epic_list["Dec (J2000)"],
                   unit=(u.hourangle,u.degree))

    idx, sep, _ = hpos.match_to_catalog_sky(mpos)
    print(len(idx),len(hpos),len(mpos))

    good_match = np.where(sep<(5*u.arcsec))[0]
    good_idx = idx[good_match]
    print(len(good_match),len(good_idx),len(np.unique(good_idx)))

    hdat["EPIC_ID"][good_match] = epic_list["K2 ID"][good_idx]


    # Read in the EqWs list
    eqws = Table(at.read("halpha_equivalent_widths.csv"))
    eqws["HYADES_IDX"] = np.zeros(len(eqws),int)*-99
    eqws["MDM_NO"] = np.zeros(len(eqws),int)

    for row in eqws:
        fname = row["filename"].split("/")[-1]
        cat_name = fname.replace("_update","")

        # First try just matching the raw filename, for old data
        cat_loc = np.where(cat_name==cat_files)[0]
        match = False
        cat_loc2 = np.where(cat_name==cat_files2)[0]
        if len(cat_loc)==1:
            row["HYADES_IDX"] = cat_loc[0]
            row["MDM_NO"] = 1
            # print("Match! 1",cat_name)
            match = True
        elif len(cat_loc2)==1:
            row["HYADES_IDX"] = cat_loc2[0]
            row["MDM_NO"] = 2
            # print("Match! 2",cat_name)
            match = True
        elif len(cat_loc)>1:
            print("multiple matches! 1",cat_name)
        elif len(cat_loc2)>1:
            print("multiple matches! 2",cat_name)
        else:
            # print("no match",cat_name)
            match = False

        # Now try matching by 2MASS Names, for Non-K2 stars and others
        if match is False:
            name = fname.split(".")[1].replace("_update","").replace("_20161208","")
            if "EP" in name:
                epic = name.replace("EP","").replace("A","").replace("B",""
                                    ).replace("_comb","")
                epic = np.int32(epic)
                loc = np.where(epic==hdat["EPIC_ID"])[0]
                if len(loc)==1:
                    row["HYADES_IDX"] = loc[0]
                    match = True
                elif len(loc)>1:
                    print("multiple matches! EP",epic)
                elif len(loc)==1:
                    # match = False
                    print("no match! EP",epic)
            else:
                tmass = name.replace("EP","").replace("A","").replace("B",""
                                    ).replace("_comb","")
                loc = np.where(tmass==hdat["TWOMASSNAME"])[0]
                if len(loc)==1:
                    row["HYADES_IDX"] = loc[0]
                    match = True
                elif len(loc)>1:
                    print("multiple matches! 2MASS",tmass)
                elif len(loc)==1:
                    print("no match! 2MASS",tmass)

            if match is False:
                for i,cat_tmass in enumerate(hdat["TWOMASSNAME"]):
                    if tmass in cat_tmass:
                        row["HYADES_IDX"] = i
                        match = True
                        print("found!",tmass)
                        break

            if match is False:
                print("No match at all!",row['filename'])

    print(len(np.where(row["HYADES_IDX"]<0)[0]))

    at.write(eqws,"halpha_equivalent_widths_matched.csv",delimiter=",")

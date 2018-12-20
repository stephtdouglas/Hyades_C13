from __future__ import division, print_function
import sys, os

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as at
import astropy.io.fits as fits
from astropy.table import Table, hstack, join
from astropy.time import Time
from astropy.coordinates import SkyCoord, Distance
import astropy.units as u
import matplotlib.ticker
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.font_manager import FontProperties
import palettable

from hypra.utils import cat_match, cat_io
import convertmass


def match_gaia(to_plot=False):

    # # Just crossmatching against Gaia

    # In[2]:


    gaia = at.read("gaia_vizier_allhyads_full.csv",delimiter="|",
                  data_start=3)
    # print(gaia.dtype)
    gpos = SkyCoord(gaia["RAJ2000"],gaia["DEJ2000"],unit=u.degree)


    # In[3]:

    hdat,_,_,_ = cat_io.get_data("H")
    hpos = SkyCoord(hdat["RA"],hdat["DEC"],unit=u.degree)


    # In[4]:


    # Match EPIC IDs
    epic_list = at.read("mast_search_13064.csv",data_start=2)
    print(epic_list.dtype)
    mpos = SkyCoord(epic_list["RA (J2000)"], epic_list["Dec (J2000)"],
                   unit=(u.hourangle,u.degree))
    idx, sep, _ = hpos.match_to_catalog_sky(mpos)
    print(len(idx),len(hpos),len(mpos))

    good_match = np.where(sep<(5*u.arcsec))[0]
    good_idx = idx[good_match]
    print(len(good_match),len(good_idx),len(np.unique(good_idx)))

    hdat["EPIC_ID"][good_match] = epic_list["K2 ID"][good_idx]


    # In[5]:


    hdat = Table(hdat)
    hdat["HYADES_IDX"] = np.arange(len(hdat))


    # In[6]:


    idx, sep, _ = hpos.match_to_catalog_sky(gpos)
    print(len(idx),len(sep),len(hpos),len(gpos))


    # In[7]:


    good_match = np.where(sep<(5*u.arcsec))[0]
    good_idx = idx[good_match]

    print(len(good_idx),len(np.unique(good_idx)))
    print(max(good_idx))


    # In[8]:


    gaia["HYADES_IDX"] = np.ones(len(gaia),int)*-99
    gaia["HYADES_IDX"][good_idx] = good_match


    # In[9]:


    gaia = gaia[gaia["HYADES_IDX"]>=0]


    # In[10]:


    # gaia["HYADES_IDX"]


    # # Join Hyades and Gaia tables

    # In[11]:


    joint_tab = join(hdat,gaia,join_type="left",keys=["HYADES_IDX"],
                     table_names=["Douglas","GaiaDR2"])


    # In[12]:


    # joint_tab.dtype.names


    # In[13]:


    # Gaia filtering (from the DR2 HRD paper):
    #    RPlx > 10
    #    RFG>50
    #    RFBP>20
    #    RFRP>20
    #    E(BR/RP) < 1.3+0.06*(BP-RP)**2
    #    E(BR/RP) > 1.0+0.015*(BP-RP)**2 # phot_bp_rp_excess_factor
    #    Nper>8 # (visibility periods)
    #    chi2AL/(NgAL-5)<1.44*greatest(1,exp(-0.4*(Gmag-19.5)))

    #    astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))

    # - astrometric_excess_noise<1 criterion,
    #   but this is less optimised for the bright stars
    #   because of the degrees of freedom (DOF) issue (Lindegren et al. 2018, Appendix A).
    # photometric cuts may remove variable stars!


    # In[14]:


    phot_excess = ((joint_tab["E(BR/RP)"] < (1.3+0.06*joint_tab["BP-RP"]**2)) &
                   (joint_tab["E(BR/RP)"] > (1.0+0.015*joint_tab["BP-RP"]**2)))

    gexp = np.exp(-0.4*(joint_tab["Gmag"]-19.5))
    gexp[gexp<1] = 1
    astrom = (joint_tab["chi2AL"]/(joint_tab["NgAL"]-5)) < (1.44*gexp)

    hrd_filter = ((joint_tab["RPlx"]>10) & (joint_tab["RFG"]>50)
                 & (joint_tab["RFBP"]>20) & (joint_tab["RFRP"]>20)
                 & phot_excess & (joint_tab["Nper"]>8)
                 & astrom)


    # In[15]:


    joint_tab["GAIA_QUAL"] = hrd_filter


    # ## Also compare to the Gaia HRD cluster list

    # In[16]:


    gaia_hrd = Table(at.read("gaia_HRD_Hyades.tsv",delimiter="|",data_start=3))
    gaia_hrd.dtype

    print(len(gaia_hrd))
    print(len(np.intersect1d(np.asarray(joint_tab["Source"]),
                             np.asarray(gaia_hrd["Source"]))))

    joint_tab["HRD"] = np.zeros(len(hdat),bool)

    for sid in gaia_hrd["Source"]:
        loc = np.where(joint_tab["Source"]==sid)[0]
        if len(loc)==1:
            joint_tab["HRD"][loc] = True
        else:
            continue


    # # Comparing earlier values with Gaia

    # In[17]:

    if to_plot is True:
        good_hip = joint_tab["HIP_PARALLAX"]>-999
        plt.errorbar(joint_tab["Plx"][good_hip],joint_tab["HIP_PARALLAX"][good_hip],
                    yerr=joint_tab["HIP_PAR_ERR"][good_hip],xerr=joint_tab["e_Plx"][good_hip],
                    linewidth=0,elinewidth=1)
        good_gaia = good_hip & joint_tab["GAIA_QUAL"]
        plt.plot(joint_tab["Plx"][good_gaia],joint_tab["HIP_PARALLAX"][good_gaia],'.')

        x = np.linspace(0,35)
        plt.plot(x,x)
        plt.xlim(10,35)
        plt.ylim(10,35)

        plt.xlabel("Gaia Parallax")
        plt.ylabel("HIP Parallax")

        plt.show()


        # In[18]:


        good_g = joint_tab["GOLDMAN_PLX"]>-999
        plt.errorbar(joint_tab["Plx"][good_g],
                     joint_tab["GOLDMAN_PLX"][good_g],
                    yerr=joint_tab["GOLDMAN_E_PLX"][good_g],
                     xerr=joint_tab["e_Plx"][good_g],
                    linewidth=0,elinewidth=1)
        good_gaia = good_g & joint_tab["GAIA_QUAL"]
        plt.plot(joint_tab["Plx"][good_gaia],
                 joint_tab["GOLDMAN_PLX"][good_gaia],'.',alpha=0.25)

        x = np.linspace(0,60)
        plt.plot(x,x)
        plt.xlim(0,60)
        plt.ylim(0,60)

        plt.xlabel("Gaia Parallax")
        plt.ylabel("Goldman Parallax")

        plt.show()


        # # Check cross-matching with photometry

        # ## Tycho photometry (from 2MASS), take 1

        # In[19]:


        g_calc = np.ones(len(hdat))*np.nan


        # In[20]:


        tycho_b = joint_tab["TWOMASS_B"]
        tycho_v = joint_tab["TWOMASS_VR"]
        tycho_v[joint_tab["TWOMASS_PHOT_FLAG"]!="T"] = np.nan
        tycho_b[joint_tab["TWOMASS_PHOT_FLAG"]!="T"] = np.nan


        # In[21]:


        tgood = np.isfinite(tycho_b) & (joint_tab["TWOMASS_PHOT_FLAG"]=="T")


        # In[22]:


        bv = tycho_b - tycho_v
        # https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html#Ch5.T7
        g_m_v = -0.02051 - 0.2706*bv + 0.03394*bv**2 -0.05937*bv**3
        g_calc[tgood] = tycho_v[tgood] + g_m_v[tgood]
        g_tycho = np.copy(g_calc)


        # In[23]:


        plt.plot(joint_tab["Gmag"],g_calc,'.')
        plt.plot(joint_tab["Gmag"][joint_tab["GAIA_QUAL"]],
                 g_calc[joint_tab["GAIA_QUAL"]],
                '.')

        plt.xlabel("Gaia G")
        plt.ylabel("Calculated G")
        plt.show()


        # ## UCAC r, i photometry

        # In[24]:


        g_calc = np.ones(len(joint_tab))*np.nan

        ucac_r = joint_tab["UCAC_R"]
        ucac_i = joint_tab["UCAC_I"]
        ugood = (ucac_r>0) & (ucac_i>0)


        # In[25]:


        ri = ucac_r - ucac_i
        # https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html#Ch5.T7
        g_m_r = 0.0014891 + 0.36291*ri -0.81282*ri**2 + 0.0060376*ri**3
        g_calc[ugood] = ucac_r[ugood] + g_m_r[ugood]
        g_ucac = np.copy(g_calc)


        # In[26]:


        plt.plot(joint_tab["Gmag"],g_calc,'.')
        plt.plot(joint_tab["Gmag"][joint_tab["GAIA_QUAL"]],
                 g_calc[joint_tab["GAIA_QUAL"]],
                '.')
        plt.xlabel("Gaia G")
        plt.ylabel("Calculated G")
        plt.show()

    # In[27]:


        len(np.where(np.isnan(g_calc))[0])


        # ## SDSS r,i photometry

        # In[28]:


        g_calc = np.ones(len(joint_tab))*np.nan

        sdss_r = joint_tab["SDSS_R"]
        sdss_i = joint_tab["SDSS_I"]
        sgood = (sdss_r>14) & (sdss_i>14)


        # In[29]:


        ri = sdss_r - sdss_i
        # https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html#Ch5.T7
        g_m_r = 0.0014891 + 0.36291*ri -0.81282*ri**2 + 0.0060376*ri**3
        g_calc[sgood] = sdss_r[sgood] + g_m_r[sgood]
        g_sdss = g_calc


        # In[30]:


        plt.plot(joint_tab["Gmag"],g_calc,'.')
        plt.xlabel("Gaia G")
        plt.ylabel("Calculated G")
        plt.show()

        # In[31]:


        len(np.where(np.isnan(g_calc))[0])


        # In[32]:


        plt.plot(joint_tab["Gmag"],joint_tab["RPRIME"],'.')
        plt.xlabel("Gaia G")
        plt.ylabel("r'")
        plt.ylim(0,25)


        # ## 2MASS

        # In[33]:


        jmk = joint_tab["TWOMASS_J"] - joint_tab["TWOMASS_K"]
        need_jmk = joint_tab["TWOMASS_J"]<=0
        jmk[need_jmk] = joint_tab["GOLDMAN_JMAG"][need_jmk] - joint_tab["GOLDMAN_KMAG"][need_jmk]
        hmk = joint_tab["TWOMASS_H"] - joint_tab["TWOMASS_K"]
        hmk[need_jmk] = joint_tab["GOLDMAN_HMAG"][need_jmk] - joint_tab["GOLDMAN_KMAG"][need_jmk]


        # In[34]:


        print(max(jmk),max(hmk))


        # In[35]:


        calc_bprp = 0.20215 + 1.9561*jmk -0.69629*jmk**2 + 0.47633*jmk**3


        # In[36]:


        plt.plot(calc_bprp,joint_tab["BP-RP"],'.')


        # In[37]:


        calc_GmK = 0.23587 + 4.0548*jmk - 2.5608*jmk**2 + 2.2228*jmk**3 -0.54944*jmk**4


        # In[38]:


        kmag = hdat["TWOMASS_K"]
        kmag[kmag<=0] = hdat["GOLDMAN_KMAG"][kmag<=0]


        # In[39]:


        calc_G = calc_GmK + kmag


        # In[40]:


        plt.plot(joint_tab["Gmag"],calc_G,'.')
        plt.xlabel("Gaia G")
        plt.ylabel("Calculated G")
        plt.ylim(0,20)


        # In[41]:


        ucac_diff = abs(joint_tab["Gmag"] - g_ucac)
        sdss_diff = abs(joint_tab["Gmag"] - g_sdss)
        tycho_diff = abs(joint_tab["Gmag"] - g_tycho)
        gK_diff = abs(joint_tab["Gmag"]-calc_G)


        # In[42]:


        print(np.median(ucac_diff[np.isfinite(ucac_diff)]),np.median(sdss_diff[np.isfinite(sdss_diff)]))
        print(np.median(tycho_diff[np.isfinite(tycho_diff)]),np.median(gK_diff[gK_diff<10000]))
        print(np.std(ucac_diff[np.isfinite(ucac_diff)]),np.std(sdss_diff[np.isfinite(sdss_diff)]))
        print(np.std(tycho_diff[np.isfinite(tycho_diff)]),np.std(gK_diff[gK_diff<10000]))


        # In[43]:


        ucac_lim = np.median(ucac_diff[np.isfinite(ucac_diff)]) + np.std(ucac_diff[np.isfinite(ucac_diff)])#*2
        sdss_lim = np.median(sdss_diff[np.isfinite(ucac_diff)]) + np.std(sdss_diff[np.isfinite(sdss_diff)])#*2
        tycho_lim = np.median(tycho_diff[np.isfinite(ucac_diff)]) + np.std(tycho_diff[np.isfinite(tycho_diff)])#*2
        tmass_lim = np.median(gK_diff[gK_diff<10000]) + np.std(gK_diff[gK_diff<10000])*2

        print(len(np.where(ucac_diff>ucac_lim)[0]))
        print(len(np.where(sdss_diff>sdss_lim)[0]))
        print(len(np.where(tycho_diff>tycho_lim)[0]))
        print(len(np.where(gK_diff>tmass_lim)[0]))
        # print(np.where(((ucac_diff>0.86) & np.isfinite(ucac_diff))
        #                & ((sdss_diff>1.46) & np.isfinite(sdss_diff))
        #                & ((tycho_diff>0.06)  & np.isfinite(tycho_diff))
        #                & ((gK_diff>1.6) & (gK_diff<10000)))[0])
        bad_match = np.where(((ucac_diff>ucac_lim) | np.isnan(ucac_diff))
                       & ((sdss_diff>sdss_lim) | np.isnan(sdss_diff))
                       & ((tycho_diff>tycho_lim)  | np.isnan(tycho_diff))
                       & ((gK_diff>tmass_lim) | (gK_diff>10000)))[0]

        no_match = np.where(np.isnan(ucac_diff) & np.isnan(sdss_diff)
                            & np.isnan(tycho_diff) & (gK_diff>10000))[0]
        print(no_match)
        print(joint_tab["Gmag"][no_match])

        print(bad_match)
        for i in bad_match:
            if i in no_match:
        #         continue
                print("\nNo phot",joint_tab["RMAG"][i],joint_tab["RMAG_FLAG"][i])
                print("UCAC",joint_tab["UCAC_R"][i],joint_tab["UCAC_I"][i])
                print(i,joint_tab["GOLDMAN_SEQ"][i],joint_tab["HIP_ID"][i])
            else:
                print("\nUCAC",ucac_diff[i],"SDSS",sdss_diff[i])
                print("Tycho2",tycho_diff[i],"2MASS",gK_diff[i])
        #         print(i,joint_tab["GOLDMAN_SEQ"][i],joint_tab["GOLDMAN_RMAG"][i])
                print("G=",joint_tab["Gmag"][i],joint_tab["GAIA_QUAL"][i],joint_tab["Source"][i])
                print("Prot=",joint_tab["PERIOD"][i],"EPIC",joint_tab["EPIC_ID"][i])


        # In[44]:


        # but it appears some stars lack Gaia matches?
        print(np.where(joint_tab["Gmag"].mask==True)[0])
        no_match = np.where(joint_tab["Gmag"].mask==True)[0]
        print(joint_tab["RA"][no_match],joint_tab["DEC"][no_match])
        print(joint_tab["PERIOD"][no_match],joint_tab["EPIC_ID"][no_match])


        # So every rotator and possible K2 target has a Gaia match that looks reasonable. 3 stars are missing the photometry to confirm their matches, 5 do fail the photometry cross-check, and 2 stars do not have Gaia crossmatches, but none of those 10 stars will be included in the paper anyway.

    # In[45]:


    # joint_tab.dtype.names


    # In[88]:
    return joint_tab

def calc_mass_unc(mK,mK_err,dist,dist_err):
    # dmod = m-M = 5*log10(d) - 5
    # M = m-dmod = m - 5*log10(d) + 5
    # dM/dd = -5/(d*ln(10))
    dmod = 5.0*np.log10(dist) - 5
    MK = mK - dmod

    bad = np.where(mK<-99)[0]
    MK[bad] = -9999

    if type(mK_err)==np.ndarray:
        mK_err_sq = mK_err**2
    else:
        # actually a fraction, but real errors
        mK_err_sq = (mK_err*mK)**2

    dist_err_sq = (dist_err * 5.0 / (dist * np.log(10.0)))**2

    sigma_MK = np.sqrt(mK_err_sq + dist_err_sq)
    #print(sigma_MK)

    MK_low =  MK + sigma_MK #dimmer/lower mass
    MK_high = MK - sigma_MK #brighter/higher mass
    #print('low',MK_low)
    #print('high',MK_high)

    good = np.where((MK>0) & (MK_low>0) & (MK_high>0)  &
                    (MK.mask==False))[0]

    mass, mass_low, mass_high = np.ones(len(MK))*-9999.,np.ones(len(MK))*-9999.,np.ones(len(MK))*-9999.

    mass[good] = convertmass.kraus(np.asarray(MK[good]),'K','None')
    mass_low[good] = convertmass.kraus(np.asarray(MK_low[good]),'K','None')
    mass_high[good] = convertmass.kraus(np.asarray(MK_high[good]),'K','None')

    mass_errs = np.zeros(2*len(mK)).reshape((2,-1))

    mass_errs[1][good] = abs(mass_low-mass)[good]
    mass_errs[0][good] = abs(mass_high-mass)[good]

    return mass,mass_errs

def hyades_mass_unc(dist, e_dist, mK, e_mK):

    mass,mass_errs = calc_mass_unc(mK, e_mK, dist, e_dist)

    good_mass = np.where(mass>0)[0]
    good_errs = np.zeros(2*len(good_mass)).reshape((2,-1))
    good_errs[0] = mass_errs[0][good_mass]
    good_errs[1] = mass_errs[1][good_mass]

    return mass, mass_errs

def calc_gaia_masses(joint_tab):
    gdist = 1000/joint_tab["Plx"]
    e_gdist = abs(1000.0 * joint_tab["e_Plx"] /
                  (joint_tab["Plx"]**2))

    gmass, gmass_err = hyades_mass_unc(gdist, e_gdist, joint_tab["TWOMASS_K"],
                                       joint_tab["TWOMASS_KERR"])

    return gmass, gmass_err


def print_gaia(joint_tab):

    out_names = ["HYADES_IDX","RA","DEC","GOLDMAN_PMRA","GOLDMAN_PMDE",
                 "ROESER_PMRA","ROESER_PMDEC",
                 "HIP_ID","TWOMASSNAME",
                "TWOMASS_J","TWOMASS_JERR","TWOMASS_H","TWOMASS_HERR",
                 "TWOMASS_K","TWOMASS_KERR","TWOMASS_PHOT_FLAG","TWOMASS_B",
                "TWOMASS_VR","SWASP_ID","UCAC_ID",
                 "UCAC_R", "UCAC_RERR","UCAC_I","UCAC_IERR",
                "SDSS_R","SDSS_RERR","SDSS_I","SDSS_IERR","HIP_PARALLAX",
                "HIP_PAR_ERR","RPRIME","RMAG_FLAG","GOLDMAN_SEQ",
                "GOLDMAN_PLX","GOLDMAN_E_PLX","GOLDMAN_RMED","GOLDMAN_E_RMED",
                "GOLDMAN_IMED","GOLDMAN_E_IMED",
                "EPIC_ID",
                ## Don't actually include these here, make a separate table for periods
                ## "K2_PERIOD","PROSSER_PERIOD","HARTMAN_PERIOD","DELORME_LITP",
                # Now for the Gaia DR2 Columns
                "DR2Name","RA_ICRS","e_RA_ICRS","DE_ICRS",
                "e_DE_ICRS","Source","Epoch","Plx","e_Plx","pmRA","e_pmRA",
                "pmDE","e_pmDE","Gmag","e_Gmag","BPmag","e_BPmag","RPmag",
                "e_RPmag","E(BR/RP)","BP-RP","RV","e_RV",
                "RPlx","RFG","NgAL","chi2AL","RFBP","RFRP","Nper",
                # And my additions
                "KH_MASS_GAIA","e_KH_MASS_GAIA",
                "e_KH_MASS_GAIA0","e_KH_MASS_GAIA1",
                "HRD","GAIA_QUAL"]


    # In[89]:


    out_tab = joint_tab[out_names]
    out_tab.show_in_notebook()


    # In[91]:


    out_tab = joint_tab[out_names]

    out_tab["RPRIME"].mask[out_tab["RPRIME"]<-98] = True

    out_tab.rename_column("RA","RAJ2000")
    out_tab.rename_column("DEC","DEJ2000")

    out_tab.rename_column("TWOMASSNAME","2MASS")
    out_tab.rename_column("HIP_ID","HIP")
    out_tab.rename_column("SWASP_ID","SWASP")
    out_tab.rename_column("UCAC_ID","UCAC")
    out_tab.rename_column("EPIC_ID","EPIC")

    out_tab.rename_column("HIP_PARALLAX","HIP_PLX")
    out_tab.rename_column("HIP_PAR_ERR","HIP_E_PLX")

    out_tab.rename_column("TWOMASS_J","Jmag")
    out_tab.rename_column("TWOMASS_JERR","e_Jmag")
    out_tab.rename_column("TWOMASS_H","Hmag")
    out_tab.rename_column("TWOMASS_HERR","e_Hmag")
    out_tab.rename_column("TWOMASS_K","Kmag")
    out_tab.rename_column("TWOMASS_KERR","e_Kmag")
    out_tab.rename_column("TWOMASS_B","Bmag")
    out_tab.rename_column("TWOMASS_VR","VRmag")

    out_tab.rename_column("RMAG_FLAG","RPRIME_SOURCE")
    out_tab.rename_column("GOLDMAN_SEQ","[RSP2011]")

    ## Don't actually include these here, make a separate table for periods
    # out_tab.rename_column("K2_PERIOD","Prot5")
    # out_tab.rename_column("DELORME_LITP","ProtD")
    # out_tab.rename_column("PROSSER_PERIOD","ProtP")
    # out_tab.rename_column("HARTMAN_PERIOD","ProtH")
    # out_tab["Prot5"][out_tab["Prot5"]<=0] = -9999

    out_tab.rename_column("HRD","HRD_MEMBER")

    out_tab.meta = {}

    fill_values = [#(at.masked, ""),
                  ("X",""),
                  ("0",""),
                  ("-",""),
                  ("nan",""),
                  ("-9999.0000",""),
                  ("-9999.000",""),
                  ("-9999.00",""),
                  ("-9999.0",""),
                  ("-9999",""),
                  ("-99",""),
                  ("no_string",""),
                  ("no string",""),
                  ("no st",""),
                  ("HIP0","")]


    # In[123]:


    formats = {"RAJ2000":"%.6f","DEJ2000":"%.6f",
               "GOLDMAN_PMRA":"%.4f","GOLDMAN_PMDE":"%.4f",
                 "ROESER_PMRA":"%.4f","ROESER_PMDEC":"%.4f",
                "Jmag":"%.3f","e_Jmag":"%.3f",
                "Hmag":"%.3f","e_Hmag":"%.3f",
                "Kmag":"%.3f","e_Kmag":"%.3f",
                "Bmag":"%.3f","VRmag":"%.3f",
                 "UCAC_R":"%.3f", "UCAC_RERR":"%.3f",
               "UCAC_I":"%.3f","UCAC_IERR":"%.3f",
                "SDSS_R":"%.3f","SDSS_RERR":"%.3f",
               "SDSS_I":"%.3f","SDSS_IERR":"%.3f",
               "RPRIME":"%.3f",
                "HIP_PLX":"%.4f","HIP_E_PLX":"%.4f",
                "GOLDMAN_PLX":"%.4f","GOLDMAN_E_PLX":"%.4f",
               "GOLDMAN_RMED":"%.3f","GOLDMAN_E_RMED":"%.3f",
                "GOLDMAN_IMED":"%.3f","GOLDMAN_E_IMED":"%.3f",
                "Prot5":"%.2f","ProtP":"%.2f",
               "ProtH":"%.2f","ProtD":"%.2f",
                 # Now for the Gaia DR2 Columns
               "RA_ICRS":"%.6f","e_RA_ICRS":"%.6f",
               "DE_ICRS":"%.6f","e_DE_ICRS":"%.6f",
               "Plx":"%.4f","e_Plx":"%.4f",
               "pmRA":"%.4f","e_pmRA":"%.4f",
                "pmDE":"%.4f","e_pmDE":"%.4f",
               "Gmag":"%.4f","e_Gmag":"%.4f","BPmag":"%.4f",
               "e_BPmag":"%.4f","RPmag":"%.4f",
                 "e_RPmag":"%.4f","E(BR/RP)":"%.4f","BP-RP":"%.4f",
               "RV":"%.2f","e_RV":"%.2f",
                 "RPlx":"%.2f","RFG":"%.2f","chi2AL":"%.2f",
               "RFBP":"%.2f","RFRP":"%.2f",
               "KH_MASS_GAIA":"%.2f","e_KH_MASS_GAIA":"%.2f",
               "KH_MASS_GAIA0":"%.2f","e_KH_MASS_GAIA1":"%.2f"
               }


    # In[124]:


    at.write(out_tab,"Gaia_Comb_Table.csv",formats=formats,
            delimiter=",",overwrite=True,fill_values=fill_values)

def compare_jasons_matches(joint_tab):

    jdir = os.path.expanduser("~/my_papers/HyPra-JC/Hyades")
    jfile1 = os.path.join(jdir,"HyadesDR2-FirstList.txt")
    jfile2 = os.path.join(jdir,"HyadesDR2-C13.txt")

    for filename in [jfile1,jfile2]:
        print(filename.split("/")[-1])
        jtab = at.read(filename)

        jpos = SkyCoord(jtab["RA"], jtab["Dec"], unit=u.degree)
        idx, sep, _ = hpos.match_to_catalog_sky(jpos)
        print(len(idx),len(hpos),len(jpos))

        good_match = np.where(sep<(5*u.arcsec))[0]
        good_idx = idx[good_match]
        print(len(good_match),len(good_idx),len(np.unique(good_idx)))

        j_gaiaid = np.zeros_like(joint_tab["DR2Name"])
        j_gaiaid[good_match] = jtab["DR2Name"][good_idx]


        for i in good_match:
            jname = "Gaia DR2 {0}".format(j_gaiaid[i])
            if jname!=joint_tab["DR2Name"][i]:
                print("Uh oh!")
                print(joint_tab["DR2Name"][i],j_gaiaid[i])
                print(jpos[jtab["DR2Name"]==j_gaiaid[i]],"J")
                print(hpos[i],"H")
                print(joint_tab["Source"][i])


if __name__=="__main__":


    hdat,_,_,_ = cat_io.get_data("H")
    hpos = SkyCoord(hdat["RA"],hdat["DEC"],unit=u.degree)

    joint_tab = match_gaia(to_plot=False)

    gmass, gmass_errs = calc_gaia_masses(joint_tab)

    # print(np.shape(gmass_errs))
    # print(np.sum(gmass_errs,0))
    # print(np.sum(gmass_errs,1))
    # print(np.sum(gmass_errs))
    avg_errs = np.sum(gmass_errs,0)/2

    joint_tab["KH_MASS_GAIA"] = gmass
    joint_tab["e_KH_MASS_GAIA"] = avg_errs
    joint_tab["e_KH_MASS_GAIA0"] = gmass_errs[0]
    joint_tab["e_KH_MASS_GAIA1"] = gmass_errs[1]

    print_gaia(joint_tab)

    compare_jasons_matches(joint_tab)

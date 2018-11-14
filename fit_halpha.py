import os

import matplotlib
matplotlib.use("agg")
import numpy as np
import astropy.io.ascii as at
from astropy.table import Table

from PHEW import EqW
import pyspeckit as p

if __name__=="__main__":

    filelist = at.read("fixed_spectra.lst")
    files = filelist["filename"]

    dtype = [("filename","U150"),
             ("eqw_med","f8"),("eqw16","f8"), ("eqw84","f8")]
    tab = Table(np.zeros(len(files),dtype=dtype))

#    mdm_path = os.path.expanduser("~/Dropbox/data/MDM_Hyades/")
    mdm_path = os.path.expanduser("/pool/cfsi04_0/data/MDM_Hyades/")
    for i,filename in enumerate(files):
        ssplit = filename.split("/")
        name = ssplit[-1][:-12]
        print(name)
        tab["filename"][i] = filename

        if os.path.exists(filename):
#            filebase = os.path.join(mdm_path,"figures/{0}.png".format(name))
            filebase = "figures/{0}.png".format(name)
            eperc = EqW.measure_equivalent_width(str(filename),
                                6550,6580,6560,6566,
                                1000,xunit="Angstrom",to_plot=True,
                                filebase=filebase)
            tab["eqw16"][i], tab["eqw_med"][i], tab["eqw84"][i] = eperc
#        if i>3:
#            break

    # fmt = {"eqw_med":"%.4f",
    #         "eqw16":"%.4f",
    #             "eqw84":"%.4f"}
    at.write(tab,"halpha_equivalent_widths.csv",#format=fmt,
            delimiter=",",overwrite=True)

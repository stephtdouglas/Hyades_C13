import os

import astropy.io.ascii as at
import astropy.io.fits as fits
import numpy as np

from hypra.utils import cat_io

if __name__=="__main__":
    out_path = os.path.expanduser("~/Dropbox/data/MDM_Hyades/")

    filelist = at.read("new_spectra.lst")
    files = filelist["filename"]

    for file in files:
        print("\n",file)
        with fits.open(file) as hdu:
            hdr = hdu[0].header
            hdr["WAT1_001"] = hdr["WAT1_001"].replace("Angstroms", "Angstrom")
            outfilename = file.split("/")[-1].replace(".fits","_update.fits")
            outfilename = os.path.join(out_path,outfilename)
            print(outfilename)
            fits.writeto(outfilename,data=hdu[0].data,header=hdr)

    hdat,_,_,_ = cat_io.get_data("H")
    filelist1 = hdat["MDM_SPEC_FILE"][hdat["MDM_SPEC_FILE"]!="no string"]
    filelist2 = hdat["MDM_SPEC_FILE2"][hdat["MDM_SPEC_FILE2"]!="no string"]
    files = np.append(filelist1,filelist2)
    mdm_path = os.path.expanduser("~/Dropbox/data/Praesepe/MODspec/")
    for file in files:
        print("\n")
        print(os.path.join(mdm_path,file))
        with fits.open(os.path.join(mdm_path,file)) as hdu:
            hdr = hdu[0].header
            hdr["WAT1_001"] = hdr["WAT1_001"].replace("Angstroms", "Angstrom")
            outfilename = file.split("/")[-1].replace(".fits","_update.fits")
            outfilename = os.path.join(out_path,outfilename)
            print(outfilename)
            fits.writeto(outfilename,data=hdu[0].data,header=hdr)

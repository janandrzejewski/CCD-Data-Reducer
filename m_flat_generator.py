from astropy import units as u
from pathlib import Path
from astropy.nddata import CCDData
import numpy as np
from astropy.stats import mad_std
import ccdproc as ccdp
from astropy.io import fits


def inv_median(a):
    return 1 / np.median(a)


    
flat1_p = 'flat1.fits'
flat2_p = 'flat2.fits'
flat3_p = 'flat3.fits'
flat1= CCDData(fits.getdata(flat1_p), meta=fits.getheader(flat1_p),unit = 'adu')
flat2= CCDData(fits.getdata(flat2_p), meta=fits.getheader(flat2_p),unit = 'adu')
flat3= CCDData(fits.getdata(flat3_p), meta=fits.getheader(flat3_p),unit = 'adu')


to_combine = [flat1, flat2, flat3]


combined_flat = ccdp.combine(to_combine,
                                 method='median',
                                 unit=u.adu, scale=inv_median,
                                 sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                 sigma_clip_func=np.ma.median, signma_clip_dev_func=mad_std,
                                )
combined_flat.write('master_flat.fits', overwrite=True)
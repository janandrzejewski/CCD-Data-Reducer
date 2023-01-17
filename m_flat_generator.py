from astropy import units as u
import os
from pathlib import Path
from astropy.nddata import CCDData
import numpy as np
from astropy.stats import mad_std
import ccdproc as ccdp
from astropy.io import fits


def inv_median(a):
    return  1/np.median(a)
main_p = '2023_01_05'
flat_p = '2023_01_05/flat'
flat_image_type = 'flat'
m_dark_p = 'master_dark.fits'
m_bias_p = 'master_bias.fits'
ifc_flat = ccdp.ImageFileCollection(flat_p)
print("LOG:zaladowano flaty:",ifc_flat)
master_bias  = CCDData(fits.getdata(m_bias_p), meta=fits.getheader(m_bias_p),unit = 'adu')
print("LOG:otworzono m_biasa",m_bias_p)
master_dark  = CCDData(fits.getdata(m_dark_p), meta=fits.getheader(m_dark_p),unit = 'adu')
print("LOG:otworzono m_darka",m_dark_p)
master_dark = ccdp.subtract_bias(master_dark, master_bias)
print("LOG:bias usuniety z master_darka")
flat_filters = set(h['filter'] for h in ifc_flat.headers(imagetyp=flat_image_type))

dark_expo = master_dark.header['EXPTIME']
for filt in flat_filters:
    flat_list = ifc_flat.files_filtered(imagetyp = flat_image_type, filter = filt,include_path=True)
    print("LOG:W filtrze",filt, "mamy tyle flatow", len(flat_list))
    CCDDATA_flat_list = []
    for path in flat_list:
        flat = CCDData(fits.getdata(path), meta=fits.getheader(path),unit = 'adu')
        bias_subtracted = ccdp.subtract_bias(flat, master_bias)
        dark_subtracted = ccdp.subtract_dark(bias_subtracted,master_dark,                                      
                                            exposure_time='EXPTIME',
                                            exposure_unit=u.adu,
                                            scale=True)
        #verscan_subtracted = ccdp.subtract_overscan(dark_subtracted, overscan=dark_subtracted[:, -1:], median=True)
        CCDDATA_flat_list.append(dark_subtracted)
    combined_flat = ccdp.combine(CCDDATA_flat_list,
                                     method='median',
                                     unit=u.adu, scale=inv_median,
                                     sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                     sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,maxiters=2
                                    )

    combined_flat.meta['combined'] = True
    flat_file_name = 'master_flat_{}.fits'.format(filt)
    combined_flat.write(os.path.join(main_p, flat_file_name), overwrite=True)

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

i=0
main_p = '2023_01_05'
flat_p = '2023_01_05/flat'
flat_image_type = 'flat'
m_dark_p = 'master_dark.fit'
ifc_flat = ccdp.ImageFileCollection(flat_p)


master_dark  = CCDData(fits.getdata(m_dark_p), meta=fits.getheader(m_dark_p),unit = 'adu')
flat_filters = set(h['filter'] for h in ifc_flat.headers(imagetyp=flat_image_type))
flat_expo = list(h['exposure'] for h in ifc_flat.headers(imagetyp=flat_image_type))
for filt in flat_filters:
    flat_list = ifc_flat.files_filtered(imagetyp = flat_image_type, filter = filt,include_path=True)
    CCDDATA_flat_list = [CCDData(fits.getdata(path), meta=fits.getheader(path),unit = 'adu')for path in flat_list]
    
    combined_flat = ccdp.combine(CCDDATA_flat_list,
                                     method='median',
                                     unit=u.adu, scale=inv_median,
                                     sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                     sigma_clip_func=np.ma.median, signma_clip_dev_func=mad_std,
                                    )
    combined_flat.meta['combined'] = True
    i=i+1
    flat_file_name = 'master_flat_{}.fits'.format(filt.replace("''", "p"))
    combined_flat.write(os.path.join(main_p, flat_file_name), overwrite=True)
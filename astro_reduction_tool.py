from astropy import units as u
import os
from astropy.nddata import CCDData
import numpy as np
from astropy.stats import mad_std
import ccdproc as ccdp
from astropy.io import fits
import logging
from datetime import datetime
import argparse


logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
logger.addHandler(console_handler)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

def read_path_into_CCDData(data_path, keys, _filter=None):
    ccd_list = []

    if not os.path.isdir(data_path): # better to add pattern  glob_include and ext parameters
        logger.error(f"No directory {data_path}")
        raise ValueError(f"No directory {data_path}")
    ccd_collection = ccdp.ImageFileCollection(data_path, keywords=keys)
    if _filter:
        ccd_collection = ccd_collection.files_filtered(filter=_filter)
    if len(ccd_collection) != 0:
        logger.info(f"Loading data from: {data_path}")
        for ccd_path in ccd_collection:
            dir = os.path.join(data_path, ccd_path)
            
            ccd_file = CCDData(
                fits.getdata(dir), meta=fits.getheader(dir, keywords=keys), unit="adu"
                )

            ccd_list.append(ccd_file)   
            now = datetime.now()
            formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
            logging.info(
                "%s  Loading %s file %s (1 of 1)",
                formatted_time,
                ccd_file.header["IMAGETYP"],
                ccd_path,
            )

    return ccd_list 


def debias_dedark(file_list, master_dark, master_bias):
    if isinstance(file_list, list):
        for h in range(0, len(file_list)):
            bias_subtracted = ccdp.subtract_bias(file_list[h], master_bias)
            file_list[h] = ccdp.subtract_dark(
                bias_subtracted,
                master_dark,
                exposure_time="EXPTIME",
                exposure_unit=u.adu,
                scale=True,
            )
    else:
        bias_subtracted = ccdp.subtract_bias(file_list, master_bias)
        file_list = ccdp.subtract_dark(
            bias_subtracted,
            master_dark,
            exposure_time="EXPTIME",
            exposure_unit=u.adu,
            scale=True,
        )
    return file_list


def inv_median(a):
    return 1 / np.median(a)


def master_flat_generator(flat_p, master_dark, master_bias, filt=""):
    m_flat = []
    for f in filt:
        keys = ["imagetyp", "filter"]
        flat_list = read_path_into_CCDData(flat_p, keys, f)
        flat_subtracted_list = debias_dedark(flat_list, master_dark, master_bias)
        combined_flat = ccdp.combine(
            flat_subtracted_list,
            method="median",
            unit=u.adu,
            scale=inv_median,
            sigma_clip=True,
            sigma_clip_low_thresh=5.0,
            sigma_clip_high_thresh=5.0,
            sigma_clip_func=np.ma.median,
            sigma_clip_dev_func=mad_std,
            flat=2,
        )
        combined_flat.meta["combined"] = True
        now = datetime.now()
        formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
        logging.info("%s  Created master_flat in %s filtr", formatted_time, f)
        m_flat.append(combined_flat)

    return m_flat


def folders_reduction(dir_folder, m_dark, m_bias, m_flat_list, filt, d):
    for f in filt:
        
        ccd_folder = read_path_into_CCDData(
            os.path.join(dir_folder, d), ["imagetyp", "filter"], f
        )
        if len(ccd_folder) != 0:
            for ccd_file in ccd_folder:
                bias_subtracted = ccdp.subtract_bias(ccd_file, m_bias)
                dark_subtracted = ccdp.subtract_dark(
                    bias_subtracted,
                    m_dark,
                    exposure_time="EXPTIME",
                    exposure_unit=u.adu,
                    scale=True,
                )
                ccd_file = ccdp.flat_correct(dark_subtracted, m_flat_list[f])
            now = datetime.now()
            formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
            logging.info("%s  Subtracted %s folder in %s filtr", formatted_time, d, f)
    return ccd_folder


if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Reduce astro images")

    # Add arguments
    parser.add_argument(
        "-d",
        "--dir",
        type=str,
        help="directory containing FITS images or a PHOT filepath",
        required=True,
    )
    parser.add_argument(
        "-md",
        "--master_dark",
        type=str,
        default="master_dark.fits",
        help="The Master dark",
    )
    parser.add_argument(
        "-mb",
        "--master_bias",
        type=str,
        default="master_bias.fits",
        help="The Master bias",
    )
    args = parser.parse_args()
    folders = [
        f
        for f in os.listdir(args.dir)
        if os.path.isdir(os.path.join(args.dir, f)) and f != "flat"
    ]
    m_dark = CCDData(
            fits.getdata(args.master_dark), meta=fits.getheader(args.master_dark), unit="adu"
            )
    m_bias = CCDData(
            fits.getdata(args.master_bias), meta=fits.getheader(args.master_bias), unit="adu"
            )
    flat_image_type = "flat"
    m_dark = ccdp.subtract_bias(m_dark, m_bias)
    flat_p = os.path.join(args.dir + "flat")
    ifc_flat = ccdp.ImageFileCollection(flat_p, keywords=["imagetyp", "filter"])
    flat_filters = set(h["filter"] for h in ifc_flat.headers(imagetyp=flat_image_type))
    m_flat_list = master_flat_generator(flat_p, m_dark, m_bias, flat_filters)
    m_flat_dictionary= {"gs":[],"rs":[],"in":[],"Luminance":[],"B":[],"V":[],"R":[]};
    for ccd_file in m_flat_list:
        filter = ccd_file.header['filter']
        m_flat_dictionary[filter]= ccd_file
    # for d in folders:
    #     folders_reduction(
    #         args.dir, m_dark, m_bias, m_flat_list, flat_filters, d
    #     )
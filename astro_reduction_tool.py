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


def read_path_into_CCDData(data_p, keys, filt=""):
    ccd_list = []
    if os.path.isdir(data_p):
        logger.info("Loading data from: %s", data_p)
        ccd_p = ccdp.ImageFileCollection(data_p, keywords=keys)
        if filt != "":
            ccd_p = ccd_p.files_filtered(filter=filt)
        for h in range(0, len(ccd_p)):
            dir = os.path.join(data_p, ccd_p[h])
            now = datetime.now()
            formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
            ccd_file = CCDData(
                fits.getdata(dir), meta=fits.getheader(dir, keywords=keys), unit="adu"
            )
            ccd_list.append(ccd_file)
            logger.info(
                "%s  Loading %s file %s (%d of %d)",
                formatted_time,
                ccd_file.header["IMAGETYP"],
                dir,
                h + 1,
                len(ccd_p),
            )
        return ccd_list
    elif os.path.isfile(data_p):
        now = datetime.now()
        formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
        ccd_file = CCDData(
            fits.getdata(data_p), meta=fits.getheader(data_p, keywords=keys), unit="adu"
        )
        logging.info(
            "%s  Loading %s file %s (1 of 1)",
            formatted_time,
            ccd_file.header["IMAGETYP"],
            data_p,
        )
        return ccd_file
    else:
        logging.error(
            "Wrong path. %s should be dir with .fits files or a path to a single .fits file",
            data_p,
        )


def debias_dedark(file_list, master_dark, master_bias):
    # dark_expo = master_dark.header["EXPTIME"]
    # if isinstance(file_list,list):
    #     print(file_list[1])
    #     for h in range(0, len(file_list)):
    #         ccd_file = file_list[h]
    #         bias_subtracted = ccdp.subtract_bias(ccd_file, master_bias)
    #         file_list[h] = ccdp.subtract_dark(
    #             bias_subtracted,
    #             master_dark,
    #             exposure_time=dark_expo,
    #             exposure_unit=u.adu,
    #             scale=True,
    #         )
    # else:
    #     bias_subtracted = ccdp.subtract_bias(file_list, master_bias)
    #     file_list = ccdp.subtract_dark(
    #         bias_subtracted,
    #         master_dark,
    #         exposure_time=dark_expo,
    #         exposure_unit=u.adu,
    #         scale=True,
    #     )
    return file_list


def inv_median(a):
    return 1 / np.median(a)


def master_flat_generator(flat_p, master_dark, master_bias, filt=""):
    m_flat = {}
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
        logging.info(
            "%s  Created master_flat in %s filtr",
            formatted_time,
            f
        )
        m_flat[f] = combined_flat

    return m_flat


def folders_reduction(dir_folder, m_dark, m_bias, m_flat_list, filt, folders):

    for h in range(0, len(folders)):
        for f in filt:
            print(os.path.join(dir_folder, folders[h]))
            ccd_folder = read_path_into_CCDData(
                os.path.join(dir_folder, folders[h]), ["imagetyp"], f
            )
            for x in range(0, len(ccd_folder)):
                ccd_file = ccd_folder[x]
                # ccd_file = ccdp.subtract_bias(ccd_file, m_bias)
                # ccd_file = ccdp.subtract_dark(
                #     ccd_file,
                #     m_dark,
                #     exposure_time="EXPTIME",
                #     exposure_unit=u.adu,
                #     scale=True,
                # )
                ccd_file = ccdp.flat_correct(ccd_file, m_flat_list[f])
                ccd_folder[x] = ccd_file
        now = datetime.now()
        formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
        logging.info(
            "%s  Subtracted %s folder in %s filtr",
            formatted_time,
            folders[h],
            f
        )

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
    m_dark = read_path_into_CCDData(args.master_dark, ["imagetyp","filter"])
    m_bias = read_path_into_CCDData(args.master_bias, ["imagetyp","filter"])
    flat_image_type = "flat"
    m_dark = ccdp.subtract_bias(m_dark, m_bias)
    flat_p = os.path.join(args.dir + "flat")
    ifc_flat = ccdp.ImageFileCollection(flat_p, keywords=["imagetyp", "filter"])
    flat_filters = set(h["filter"] for h in ifc_flat.headers(imagetyp=flat_image_type))
    m_flat_list = master_flat_generator(flat_p, m_dark, m_bias, flat_filters)
    for d in folders:
        folders_reduction(args.dir, args.master_dark, args.master_bias, m_flat_list,flat_filters,folders)

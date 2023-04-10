import argparse
import configparser
import logging
import os
from datetime import datetime

import ccdproc as ccdp
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import mad_std

logging.basicConfig(level=logging.INFO, handlers=[logging.StreamHandler()])
logger = logging.getLogger()

config = configparser.ConfigParser()
config.read("config.ini")

flat_image_type = config["data"]["flat_image_type"]
flat_dir = config["data"]["flat_dir"]
flat_keywords = config["processing"]["flat_keywords"]
combine_method = config["processing"]["combine_method"]
sigma_clip_low_thresh_config = float(config["processing"]["sigma_clip_low_thresh"])
sigma_clip_high_thresh_config = float(config["processing"]["sigma_clip_high_thresh"])


class ccdData:
    def __init__(self, data_list, name_list):
        self.data_list = data_list
        self.name_list = name_list

    def __str__(self):
        return "List 1: {}\nList 2: {}".format(self.data_list, self.name_list)


def inv_median(a):
    return 1 / np.median(a)


def log(message, *args):
    now = datetime.now()
    formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
    logging.info("%s  " + message, formatted_time, *args)


def read_path_into_CCDData(data_path, keys, _filter=None, ifc_master_flat=None):
    if ifc_master_flat is None:
        ifc_master_flat = []
    data_list = []
    name_list = []
    if not os.path.isdir(data_path):
        logger.error(f"No directory {data_path}")
        raise ValueError(f"No directory {data_path}")
    if len(ifc_master_flat) == 0:
        ccd_collection = ccdp.ImageFileCollection(data_path, keywords=keys)
        if _filter:
            ccd_collection = ccd_collection.files_filtered(filter=_filter)
        else:
            ccd_collection = ccd_collection.files
    else:
        ccd_collection = ifc_master_flat
    max = len(ccd_collection)
    if max != 0:
        logger.info(f"Loading data from: {data_path}")
        for ccd_path_idx, ccd_path in enumerate(ccd_collection):
            dir = os.path.join(data_path, ccd_path)
            ccd_file = CCDData(
                fits.getdata(dir), meta=fits.getheader(dir, keywords=keys), unit="adu"
            )
            header_value = ccd_file.meta.get("combined")
            if header_value is None or len(ifc_master_flat) != 0:
                data_list.append(ccd_file)
                name_list.append(ccd_path)
                log(
                    "Loading %s file %s (%d of %d)",
                    ccd_file.header["IMAGETYP"],
                    ccd_path,
                    ccd_path_idx + 1,
                    max,
                )
    ccd_data_list = ccdData(data_list, name_list)
    return ccd_data_list


def debias_dedark(file_list, master_dark, master_bias):
    for ccd_file_idx, ccd_file in enumerate(file_list.data_list):
        ccd_name = file_list.name_list[ccd_file_idx]
        bias_subtracted = ccdp.subtract_bias(ccd_file, master_bias)
        file_list.data_list[ccd_file_idx] = ccdp.subtract_dark(
            bias_subtracted,
            master_dark,
            exposure_time="EXPTIME",
            exposure_unit=u.adu,
            scale=True,
        )
        log(
            "Debias and dedark %s file %s (%d of %d)",
            ccd_file.header["IMAGETYP"],
            ccd_name,
            ccd_file_idx + 1,
            len(file_list.data_list),
        )
    return file_list.data_list


def master_flat_generator(flat_p, master_dark, master_bias, filt):
    m_flat = []
    m_flat_n = []
    keys = flat_keywords
    for f in filt:
        flat_data = read_path_into_CCDData(flat_p, keys, f)
        flat_subtracted_list = debias_dedark(flat_data, master_dark, master_bias)
        combined_flat = ccdp.combine(
            flat_subtracted_list,
            method=combine_method,
            unit=u.adu,
            scale=inv_median,
            sigma_clip=True,
            sigma_clip_low_thresh=sigma_clip_low_thresh_config,
            sigma_clip_high_thresh=sigma_clip_high_thresh_config,
            sigma_clip_func=np.ma.median,
            sigma_clip_dev_func=mad_std,
            flat=2,
        )
        combined_flat.meta["combined"] = True
        name = os.path.join(flat_p, f"master_flat_{f}.fits")
        combined_flat.write(name, overwrite=True)
        log("Created master_flat in %s filtr", f)
        m_flat.append(combined_flat)
        m_flat_n.append(name)
        m_flat_list = ccdData(m_flat, m_flat_n)
    return m_flat_list


def folders_reduction(dir_folder, m_dark, m_bias, m_flat_dictionary, folder):
    path = os.path.join(dir_folder, folder)
    ifc_folder = ccdp.ImageFileCollection(path, keywords=flat_keywords)
    filt = set(h["filter"] for h in ifc_folder.headers())
    out_dir = os.path.join(path, "pipeline_out")
    os.makedirs(out_dir, exist_ok=True)
    for f in filt:
        ccd_folder = read_path_into_CCDData(path, flat_keywords, f)
        max = len(ccd_folder.data_list)
        if max != 0:
            ccd_file_idx = 0
            for ccd_file_idx, ccd_file in enumerate(ccd_folder.data_list):
                ccd_name = ccd_folder.name_list[ccd_file_idx]
                bias_subtracted = ccdp.subtract_bias(ccd_file, m_bias)
                dark_subtracted = ccdp.subtract_dark(
                    bias_subtracted,
                    m_dark,
                    exposure_time="EXPTIME",
                    exposure_unit=u.adu,
                    scale=True,
                )
                ccd_file = ccdp.flat_correct(
                    dark_subtracted, m_flat_dictionary[f], min_value=0.01
                )
                ccd_file.meta["HISTORY"] = "Bias corrected"
                ccd_file.meta["HISTORY"] = "Dark corrected"
                ccd_file.meta["HISTORY"] = "Flat corrected"
                ccd_file.data = ccd_file.data.astype(np.uint16)
                ccd_file.write(os.path.join(out_dir, ccd_name), overwrite=True)
                log(
                    "Saved %s (%d in %d)",
                    os.path.join(out_dir, ccd_name),
                    ccd_file_idx + 1,
                    max,
                )
            log("Subtracted %s folder in %s filtr", folder, f)


def astro_reduction(args):
    folders = [
        f
        for f in os.listdir(args.dir)
        if os.path.isdir(os.path.join(args.dir, f)) and f != "flat"
    ]
    m_dark = CCDData(
        fits.getdata(args.master_dark),
        meta=fits.getheader(args.master_dark),
        unit="adu",
    )
    m_bias = CCDData(
        fits.getdata(args.master_bias),
        meta=fits.getheader(args.master_bias),
        unit="adu",
    )
    flat_image_type = "flat"
    flat_p = os.path.join(args.dir, flat_dir)
    if args.master_flat is None:
        ifc_flat = ccdp.ImageFileCollection(flat_p, keywords=flat_keywords)
        flat_filters = set(
            h["filter"] for h in ifc_flat.headers(imagetyp=flat_image_type)
        )
        m_flat_list = master_flat_generator(flat_p, m_dark, m_bias, flat_filters)
    else:
        log("Loading master flats from %s", args.master_flat)
        ifc_flat = ccdp.ImageFileCollection(
            args.master_flat, keywords=["imagetyp", "filter", "combined"]
        )
        m_flat_list = read_path_into_CCDData(
            args.master_flat,
            ["imagetyp", "filter", "combined"],
            None,
            ifc_flat.files_filtered(combined=True),
        )
    m_flat_dictionary = {}
    for f in flat_filters:
        m_flat_dictionary[f] = None
    for ccd_file in m_flat_list.data_list:
        filter = ccd_file.header["filter"]
        m_flat_dictionary[filter] = ccd_file
    for d in folders:
        folders_reduction(args.dir, m_dark, m_bias, m_flat_dictionary, d)


def parse_args():
    parser = argparse.ArgumentParser(description="Reduce astro images")
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
    parser.add_argument(
        "-mf",
        "--master_flat",
        type=str,
        default=None,
        help="The master flat",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    astro_reduction(args)


if __name__ == "__main__":
    main()

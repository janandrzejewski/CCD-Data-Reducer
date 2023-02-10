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


class ccd_data:
    def __init__(self, data_list, name_list):
        self.data_list = data_list
        self.name_list = name_list

    def __str__(self):
        return "List 1: {}\nList 2: {}".format(self.data_list, self.name_list)

def inv_median(a):
    return 1 / np.median(a)

def read_path_into_CCDData(data_path, keys, _filter=None, _ifc_mf = [] ):
    data_list = []
    name_list = []
    if not os.path.isdir(
        data_path
    ):
        logger.error(f"No directory {data_path}")
        raise ValueError(f"No directory {data_path}")
    if len(_ifc_mf) == 0:
        ccd_collection = ccdp.ImageFileCollection(data_path, keywords=keys)
        if _filter:
            ccd_collection = ccd_collection.files_filtered(filter=_filter)
        else:
            ccd_collection = ccd_collection.files
    else:
        ccd_collection = _ifc_mf
    if len(ccd_collection) != 0:
        logger.info(f"Loading data from: {data_path}")
        for ccd_path in ccd_collection:
            dir = os.path.join(data_path, ccd_path)
            ccd_file = CCDData(
                fits.getdata(dir), meta=fits.getheader(dir, keywords=keys), unit="adu"
            )
            header_value = ccd_file.meta.get("combined")
            if header_value is None or len(_ifc_mf) != 0:
                data_list.append(ccd_file)
                name_list.append(ccd_path)
                now = datetime.now()
                formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
                logging.info(
                    "%s  Loading %s file %s (1 of 1)",
                    formatted_time,
                    ccd_file.header["IMAGETYP"],
                    ccd_path,
                )
    ccd_data_list = ccd_data(data_list,name_list)
    return ccd_data_list

def debias_dedark(file_list, master_dark, master_bias):
    for h in range(0, len(file_list.data_list)):
        ccd_file = file_list.data_list[h]
        ccd_name = file_list.name_list[h]
        bias_subtracted = ccdp.subtract_bias(ccd_file, master_bias)
        file_list.data_list[h] = ccdp.subtract_dark(
            bias_subtracted,
            master_dark,
            exposure_time="EXPTIME",
            exposure_unit=u.adu,
            scale=True,
        )
        now = datetime.now()
        formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
        logging.info(
            "%s  Debias and dedark %s file %s (%d of %d)",
            formatted_time,
            ccd_file.header["IMAGETYP"],
            ccd_name,
            h+1,
            len(file_list.data_list)
        )
    return file_list.data_list


def master_flat_generator(flat_p, master_dark, master_bias, filt):
    m_flat = []
    m_flat_n = []
    keys = ["imagetyp", "filter"]
    for f in filt:
        flat_data = read_path_into_CCDData(flat_p, keys, f)
        flat_subtracted_list = debias_dedark(flat_data, master_dark, master_bias)
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
        name = os.path.join(flat_p,"master_flat_{}.fits".format(f))
        combined_flat.write(name, overwrite=True)
        now = datetime.now()
        formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
        logging.info("%s  Created master_flat in %s filtr", formatted_time, f)
        m_flat.append(combined_flat)
        m_flat_n.append(name)
        m_flat_list = ccd_data(m_flat,m_flat_n)
    return m_flat_list


def folders_reduction(dir_folder, m_dark, m_bias, m_flat_dictionary, d):
    path = os.path.join(dir_folder,d)
    ifc_folder = ccdp.ImageFileCollection(path, keywords=["imagetyp", "filter"])
    filt = set(h["filter"] for h in ifc_folder.headers())
    if os.path.exists(os.path.join(path, "pipeline_out")) is not True:
        os.mkdir(os.path.join(path, "pipeline_out"))
    for f in filt:
        ccd_folder = read_path_into_CCDData(
            path, ["imagetyp", "filter"], f
        )
        max = len(ccd_folder.data_list)
        if max != 0:
            x=0
            for x in range(0,len(ccd_folder.data_list)):
                ccd_name = ccd_folder.name_list[x]
                ccd_file = ccd_folder.data_list[x]
                bias_subtracted = ccdp.subtract_bias(ccd_file, m_bias)
                dark_subtracted = ccdp.subtract_dark(
                    bias_subtracted,
                    m_dark,
                    exposure_time="EXPTIME",
                    exposure_unit=u.adu,
                    scale=True,
                )
                ccd_file = ccdp.flat_correct(dark_subtracted, m_flat_dictionary[f],min_value=0.01)
                ccd_file.write(os.path.join(path,"pipeline_out", ccd_name), overwrite=True)
                now = datetime.now()
                formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
                logging.info("%s  Saved %s (%d in %d)", formatted_time, os.path.join(path,"pipeline_out", ccd_name),x,max)
            now = datetime.now()
            formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
            logging.info("%s  Subtracted %s folder in %s filtr", formatted_time, d, f)


if __name__ == "__main__":
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
    flat_p = os.path.join(args.dir, "flat")
    if args.master_flat is None:
        ifc_flat = ccdp.ImageFileCollection(flat_p, keywords=["imagetyp", "filter"])
        flat_filters = set(h["filter"] for h in ifc_flat.headers(imagetyp=flat_image_type))
        m_flat_list = master_flat_generator(flat_p, m_dark, m_bias, flat_filters)
    else:
        now = datetime.now()
        formatted_time = now.strftime("[%Y-%m-%dT%H:%M:%S.%f]")
        logging.info(
            "%s  Loading master flats from %s",
            formatted_time,
            args.master_flat
        )
        ifc_flat = ccdp.ImageFileCollection(args.master_flat, keywords=["imagetyp", "filter","combined"])
        m_flat_list = read_path_into_CCDData(args.master_flat,["imagetyp", "filter","combined"], None ,ifc_flat.files_filtered(combined = True))
    m_flat_dictionary= {"gs":[],"rs":[],"in":[],"Luminance":[],"B":[],"V":[],"R":[]};
    for ccd_file in m_flat_list.data_list:
        filter = ccd_file.header['filter']
        m_flat_dictionary[filter]= ccd_file
    for d in folders:
        folders_reduction(
            args.dir, m_dark, m_bias, m_flat_dictionary, d
        )
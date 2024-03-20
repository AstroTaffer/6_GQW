import json
from os import listdir, makedirs
from os.path import exists

import astropy.io.fits as fits


# region FITS files
def get_fits_file_paths_from_dir(dir_path):
    ff_list = listdir(dir_path)
    _ = 0
    while _ < len(ff_list):
        if not ff_list[_].count(".fits"):
            ff_list.pop(_)
            continue
        ff_list[_] = dir_path + ff_list[_]
        _ += 1
    print(f"{len(ff_list)} fits files found in {dir_path}\n")
    return ff_list


def read_fits_file(file_path):
    with fits.open(file_path) as hdu_list:
        hdu_list.verify("fix")
        header = hdu_list[0].header
        data = hdu_list[0].data
    return header, data


def check_out_directory(out_dir):
    if not exists(out_dir):
        makedirs(out_dir)
# endregion


# region configs
def read_json_config(file_name):
    with open(file_name, "r") as config_file:
        return json.load(config_file)
# endregion

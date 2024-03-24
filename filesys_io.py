import json
from os import makedirs
from os.path import exists

import astropy.io.fits as fits


def read_fits_file(ff_name):
    with fits.open(ff_name) as hdu_list:
        hdu_list.verify("fix")
        header = hdu_list[0].header
        data = hdu_list[0].data
    return header, data


def read_json_config(file_name):
    with open(file_name, 'r') as config_file:
        return json.load(config_file)


def check_out_dir(out_dir):
    if not exists(out_dir):
        makedirs(out_dir)

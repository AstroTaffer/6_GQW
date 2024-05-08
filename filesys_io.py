import json
from os import makedirs
from os.path import exists

import astropy.io.fits as fits
from astropy.io.ascii import masked
from astropy.table import Table
import numpy as np


def read_fits_file(ff_name):
    with fits.open(ff_name) as hdu_list:
        hdu_list.verify('fix')
        header = hdu_list[0].header
        data = hdu_list[0].data
    return header, data


def read_json_config(file_name):
    with open(file_name, 'r') as config_file:
        config = json.load(config_file)
        check_out_dir(f"{config['OUT_DIR']}r\\")
        check_out_dir(f"{config['OUT_DIR']}i\\")
        return config


def check_out_dir(out_dir):
    if not exists(out_dir):
        makedirs(out_dir)


def write_ff_list(ff_list, out_dir, prefix='unknw'):
    with open(f"{out_dir}{prefix}_ff_list.txt", 'w') as ff_file:
        ff_file.write('\n'.join(ff_list))


def read_ff_list(out_dir, prefix='unknw'):
    with open(f"{out_dir}{prefix}_ff_list.txt", 'r') as ff_file:
        return ff_file.read().splitlines()


def write_sel_src_cat(src_cat, out_dir, prefix='unknw'):
    src_cat.write(f"{out_dir}{prefix}_src_cat.txt", overwrite=True, delimiter='\t',
                  format='ascii.commented_header', fill_values=[(masked, 'nan')])


def read_sel_src_cat(out_dir, prefix='unknw'):
    return Table.read(f"{out_dir}{prefix}_src_cat.txt", delimiter='\t', format='ascii.commented_header',
                      fill_values=[(masked, 'nan')])


def write_phot_res(flux, magn, merr, out_dir, prefix='unknw'):
    if flux is not None:
        np.save(f"{out_dir}{prefix}_flux.npy", flux)
    if magn is not None:
        np.save(f"{out_dir}{prefix}_magn.npy", magn)
    if merr is not None:
        np.save(f"{out_dir}{prefix}_merr.npy", merr)


def read_phot_res(out_dir, prefix='unknw'):
    try:
        flux = np.load(f"{out_dir}{prefix}_flux.npy")
    except OSError:
        flux = None
    try:
        magn = np.load(f"{out_dir}{prefix}_magn.npy")
    except OSError:
        magn = None
    try:
        merr = np.load(f"{out_dir}{prefix}_merr.npy")
    except OSError:
        merr = None
    return flux, magn, merr

from os import listdir

from filesys_io import *


def rm_select_fits_files(cfg):
    # ff_list = _list_ff_in_dir(cfg['IN_FITS_DIR'])
    # ff_list = _sel_ff_with_wcs(ff_list)
    # print(f"\n{len(ff_list)} applicable fits files found in {cfg['IN_FITS_DIR']}")

    # write_ff_list(ff_list, cfg['OUT_DIR'])

    ff_list = read_ff_list(cfg['OUT_DIR'])

    return ff_list


def _list_ff_in_dir(dir_path):
    ff_list = listdir(dir_path)
    _ = 0
    while _ < len(ff_list):
        if not ff_list[_].count(".fits"):
            ff_list.pop(_)
            continue
        ff_list[_] = dir_path + ff_list[_]
        _ += 1
    return ff_list


def _sel_ff_with_wcs(ff_list):
    _ = 0
    while _ < len(ff_list):
        header = read_fits_file(ff_list[_])[0]
        if 'WCSAXES' not in header:
            ff_list.pop(_)
            continue
        _ += 1
    return ff_list

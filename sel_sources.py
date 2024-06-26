import numpy as np
from astropy.wcs import WCS
from astropy.table import Table
from regions import PixCoord

from filesys_io import read_fits_file


def rm_select_sources(ff_list, cfg, flt_cname, mag_lim):
    cat = _read_hc_from_file(cfg['HC_CAT_PATH'], flt_cname, mag_lim)
    cat = _sel_hc_sources(cat, ff_list[0], cfg['IMAGE_EDGE'], cfg['SRC_PROX_LIM'])
    print(f"{len(cat)} sources selected in {flt_cname}")

    return cat


def _read_hc_from_file(cat_path, cat_filter, cat_mag_lim):
    cat = Table.read(cat_path, format='ascii.fixed_width_no_header',
                     names=('ID', 'RAh', 'RAm', 'RAs', 'DE-', 'DEd', 'DEm', 'DEs',
                            'rmag', 'imag', 'e_rmag', 'e_imag'),
                     col_starts=(0, 7, 10, 13, 19, 20, 23, 26, 38, 45, 58, 64),
                     col_ends=(5, 8, 11, 17, 19, 21, 24, 29, 43, 50, 62, 68))

    cat.remove_rows(np.where(cat[cat_filter] > cat_mag_lim))

    cat.add_column(np.zeros(len(cat)), index=1, name='RAJ2000')
    cat['RAJ2000'] = (cat['RAh'] + cat['RAm'] / 60 + cat['RAs'] / 3600) * 15
    cat.remove_columns(['RAh', 'RAm', 'RAs'])

    cat.add_column(np.zeros(len(cat)), index=2, name='DEJ2000')
    cat['DEJ2000'] = cat['DEd'] + cat['DEm'] / 60 + cat['DEs'] / 3600
    cat['DEJ2000'][np.where(cat['DE-'] == '-')] *= -1
    cat.remove_columns(['DE-', 'DEd', 'DEm', 'DEs'])

    # ID              Object identification number
    # RAJ2000  [deg]  Right Ascension (J2000)
    # DEJ2000  [deg]  Declination (J2000)
    # rmag     [mag]  The r band magnitude
    # imag     [mag]  The i band magnitude
    # e_rmag   [mag]  The 1{sigma} error in rmag
    # e_imag   [mag]  The 1{sigma} error in imag
    return cat


def _sel_hc_sources(cat, ff_name, img_edge, src_prox_lim):
    ff_header = read_fits_file(ff_name)[0]
    ff_wcs = WCS(ff_header)

    src_xy_coords = ff_wcs.all_world2pix(cat['RAJ2000'], cat['DEJ2000'], 0)
    bad_src_mask = ((src_xy_coords[0] < img_edge) |
                    (src_xy_coords[0] > ff_header['NAXIS1'] - img_edge) |
                    (src_xy_coords[1] < img_edge) |
                    (src_xy_coords[1] > ff_header['NAXIS2'] - img_edge))
    cat.remove_rows(bad_src_mask)
    src_xy_coords[0] = src_xy_coords[0][~bad_src_mask]
    src_xy_coords[1] = src_xy_coords[1][~bad_src_mask]

    src_num = len(cat)
    src_pc = PixCoord(src_xy_coords[0], src_xy_coords[1])
    src_is_close = np.zeros((src_num, src_num), dtype=bool)

    for target_id in range(src_num):
        src_sep = src_pc[target_id].separation(src_pc)
        src_is_close[target_id] = src_sep <= src_prox_lim
        src_is_close[target_id, target_id] = False
    src_is_close = src_is_close.any(axis=0)
    cat.remove_rows(src_is_close)

    return cat

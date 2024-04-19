import warnings

from astropy.wcs import WCS

from filesys_io import *


def rm_select_sources(ff_list, cfg):
    catalog = _read_hc_from_file(cfg['HC_CAT_PATH'], cfg['MAG_LIMIT'])
    catalog = _sel_hc_sources_by_image_edges(catalog, ff_list[0], cfg['IMAGE_EDGE'])
    print(f"{len(catalog)} applicable sources found in {cfg['HC_CAT_PATH']}\n")
    # write_sel_src_cat(catalog, cfg['OUT_DIR'])

    # catalog = read_sel_src_cat(cfg['OUT_DIR'])
    return catalog


def _read_hc_from_file(cat_path, mag_limit):
    catalog = Table.read(cat_path, format='ascii.fixed_width_no_header',
                         names=('ID', 'RAh', 'RAm', 'RAs', 'DE-', 'DEd', 'DEm', 'DEs',
                                'rmag', 'imag', 'e_rmag', 'e_imag'),
                         col_starts=(0, 7, 10, 13, 19, 20, 23, 26, 38, 45, 58, 64),
                         col_ends=(5, 8, 11, 17, 19, 21, 24, 29, 43, 50, 62, 68))

    catalog.remove_rows(np.where((catalog['rmag'] >= mag_limit) | (catalog['imag'] >= mag_limit)))

    catalog.add_column(np.zeros(len(catalog)), index=1, name='RAJ2000')
    catalog['RAJ2000'] = (catalog['RAh'] + catalog['RAm'] / 60 + catalog['RAs'] / 3600) * 15
    catalog.remove_columns(['RAh', 'RAm', 'RAs'])

    catalog.add_column(np.zeros(len(catalog)), index=2, name='DEJ2000')
    catalog['DEJ2000'] = catalog['DEd'] + catalog['DEm'] / 60 + catalog['DEs'] / 3600
    catalog['DEJ2000'][np.where(catalog['DE-'] == '-')] *= -1
    catalog.remove_columns(['DE-', 'DEd', 'DEm', 'DEs'])

    # ID              Object identification number
    # RAJ2000  [deg]  Right Ascension (J2000)
    # DEJ2000  [deg]  Declination (J2000)
    # rmag     [mag]  The r band magnitude
    # imag     [mag]  The i band magnitude
    # e_rmag   [mag]  The 1{sigma} error in rmag
    # e_imag   [mag]  The 1{sigma} error in imag
    return catalog


def _sel_hc_sources_by_image_edges(catalog, ff_name, image_edge):
    warnings.simplefilter("ignore")
    ff_header = read_fits_file(ff_name)[0]
    ff_wcs = WCS(ff_header)

    sources_xy_coords = ff_wcs.all_world2pix(catalog['RAJ2000'], catalog['DEJ2000'], 0)
    bad_sources_mask = np.where((sources_xy_coords[0] < image_edge) |
                                (sources_xy_coords[0] > ff_header['NAXIS1'] - image_edge) |
                                (sources_xy_coords[1] < image_edge) |
                                (sources_xy_coords[0] > ff_header['NAXIS2'] - image_edge))
    if len(bad_sources_mask) > 0:
        catalog.remove_rows(bad_sources_mask)

    return catalog

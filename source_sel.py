import numpy as np
import astropy.coordinates as coord
import astropy.units as u
from astropy.wcs import WCS
from astropy.table import Table
from astroquery.vizier import Vizier

from filesys_io import read_fits_file


# region initial catalog load
def get_hartman_catalog_from_file(catalog_path):
    mag_limit = 16

    catalog = Table.read(catalog_path, format='ascii.fixed_width_no_header',
                         names=('ID', 'RAh', 'RAm', 'RAs', 'DE-', 'DEd', 'DEm', 'DEs',
                                'rmag', 'imag', 'e_rmag', 'e_imag'),
                         col_starts=(0, 7, 10, 13, 19, 20, 23, 26, 38, 45, 58, 64),
                         col_ends=(5, 8, 11, 17, 19, 21, 24, 29, 43, 50, 62, 68))

    print(len(catalog))
    catalog.remove_rows(np.where((catalog['rmag'] >= mag_limit) | (catalog['imag'] >= mag_limit)))
    print(len(catalog))

    catalog.add_column(np.zeros(len(catalog)), index=1, name='RAJ2000')
    catalog['RAJ2000'] = (catalog['RAh'] + catalog['RAm'] / 60 + catalog['RAs'] / 3600)
    catalog.remove_columns(['RAh', 'RAm', 'RAs'])

    catalog.add_column(np.zeros(len(catalog)), index=2, name='DEJ2000')
    catalog['DEJ2000'] = (catalog['DEd'] + catalog['DEm'] / 60 + catalog['DEs'] / 3600)
    catalog['DEJ2000'][np.where(catalog['DE-'] == '-')] *= -1
    catalog.remove_columns(['DE-', 'DEd', 'DEm', 'DEs'])

    # ID              Object identification number
    # RAJ2000  [h]    Right Ascension (J2000)
    # DEJ2000  [d]    Declination (J2000)
    # rmag     [mag]  The r band magnitude
    # imag     [mag]  The i band magnitude
    # e_rmag   [mag]  The 1{sigma} error in rmag
    # e_imag   [mag]  The 1{sigma} error in imag
    print(f"{len(catalog)} sources found in {catalog_path}\n")
    return catalog


def get_hartman_catalog_from_vizier(fits_file_path):
    catalog_id = 'J/ApJ/675/1233/table6'
    catalog_columns = ['ID', 'RAJ2000', 'DEJ2000', 'rmag', 'imag', 'e_rmag', 'e_imag']
    mag_limit = '16'
    catalog_filters = {'rmag': f'<{mag_limit}', 'imag': f'<{mag_limit}'}

    file_header = read_fits_file(fits_file_path)[0]
    file_wcs = WCS(file_header)

    image_center_world = file_wcs.all_pix2world(file_header['NAXIS1'] // 2,
                                                file_header['NAXIS2'] // 2, 0)
    # image_center_world = np.array([file_header['ALPHA'], file_header['DELTA']])
    image_center_skycoord = coord.SkyCoord(ra=image_center_world[0] * u.deg,
                                           dec=image_center_world[1] * u.deg,
                                           frame='icrs')
    # FIXME: 'image_center_skycoord' is JNow while 'cat' is J2000

    search_radius = 0.3
    # Here only CD1_1 value was used because this simple method provides tolerable accuracy
    # search_radius = abs(file_header["CD1_1"]) * ((file_header["NAXIS1"] ** 2 + file_header["NAXIS2"] ** 2) ** 0.5) / 2
    search_angle = coord.Angle(search_radius * u.deg)

    vizier_object = Vizier(columns=catalog_columns,
                           column_filters=catalog_filters)
    vizier_object.ROW_LIMIT = -1
    vizier_result = vizier_object.query_region(image_center_skycoord, radius=search_angle, catalog=catalog_id)[0]

    print(f"{len(vizier_result)} sources found in {catalog_id}\n")

    print(vizier_result)

    # TODO: Doesn't connect to VizieR. Try again later and finish this function:
    # TODO: - [d:m:s] and [h:m:s] to [d] and [h]
# endregion


# region select appropriate sources
def select_sources_by_image_edge(catalog, fits_file_path):
    image_edge = 100  # [pix]

    file_header = read_fits_file(fits_file_path)[0]
    file_wcs = WCS(file_header)

    sources_xy_coords = file_wcs.all_world2pix(catalog['RAJ2000'], catalog['DEJ2000'], 0)
    bad_sources_mask = np.where((sources_xy_coords[0] < image_edge) |
                                (sources_xy_coords[0] > file_header['NAXIS1'] - image_edge) |
                                (sources_xy_coords[1] < image_edge) |
                                (sources_xy_coords[0] > file_header['NAXIS2'] - image_edge))
    if len(bad_sources_mask) > 0:
        catalog.remove_rows(bad_sources_mask)

    return catalog
# endregion

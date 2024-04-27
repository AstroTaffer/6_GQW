import warnings

import numpy as np
import photutils as pht
from astropy.wcs import WCS
from astropy.stats import SigmaClip

from filesys_io import read_fits_file
from res_plot import _draw_sky_map


def rm_sources_aperture_photometry(ff_list, cat, cfg):
    warnings.simplefilter("ignore")

    raw_flux, raw_magn, raw_merr = _aper_phot_core(ff_list, cat, cfg['IMAGE_EDGE'] // 10,
                                                   cfg['APER_RADII'], cfg['OUT_DIR'])

    return raw_flux, raw_magn, raw_merr


def _aper_phot_core(ff_list, cat, img_edge, aper_radii, out_dir):
    apr_num = len(aper_radii)
    img_num = len(ff_list)
    src_num = len(cat)

    raw_flux = np.zeros((apr_num, img_num, src_num))
    raw_magn = np.zeros_like(raw_flux)
    raw_merr = np.zeros_like(raw_flux)  # 1 \sigma

    sky_flux = np.zeros(len(ff_list))
    sky_flux_rms = np.zeros_like(sky_flux)
    if 'SKY' in cat.colnames:
        cat.remove_column('SKY')
    if 'SKY_RMS' in cat.colnames:
        cat.remove_column('SKY_RMS')

    for img_id in range(img_num):
        header, data = read_fits_file(ff_list[img_id])
        wcs = WCS(header)

        src_xy_coords = wcs.all_world2pix(cat['RAJ2000'], cat['DEJ2000'], 0)
        bad_src_mask = ((src_xy_coords[0] < img_edge) |
                        (src_xy_coords[0] > header['NAXIS1'] - img_edge) |
                        (src_xy_coords[1] < img_edge) |
                        (src_xy_coords[1] > header['NAXIS2'] - img_edge))
        src_xy_coords[0][bad_src_mask] = 0
        src_xy_coords[1][bad_src_mask] = 0

        bkg = pht.Background2D(data, (100, 100),
                               filter_size=(9, 9),
                               sigma_clip=SigmaClip(sigma=3.),
                               bkg_estimator=pht.MedianBackground())
        data = data - bkg.background
        # Median sky was already subtracted, so sky_signal points only to residuals
        sky_flux[img_id] = bkg.background_median
        sky_flux_rms[img_id] = bkg.background_rms_median

        ctr_x_coords, ctr_y_coords = pht.centroids.centroid_sources(data,
                                                                    xpos=src_xy_coords[0],
                                                                    ypos=src_xy_coords[1],
                                                                    box_size=9)
        apertures = [pht.CircularAperture(np.vstack((ctr_x_coords,
                                                     ctr_y_coords)).T,
                                          r=r) for r in aper_radii]
        buff_flux = pht.aperture_photometry(data, apertures)

        for apr_id in range(apr_num):
            buff_max = pht.aperture.ApertureStats(data, apertures[apr_id]).max
            raw_flux[apr_id][img_id] = np.where((bad_src_mask | (buff_max > 60000) |
                                                 (buff_flux[f'aperture_sum_{apr_id}'] < 0)),
                                                np.nan, buff_flux[f'aperture_sum_{apr_id}'])

            raw_magn[apr_id][img_id] = -2.5 * np.log10(raw_flux[apr_id][img_id]) + 2.5 * np.log10(header['EXPTIME'])
            raw_merr[apr_id][img_id] = (1.0857 * np.sqrt(raw_flux[apr_id][img_id] * header['GAIN'] +
                                                         apertures[apr_id].area * (sky_flux_rms[img_id] * header['GAIN']
                                                                                   + np.square(header['RDNOISE']))) /
                                        (raw_flux[apr_id][img_id] * header['GAIN']))

        if img_id == 0:
            _draw_sky_map(header, data, wcs, apertures[0], img_edge, out_dir)
            print("Celestial map plotted")

        if img_id % 10 == 9 or img_id == img_num - 1:
            print(f"Aperture photometry: Processed {img_id + 1} images")

    # cat.add_columns([sky_flux, sky_flux_rms], names=['SKY', 'SKY_RMS'])

    return raw_flux, raw_magn, raw_merr

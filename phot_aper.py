import warnings

import photutils as pht
from astropy.wcs import WCS
from astropy.stats import SigmaClip

from filesys_io import *
from res_plot import _draw_sky_map


def rm_sources_aperture_photometry(ff_list, cat, cfg):
    warnings.simplefilter("ignore")

    # raw_flux, raw_magn, raw_merr = _aper_phot_core(ff_list, cat, cfg['IMAGE_EDGE'] // 10,
    #                                                cfg['APER_RADII'], cfg['OUT_DIR'])

    # write_phot_res(raw_flux, raw_magn, raw_merr, cfg['OUT_DIR'], prefix='raw')

    raw_flux, raw_magn, raw_merr = read_phot_res(cfg['OUT_DIR'], prefix='raw')

    return raw_flux, raw_magn, raw_merr


def _aper_phot_core(ff_list, cat, img_edge, aper_radii, out_dir):
    apr_num = len(aper_radii)
    img_num = len(ff_list)
    src_num = len(cat)

    raw_flux = np.zeros((apr_num, img_num, src_num))
    raw_magn = np.zeros_like(raw_flux)
    raw_merr = np.zeros_like(raw_flux)  # 1 \sigma

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
        sky_signal = bkg.background_rms_median

        ctr_x_coords, ctr_y_coords = pht.centroids.centroid_sources(data,
                                                                    xpos=src_xy_coords[0],
                                                                    ypos=src_xy_coords[1],
                                                                    box_size=15)
        apertures = [pht.CircularAperture(np.vstack((ctr_x_coords,
                                                     ctr_y_coords)).T,
                                          r=r) for r in aper_radii]
        buff_flux = pht.aperture_photometry(data, apertures)

        for apr_id in range(apr_num):
            buff_max = pht.aperture.ApertureStats(data, apertures[apr_id]).max
            raw_flux[apr_id][img_id] = np.where(bad_src_mask | (buff_max > 55000),
                                                np.nan, buff_flux[f'aperture_sum_{apr_id}'])

            buff_magn = -2.5 * np.log10(raw_flux[apr_id][img_id]) + 2.5 * np.log10(header['EXPTIME'])
            buff_merr = (1.0857 * np.sqrt(raw_flux[apr_id][img_id] * header['GAIN'] + apertures[apr_id].area *
                                          (sky_signal * header['GAIN'] + header['RDNOISE'] ** 2)) /
                         (raw_flux[apr_id][img_id] * header['GAIN']))
            raw_magn[apr_id][img_id] = np.where(np.isnan(raw_flux[apr_id][img_id]), np.nan, buff_magn)
            raw_merr[apr_id][img_id] = np.where(np.isnan(raw_flux[apr_id][img_id]), np.nan, buff_merr)

        if img_id == 0:
            _draw_sky_map(header, data, wcs, apertures[0], img_edge, out_dir)
            print("Celestial map plotted")

        if img_id % 10 == 9 or img_id == img_num - 1:
            print(f"Aperture photometry: Processed {img_id + 1} images")

    return raw_flux, raw_magn, raw_merr

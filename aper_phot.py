import warnings

import numpy as np
import photutils as pht
from astropy.wcs import WCS
from astropy.stats import SigmaClip

from filesys_io import read_fits_file
from plot_res import draw_sky_map


def rm_sources_aperture_photometry(ff_list, cat, cfg):
    warnings.simplefilter("ignore")
    a, b, c = _phot_core(ff_list, cat, cfg['IMAGE_EDGE'] // 10, cfg['APER_RADII'])
    pass


def _phot_core(ff_list, cat, img_edge, aper_radii):
    apr_num = len(aper_radii)
    img_num = len(ff_list)
    src_num = len(cat)

    raw_flux = np.zeros((apr_num, img_num, src_num))
    raw_magn = np.zeros_like(raw_flux)
    raw_merr = np.zeros_like(raw_flux)  # 1 sigma

    # for _ in range(img_num):
    for _ in range(1):
        header, data = read_fits_file(ff_list[_])
        wcs = WCS(header)

        sources_xy_coords = wcs.all_world2pix(cat['RAJ2000'], cat['DEJ2000'], 0)
        bad_sources_mask = np.where((sources_xy_coords[0] < img_edge) |
                                    (sources_xy_coords[0] > header['NAXIS1'] - img_edge) |
                                    (sources_xy_coords[1] < img_edge) |
                                    (sources_xy_coords[1] > header['NAXIS2'] - img_edge))
        if len(bad_sources_mask) > 0:
            sources_xy_coords[0][bad_sources_mask] = 0
            sources_xy_coords[1][bad_sources_mask] = 0

        bkg = pht.Background2D(data, (100, 100),
                               filter_size=(9, 9),
                               sigma_clip=SigmaClip(sigma=3.),
                               bkg_estimator=pht.MedianBackground())
        data = data - bkg.background
        sky_signal = bkg.background_rms_median

        centroids_x_coords, centroids_y_coords = pht.centroids.centroid_sources(data,
                                                                                xpos=sources_xy_coords[0],
                                                                                ypos=sources_xy_coords[1],
                                                                                box_size=15)
        apertures = [pht.CircularAperture(np.vstack((centroids_x_coords,
                                                     centroids_y_coords)).T,
                                          r=r) for r in aper_radii]

        buff_flux = pht.aperture_photometry(data, apertures)

        # TODO: Delete overexposed sources

        for __ in range(len(aper_radii)):
            raw_flux[__][_] = buff_flux[f'aperture_sum_{__}']
            raw_magn[__][_] = -2.5 * np.log10(raw_flux[__][_]) + 2.5 * np.log10(header['EXPTIME'])
            raw_merr[__][_] = (1.0857 * np.sqrt(raw_flux[__][_] * header['GAIN'] +
                                                apertures[__].area * (sky_signal * header['GAIN'] +
                                                                      header['RDNOISE'] ** 2)) /
                               raw_flux[__][_] * header['GAIN'])

            raw_flux[__][_][bad_sources_mask] = np.nan
            raw_magn[__][_][bad_sources_mask] = np.nan
            raw_merr[__][_][bad_sources_mask] = np.nan

        if _ == 0:
            draw_sky_map(header, data, wcs, apertures[0], img_edge)
            print("Celestial map plotted")

        if _ % 10 == 0 or _ == img_num - 1:
            print(f"Aperture photometry: {_} images ready")

    return raw_flux, raw_magn, raw_merr

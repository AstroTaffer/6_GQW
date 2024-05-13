import warnings

import numpy as np
from astropy.wcs import WCS
from astropy.stats import SigmaClip, gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel, convolve
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import detect_sources, SourceCatalog
from photutils.centroids import centroid_sources
from photutils.aperture import CircularAperture, aperture_photometry, ApertureStats
from photutils.utils import circular_footprint

from filesys_io import read_fits_file
from res_plot import _draw_sky_map


def rm_sources_aperture_photometry(ff_list, cat, cfg, prefix='unknw'):
    warnings.simplefilter("ignore")
    print('')

    raw_flux, raw_magn, raw_merr = _aper_phot_core(ff_list, cat, cfg['IMAGE_EDGE'] // 10,
                                                   cfg['APER_RADII'], cfg['BEST_APER_ID'], f"{cfg['OUT_DIR']}{prefix}\\{prefix}_")

    return raw_flux, raw_magn, raw_merr


def _aper_phot_core(ff_list, cat, img_edge, aper_radii, ba_id, out_dir):
    apr_num = len(aper_radii)
    img_num = len(ff_list)
    src_num = len(cat)

    raw_flux = np.zeros((apr_num, img_num, src_num))
    raw_magn = np.zeros_like(raw_flux)
    raw_merr = np.zeros_like(raw_flux)  # 1 \sigma

    fwhm = np.zeros(img_num)

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

        bkg = Background2D(data, box_size=32, filter_size=9, sigma_clip=SigmaClip(sigma=3.),
                           bkg_estimator=MedianBackground())
        data = data - bkg.background
        # Median sky was already subtracted, so sky_signal points only to residuals
        sky_flux_rms = bkg.background_rms_median

        kernel = Gaussian2DKernel(3.0 * gaussian_fwhm_to_sigma, x_size=3, y_size=3)
        kernel.normalize()
        # I cut a corner here by using sky_flux_rms
        segm_data = detect_sources(convolve(data, kernel), 50 * sky_flux_rms, npixels=5)
        fwhm[img_id] = np.nanmedian(SourceCatalog(data, segm_data).fwhm.value)

        circ_ftpt = circular_footprint(int(fwhm[img_id] * 2))

        ctr_x_coords, ctr_y_coords = centroid_sources(data,
                                                      xpos=src_xy_coords[0],
                                                      ypos=src_xy_coords[1],
                                                      footprint=circ_ftpt)
        apertures = [CircularAperture(np.vstack((ctr_x_coords, ctr_y_coords)).T, r=r) for r in aper_radii]
        buff_flux = aperture_photometry(data, apertures)

        for apr_id in range(apr_num):
            buff_max = ApertureStats(data, apertures[apr_id]).max
            # raw_flux[apr_id][img_id] = buff_flux[f'aperture_sum_{apr_id}']
            raw_flux[apr_id][img_id] = np.where((bad_src_mask | (buff_max > 60000) |
                                                 (buff_flux[f'aperture_sum_{apr_id}'] < 0)),
                                                np.nan, buff_flux[f'aperture_sum_{apr_id}'])

            raw_magn[apr_id][img_id] = -2.5 * np.log10(raw_flux[apr_id][img_id]) + 2.5 * np.log10(header['EXPTIME'])
            raw_merr[apr_id][img_id] = (1.0857 * np.sqrt(raw_flux[apr_id][img_id] * header['GAIN'] +
                                                         apertures[apr_id].area * (sky_flux_rms * header['GAIN']
                                                                                   + np.square(header['RDNOISE']))) /
                                        (raw_flux[apr_id][img_id] * header['GAIN']))

        if img_id == 0:
            _draw_sky_map(header, data, wcs, apertures[ba_id], img_edge, out_dir)
            print("Celestial map plotted")

        if img_id % 10 == 9 or img_id == img_num - 1:
            print(f"Aperture photometry: Processed {img_id + 1} images")

    print(f'Median FWHM is {np.median(fwhm)} ({src_num})')

    return raw_flux, raw_magn, raw_merr


def _aper_phot_altfull(mycat):
    import os

    import astropy.io.fits as fits
    from astropy.stats import SigmaClip, mad_std, sigma_clip

    # disable warnings
    import warnings
    warnings.simplefilter("ignore")

    path_to_data = r'H:\FITS IMAGES\2024-03-19 M37\CALIBRATED\i'

    aper_radii = np.array([3, 5, 7, 9])

    file_list = []
    dir_content = os.listdir(path_to_data)
    for f in dir_content:
        if f.count('.fits'):
            file_list.append(path_to_data + '\\' + f)
    counter = len(file_list)

    df = open(path_to_data + '\\Result.txt', 'w')
    df.write('DATE-OBS\tEXPTIME\tSky\tS_Sky\tFWHM\n')

    outputs = [open(path_to_data + '/Phot' + str(r) + '.txt', 'a') for r in aper_radii]

    objects = mycat

    for file in file_list:
        print('----------<>-----------')
        print('Frame: ' + str(counter), file)
        counter = counter - 1
        hdulist = fits.open(file)
        header = hdulist[0].header
        data = hdulist[0].data
        hdulist.verify('fix')
        hdulist.close()
        w = WCS(header)

        df.write(header['DATE-OBS'] + '\t')
        df.write(str(header['EXPTIME']) + '\t')

        sigmaclip = SigmaClip(sigma=3.)
        bkg_estimator = MedianBackground()
        bkg = Background2D(data, box_size=32, filter_size=9, sigma_clip=sigmaclip, bkg_estimator=bkg_estimator)
        data = data - bkg.background
        sky = np.median(bkg.background)
        s_sky = sigma_clip(data, stdfunc=mad_std).filled(np.nan)
        s_sky = np.nanstd(s_sky)

        print('Sky = ', sky)
        print('S_Sky = ', s_sky)
        df.write('{:.1f}'.format(sky) + '\t')
        df.write('{:.2f}'.format(s_sky) + '\t')

        x, y = w.wcs_world2pix(objects['RAJ2000'], objects['DEJ2000'], 0)

        sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
        kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
        kernel.normalize()
        segm = detect_sources(convolve(data, kernel), 50 * s_sky, npixels=5)
        cat = SourceCatalog(data, segm)
        fwhm = round(np.nanmedian(cat.fwhm.value), 2)
        print('FWHM = ', fwhm)
        df.write(str(fwhm) + '\n')

        footprint = circular_footprint(int(fwhm * 2))
        xc, yc = centroid_sources(data, x, y, footprint=footprint)
        positions = np.transpose([xc, yc])
        stellar_aper = [CircularAperture(positions, r=r) for r in aper_radii]

        stellar_phot = aperture_photometry(data, stellar_aper, method='exact')

        if counter == len(file_list) - 1:
            _draw_sky_map(header, data, w, stellar_aper[-1], 10, path_to_data)

        for count, value in enumerate(aper_radii):
            flux = stellar_phot['aperture_sum_' + str(count)].value
            # flux_err = np.sqrt(flux + 3.14 * value * value * s_sky * s_sky)
            np.savetxt(outputs[count], flux, fmt='%1.3f', delimiter='\t', newline='\t')
            outputs[count].write('\n')

    for file in outputs:
        file.close()

    df.close()

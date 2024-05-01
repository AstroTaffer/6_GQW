import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u


def rm_sources_ensemble_photometry(raw_magn, raw_merr, cat, cfg):
    # clr_magn, clr_merr = _diff_phot_core(raw_magn, raw_merr, cat, cfg['APER_RADII'],
    #                                      cfg['ENS_MAX_MAG_DIFF'], cfg['ENS_MAX_SIGMA_CRIT'])
    clr_magn, clr_merr = _diff_phot_altcore(raw_magn, raw_merr, cat, cfg['APER_RADII'],
                                            cfg['ENS_INIT_SEARCH_R'], cfg['ENS_MAX_SEARCH_R'],
                                            cfg['ENS_MAX_MAG_DIFF'], cfg['ENS_MAX_SIGMA_CRIT'])
    return clr_magn, clr_merr


def _diff_phot_core(raw_magn, raw_merr, cat, aper_radii, ens_mmd, ens_msc):
    clr_magn = np.zeros_like(raw_magn)
    clr_merr = np.zeros_like(raw_magn)
    ens_corr = 0
    ens_src_merr = 0

    for apr_id in range(len(aper_radii)):
        success_count = 0
        # Technically, if raw_magn is NaN then raw_merr is NaN and vice versa
        # But I'm so stressed right now, that I'll double-check that
        src_is_all_fin = np.isfinite(raw_magn[apr_id]).all(axis=0) & np.isfinite(raw_merr[apr_id]).all(axis=0)
        src_is_all_nan = np.isnan(raw_magn[apr_id]).all(axis=0) | np.isnan(raw_merr[apr_id]).all(axis=0)

        for target_id in range(len(cat)):
            # If target star is "all-NaN"
            if src_is_all_nan[target_id]:
                clr_magn[apr_id, :, target_id] = np.nan
                clr_merr[apr_id, :, target_id] = np.nan
                continue

            src_imag_diff = np.abs(cat['imag'] - cat['imag'][target_id])
            # src_rmag_diff = np.abs(cat['rmag'] - cat['rmag'][target_id])

            # Pick the sources that are:
            #   2. Not too bright or too dim
            #   3. "All-finite"
            # NOTE: Target_star is included if not excluded by 3rd criteria
            ens_src_ids = np.where((src_imag_diff <= ens_mmd) &
                                   # (src_rmag_diff <= ens_mmd) &
                                   src_is_all_fin)[0]

            # Step 2: Set the sigma_crit loop
            while len(ens_src_ids) >= 10:
                ens_src_magn = raw_magn[apr_id, :, ens_src_ids].T
                ens_src_merr = raw_merr[apr_id, :, ens_src_ids].T

                # Step 2.1: Calculate mean weighted magnitude for each image
                ens_src_weight = np.ones(len(ens_src_ids)) / np.nanmean(np.square(ens_src_merr), axis=0)
                ens_image_mw_magn = np.nansum(ens_src_magn * ens_src_weight, axis=1) / np.sum(ens_src_weight)

                # Step 2.2: Apply correction to ensemble sources
                ens_corr = ens_image_mw_magn - np.nanmean(ens_image_mw_magn)
                ens_src_magn = ens_src_magn - ens_corr.reshape((-1, 1))

                # Step 2.3: Check ensemble sources magnitude sigma-criteria
                ens_src_sigma_crit = np.nanstd(ens_src_magn, axis=0) / np.nanmean(ens_src_merr, axis=0)
                if np.max(ens_src_sigma_crit) <= ens_msc:
                    break
                ens_src_ids = np.delete(ens_src_ids, np.argmax(ens_src_sigma_crit))

            # If couldn't find enough good sources
            if len(ens_src_ids) < 10:
                clr_magn[apr_id, :, target_id] = np.nan
                clr_merr[apr_id, :, target_id] = np.nan
                continue

            # Step 3: Calculate clear magnitudes and errors
            success_count += 1
            clr_magn[apr_id, :, target_id] = raw_magn[apr_id, :, target_id] - ens_corr
            clr_merr[apr_id, :, target_id] = np.sqrt(np.square(raw_merr[apr_id, :, target_id]) +
                                                     1 / np.nansum(1 / np.square(ens_src_merr), axis=1))

            if target_id % 50 == 49 or target_id == len(cat) - 1:
                print(f"Ensemble photometry (aperture #{apr_id + 1}): Processed {target_id + 1} stars, "
                      f"{success_count} successfully")

    return clr_magn, clr_merr


def _diff_phot_altcore(raw_magn, raw_merr, cat, aper_radii, ens_isr, ens_msr, ens_mmd, ens_msc):
    clr_magn = np.zeros_like(raw_magn)
    clr_merr = np.zeros_like(raw_magn)

    src_sc = SkyCoord(ra=cat['RAJ2000'] * u.deg, dec=cat['DEJ2000'] * u.deg, frame='icrs')

    for apr_id in range(len(aper_radii)):
        success_count = 0
        # Technically, if raw_magn is NaN then raw_merr is NaN and vice versa
        # But I'm so stressed right now, that I'll double-check that
        src_is_all_fin = np.isfinite(raw_magn[apr_id]).all(axis=0) & np.isfinite(raw_merr[apr_id]).all(axis=0)
        src_is_all_nan = np.isnan(raw_magn[apr_id]).all(axis=0) | np.isnan(raw_merr[apr_id]).all(axis=0)

        for target_id in range(len(cat)):
            # If target star is "all-NaN"
            if src_is_all_nan[target_id]:
                clr_magn[apr_id, :, target_id] = np.nan
                clr_merr[apr_id, :, target_id] = np.nan
                continue

            # Step 1: Set search_radius loop
            search_radius = ens_isr
            src_sep = src_sc[target_id].separation(src_sc)
            src_imag_diff = np.abs(cat['imag'] - cat['imag'][target_id])
            # src_rmag_diff = np.abs(cat['rmag'] - cat['rmag'][target_id])

            while search_radius <= ens_msr:
                # Pick the sources that are:
                #   1. Not too far
                #   2. Not too bright or too dim
                #   3. "All-finite"
                # NOTE: Target_star is included if not excluded by 3rd criteria
                ens_src_ids = np.where((src_sep <= search_radius * u.arcmin) &
                                       (src_imag_diff <= ens_mmd) &
                                       # (src_rmag_diff <= ens_mmd) &
                                       src_is_all_fin)[0]

                # If couldn't find enough sources at all
                if len(ens_src_ids) < 10:
                    search_radius += 1
                    continue

                # Step 2: Set the sigma_crit loop
                while len(ens_src_ids) >= 10:
                    ens_src_magn = raw_magn[apr_id, :, ens_src_ids].T
                    ens_src_merr = raw_merr[apr_id, :, ens_src_ids].T

                    # Step 2.1: Calculate mean weighted magnitude for each image
                    ens_src_weight = np.ones(len(ens_src_ids)) / np.nanmean(np.square(ens_src_merr), axis=0)
                    ens_image_mw_magn = np.nansum(ens_src_magn * ens_src_weight, axis=1) / np.sum(ens_src_weight)

                    # Step 2.2: Apply correction to ensemble sources
                    ens_corr = ens_image_mw_magn - np.nanmean(ens_image_mw_magn)
                    ens_src_magn = ens_src_magn - ens_corr.reshape((-1, 1))

                    # Step 2.3: Check ensemble sources magnitude sigma-criteria
                    ens_src_sigma_crit = np.nanstd(ens_src_magn, axis=0) / np.nanmean(ens_src_merr, axis=0)
                    if np.max(ens_src_sigma_crit) <= ens_msc:
                        break
                    ens_src_ids = np.delete(ens_src_ids, np.argmax(ens_src_sigma_crit))

                # If couldn't find enough good sources
                if len(ens_src_ids) < 10:
                    search_radius += 1
                    continue

                # Step 3: Calculate clear magnitudes and errors
                success_count += 1
                clr_magn[apr_id, :, target_id] = raw_magn[apr_id, :, target_id] - ens_corr
                clr_merr[apr_id, :, target_id] = np.sqrt(np.square(raw_merr[apr_id, :, target_id]) +
                                                         1 / np.nansum(1 / np.square(ens_src_merr), axis=1))
                break

            if search_radius > ens_msr:
                clr_magn[apr_id, :, target_id] = np.nan
                clr_merr[apr_id, :, target_id] = np.nan

            if target_id % 50 == 49 or target_id == len(cat) - 1:
                print(f"Ensemble photometry (aperture #{apr_id + 1}): Processed {target_id + 1} stars, "
                      f"{success_count} successfully")

    return clr_magn, clr_merr


def _diff_phot_fluxcore(raw_flux):
    clr_flux = np.zeros_like(raw_flux)
    clr_magn = np.zeros_like(clr_flux)

    for apr_id in range(raw_flux.shape[0]):
        ens_src_ids = np.where(np.isfinite(raw_flux[apr_id]).all(axis=0))[0]

        while len(ens_src_ids) >= 10:
            ens_flux = raw_flux[apr_id, :, ens_src_ids].T

            ens_trend = np.sum(ens_flux, axis=1)
            ens_trend /= np.mean(ens_trend)
            ens_flux /= ens_trend.reshape((-1, 1))

            ens_std = np.std(ens_flux, axis=0) / np.sqrt(np.mean(ens_flux))
            if np.max(ens_std) <= 5:
                break
            ens_src_ids = np.delete(ens_src_ids, np.argmax(ens_std))

        if len(ens_src_ids) < 10:
            clr_flux[apr_id] = np.nan
            clr_magn[apr_id] = np.nan
            continue

        # noinspection PyUnboundLocalVariable
        clr_flux[apr_id] = raw_flux[apr_id] / ens_trend.reshape((-1, 1))
        clr_magn[apr_id] = -2.5 * np.log10(clr_flux[apr_id]) + 2.5 * np.log10(80)

    return clr_flux, clr_magn

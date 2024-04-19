from astropy.coordinates import SkyCoord
import astropy.units as u

from filesys_io import *


def rm_sources_ensemble_photometry(raw_magn, raw_merr, cat, cfg):
    clr_magn, clr_merr = _diff_phot_classic_core(raw_magn, raw_merr, cat, cfg['APER_RADII'],
                                                 cfg['ENS_INIT_SEARCH_R'], cfg['ENS_MAX_SEARCH_R'],
                                                 cfg['ENS_MAX_MAG_DIFF'], cfg['ENS_MAX_SIGMA_CRIT'])
    # write_phot_res(None, clr_magn, clr_merr, cfg['OUT_DIR'], prefix='clr')

    # clr_flux, clr_magn, clr_merr = read_phot_res(cfg['OUT_DIR'], prefix='clr')
    return clr_magn, clr_merr


# FIXME: This function is extremely slow (0.1 star/sec)
def _diff_phot_classic_core(raw_magn, raw_merr, cat, aper_radii, ens_isr, ens_msr, ens_mmd, ens_msc):
    clr_magn = np.zeros_like(raw_magn)
    clr_merr = np.zeros_like(raw_magn)

    # Quick plug for PyCharm to stop worrying
    ens_src_merr = 0
    ens_corr = 0

    sources_sc = SkyCoord(ra=cat['RAJ2000'] * u.deg, dec=cat['DEJ2000'] * u.deg, frame='icrs')

    for aper_id in range(len(aper_radii)):
        success_counter = 0
        for target_id in range(len(cat)):
            # If target star is "all-NaN"
            if np.isnan(raw_magn[aper_id, :, target_id]).all() or np.isnan(raw_merr[aper_id, :, target_id]).all():
                clr_magn[aper_id, :, target_id] = np.nan
                clr_merr[aper_id, :, target_id] = np.nan
                continue

            # Step 1: Set search_radius loop
            search_radius = ens_isr
            while search_radius <= ens_msr:
                print(search_radius)
                ens_src_ids = []
                for _ in range(len(cat)):
                    # Skip the source if any of this is true:
                    #   1. It's too far
                    #   2. It's too different in terms of rmag or imag
                    #   3. It's "all-NaN"
                    #   4. It's the target star
                    if (sources_sc[target_id].separation(sources_sc[_]) > search_radius * u.arcmin or
                            abs(cat[target_id]['rmag'] - cat[_]['rmag']) > ens_mmd or
                            abs(cat[target_id]['imag'] - cat[_]['imag']) > ens_mmd or
                            np.isnan(raw_magn[aper_id, :, _]).all() or np.isnan(raw_merr[aper_id, :, _]).all() or
                            _ == target_id):
                        continue
                    ens_src_ids.append(_)

                # If couldn't find enough sources
                if len(ens_src_ids) < 10:
                    search_radius += 1
                    continue

                # Step 2: Set the rejection loop
                while len(ens_src_ids) >= 10:
                    ens_src_magn = np.zeros((np.shape(raw_magn)[1], len(ens_src_ids)))
                    ens_src_merr = np.zeros_like(ens_src_magn)

                    for _ in range(len(ens_src_ids)):
                        ens_src_magn[:, _] = raw_magn[aper_id, :, ens_src_ids[_]]
                        ens_src_merr[:, _] = raw_merr[aper_id, :, ens_src_ids[_]]

                    # Step 2.1: Calculate mean weighted magnitude for each image
                    ens_src_weight = np.ones(len(ens_src_ids)) / np.nanmean(np.square(ens_src_merr), axis=0)
                    ens_image_mw_magn = np.nansum(ens_src_magn * ens_src_weight, axis=1) / np.nansum(ens_src_weight)

                    # Step 2.2: Apply correction to ensemble sources
                    ens_corr = ens_image_mw_magn - np.nanmean(ens_image_mw_magn)
                    ens_src_magn -= ens_corr.reshape((-1, 1))

                    # Step 2.3: Check ensemble sources magnitude sigma-criteria
                    ens_src_sigma_crit = np.nanstd(ens_src_magn, axis=0) / np.nanmean(ens_src_merr, axis=0)
                    if np.nanmax(ens_src_sigma_crit) <= ens_msc:
                        break
                    ens_src_ids.pop(np.nanargmax(ens_src_sigma_crit))

                # If couldn't find enough good sources
                if len(ens_src_ids) < 10:
                    search_radius += 1
                    continue

                # Step 3: Calculate clear magnitudes and errors
                success_counter += 1
                clr_magn[aper_id, :, target_id] = raw_magn[aper_id, :, target_id] - ens_corr
                clr_merr[aper_id, :, target_id] = np.sqrt(np.square(raw_merr[aper_id, :, target_id]) +
                                                          1 / np.nansum(1 / np.square(ens_src_merr), axis=1))
                break

            if search_radius > ens_msr:
                clr_magn[aper_id, :, target_id] = np.nan
                clr_merr[aper_id, :, target_id] = np.nan

            # print(f"Ensemble photometry (aperture number {aper_id + 1}): Processed {target_id + 1} stars, "
            #       f"{success_counter} successfully")
            if target_id % 100 == 99 or target_id == len(cat) - 1:
                print(f"Ensemble photometry (aperture number {aper_id + 1}): Processed {target_id + 1} stars, "
                      f"{success_counter} successfully")

    return clr_magn, clr_merr

from filesys_io import *


def rm_get_results(clr_magn, clr_merr, cat, cfg):
    ba_ids = _select_best_aperture_by_magn_std(clr_magn)
    cat = _create_best_aperture_column(ba_ids, cat, cfg['APER_RADII'])
    res_magn, res_merr = _get_results_by_ba_id(ba_ids, clr_magn, clr_merr)
    # write_sel_src_cat(catalog, cfg['OUT_DIR'])

    # catalog = read_sel_src_cat(cfg['OUT_DIR'])
    return cat, res_magn, res_merr


def _select_best_aperture_by_magn_std(clr_magn):
    best_aper_ids = [0] * clr_magn.shape[2]
    for star_id in range(len(best_aper_ids)):
        try:
            best_aper_ids[star_id] = np.nanargmin(np.nanstd(clr_magn[:, :, star_id], axis=1))
        except ValueError:
            best_aper_ids[star_id] = -1
    return best_aper_ids


def _get_results_by_ba_id(best_aper_ids, clr_magn, clr_merr):
    res_magn = np.zeros(clr_magn.shape[2])
    res_merr = np.zeros_like(res_magn)

    for star_id in range(len(res_magn)):
        if best_aper_ids[star_id] == -1:
            res_magn[star_id] = np.nan
            res_merr[star_id] = np.nan
            continue

        res_magn[star_id] = np.nanmean(clr_magn[best_aper_ids[star_id], :, star_id])
        res_merr[star_id] = np.nanmean(clr_merr[best_aper_ids[star_id], :, star_id])

    return res_magn, res_merr


def _create_best_aperture_column(best_aper_ids, cat, aper_radii):
    best_apers = [-1] * len(cat)
    for star_id in range(len(cat)):
        if best_aper_ids[star_id] == -1:
            continue
        best_apers[star_id] = aper_radii[best_aper_ids[star_id]]
    cat.add_column(best_apers, name='BEAR')
    return cat

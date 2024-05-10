import numpy as np
from astropy.stats import sigma_clip


def rm_get_final_results(rb_flux, rb_magn, rb_merr, cat, aper_radii, best_aper, flt_cname):
    rb_flux, rb_magn, rb_merr = _select_best_aperture(rb_flux, rb_magn, rb_merr,
                                                      aper_radii, best_aper)

    med_rb_magn, cat, good_mask = _select_magn_by_sigma_clip(rb_magn, cat, flt_cname)
    med_rb_flux, med_rb_merr = _select_flux_and_merr_by_mask(rb_flux, rb_merr, good_mask)

    return med_rb_flux, med_rb_magn, med_rb_merr, cat


def _select_best_aperture(rb_flux, rb_magn, rb_merr, aper_radii, best_aper):
    ba_id = aper_radii.index(best_aper)
    if rb_flux is not None:
        rb_flux = rb_flux[ba_id]
    rb_magn = rb_magn[ba_id]
    if rb_merr is not None:
        rb_merr = rb_merr[ba_id]
    return rb_flux, rb_magn, rb_merr


def _select_magn_by_sigma_clip(rb_magn, cat, flt_cname):
    med_rb_magn = np.nanmedian(rb_magn, axis=0)
    cat_magn = cat[flt_cname]

    delta_magn = med_rb_magn - cat_magn
    good_mask = ~np.ma.getmaskarray(sigma_clip(delta_magn, sigma=3, maxiters=5, masked=True))
    print(f'Sigma clipping left {np.sum(good_mask)} stars in {flt_cname}')

    med_rb_magn = med_rb_magn[good_mask]
    cat = cat[good_mask]

    return med_rb_magn, cat, good_mask


def _select_flux_and_merr_by_mask(rb_flux, rb_merr, good_mask):
    if rb_flux is not None:
        med_rb_flux = np.nanmedian(rb_flux, axis=0)
        med_rb_flux = med_rb_flux[good_mask]
    else:
        med_rb_flux = None
    if rb_merr is not None:
        med_rb_merr = np.nanmedian(rb_merr, axis=0)
        med_rb_merr = med_rb_merr[good_mask]
    else:
        med_rb_merr = None
    return med_rb_flux, med_rb_merr

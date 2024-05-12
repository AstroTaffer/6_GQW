import numpy as np
from astropy.stats import sigma_clip
from astropy.modeling.fitting import LinearLSQFitter, FittingWithOutlierRemoval
from astropy.modeling.models import Linear1D

from res_plot import _plot_raw_magn, _plot_fitted_magn, _plot_std_magn


def rm_get_results(flux, magn, merr, cat, cfg, flt_cname):
    flux, magn, merr, cat = _select_any_finite_ba_slice(flux, magn, merr, cat, cfg['BEST_APER_ID'])

    inliers_mask = _get_sigma_clip_mask(magn, cat, flt_cname, cfg['OUT_DIR'])
    flux, magn, merr, cat = _select_inliers(flux, magn, merr, cat, inliers_mask)

    _calc_magn_std(magn, cat, flt_cname, cfg['OUT_DIR'])


def _select_any_finite_ba_slice(flux, magn, merr, cat, ba_id):
    finite_mask = np.ones(len(cat), dtype=bool)

    if flux is not None:
        flux = flux[ba_id]
        finite_mask = finite_mask & np.isfinite(flux).any(axis=0)
    if magn is not None:
        magn = magn[ba_id]
        finite_mask = finite_mask & np.isfinite(magn).any(axis=0)
    if merr is not None:
        merr = merr[ba_id]
        finite_mask = finite_mask & np.isfinite(merr).any(axis=0)

    if flux is not None:
        flux = flux[:, finite_mask]
    if magn is not None:
        magn = magn[:, finite_mask]
    if merr is not None:
        merr = merr[:, finite_mask]
    cat = cat[finite_mask]

    print(f"{np.sum(finite_mask)} finite stars in aperture slice [{ba_id}]")
    return flux, magn, merr, cat


def _select_inliers(flux, magn, merr, cat, inliers_mask):
    if flux is not None:
        flux = flux[:, inliers_mask]
    if magn is not None:
        magn = magn[:, inliers_mask]
    if merr is not None:
        merr = merr[:, inliers_mask]
    cat = cat[inliers_mask]

    return flux, magn, merr, cat


def _get_sigma_clip_mask(rb_magn, cat, flt_cname, out_dir):
    mrb_magn = np.nanmedian(rb_magn, axis=0)
    cat_magn = cat[flt_cname]

    _plot_raw_magn(mrb_magn, cat_magn, flt_cname[0], out_dir)

    # noinspection PyTypeChecker
    sc_fit = FittingWithOutlierRemoval(LinearLSQFitter(), sigma_clip, niter=5, sigma=3.0)
    fitted_line, inliers_mask = sc_fit(Linear1D(), mrb_magn, cat_magn)
    inliers_mask = np.invert(inliers_mask)

    # TODO: Round fitted parameters and plot 1\sigma uncertainty
    _plot_fitted_magn(mrb_magn, cat_magn, flt_cname[0], out_dir, inliers_mask, fitted_line)
    print(f"{np.sum(inliers_mask)} stars, {sc_fit.fit_info['niter']} iterations, parameters {fitted_line.parameters}")

    return inliers_mask


def _calc_magn_std(rb_magn, cat, flt_cname, out_dir):
    std_rb_magn = np.nanstd(rb_magn, axis=0)
    cat_magn = cat[flt_cname]

    otw_data = np.loadtxt(f"{out_dir}1.2M_errors_{flt_cname[0]}.txt")
    otw_magn = otw_data[:, 0]
    otw_std = otw_data[:, 1]

    _plot_std_magn(cat_magn, std_rb_magn, otw_magn, otw_std, flt_cname[0], out_dir)

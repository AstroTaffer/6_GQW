import numpy as np
from astropy.stats import sigma_clip
from astropy.modeling.fitting import LinearLSQFitter, FittingWithOutlierRemoval
from astropy.modeling.models import Linear1D

from res_plot import (_plot_raw_magn, _plot_fitted_magn, _plot_std_magn, _plot_merr,
                      _plot_color_term, _plot_alt_color_term)


def rm_get_results(flux, magn, merr, cat, cfg, flt_cname):
    print('')
    out_dir = cfg['OUT_DIR']

    flux, magn, merr, cat = _select_any_finite_ba_slice(flux, magn, merr, cat, cfg['BEST_APER_ID'])

    fitted_line,inliers_mask = _get_sigma_clip_mask(magn, cat, flt_cname, out_dir)

    flux, magn, merr, cat = _select_inliers(flux, magn, merr, cat, inliers_mask)
    clc_magn = fitted_line(magn)

    _calc_magn_std(clc_magn, flt_cname, out_dir)

    if merr is not None:
        _calc_merr(clc_magn, merr, flt_cname, out_dir)

    if flux is None:
        flux = 80 * np.power(10.0, -0.4 * magn)
    _calc_total_throughput(clc_magn, flux, flt_cname, out_dir)

    return clc_magn, cat


def rm_get_color_terms(r_clc_magn, i_clc_magn, r_cat, i_cat, cfg):
    print('')
    out_dir = cfg['OUT_DIR']

    inter, r_ids, i_ids = np.intersect1d(r_cat['ID'], i_cat['ID'], return_indices=True)

    r_clc_magn_med = np.nanmedian(r_clc_magn[:, r_ids], axis=0)
    i_clc_magn_med = np.nanmedian(i_clc_magn[:, i_ids], axis=0)
    cat = r_cat[r_ids]

    # TODO: Color = delta_RB or delta_CAT?
    cat_color = cat['rmag'] - cat['imag']
    r_mag_delta = cat['rmag'] - r_clc_magn_med
    i_mag_delta = cat['imag'] - i_clc_magn_med

    lin_fit = LinearLSQFitter()
    r_fitted_line = lin_fit(Linear1D(), cat_color, r_mag_delta)
    i_fitted_line = lin_fit(Linear1D(), cat_color, i_mag_delta)
    print(f"rmag color parameters {r_fitted_line.parameters}")
    print(f"imag color parameters {i_fitted_line.parameters}")

    _plot_color_term(cat_color, r_mag_delta, 'r', out_dir, r_fitted_line)
    _plot_color_term(cat_color, i_mag_delta, 'i', out_dir, i_fitted_line)

    clc_color = r_clc_magn_med - i_clc_magn_med

    alt_fitted_line = lin_fit(Linear1D(), cat_color, clc_color)
    print(f"alt color parameters {alt_fitted_line.parameters}")

    _plot_alt_color_term(cat_color, clc_color, out_dir, alt_fitted_line)


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


def _get_sigma_clip_mask(rb_magn, cat, flt_cname, out_dir):
    rb_magn_med = np.nanmedian(rb_magn, axis=0)
    cat_magn = cat[flt_cname]

    _plot_raw_magn(rb_magn_med, cat_magn, flt_cname[0], out_dir)

    # noinspection PyTypeChecker
    sc_fit = FittingWithOutlierRemoval(LinearLSQFitter(), sigma_clip, niter=5, sigma=3.0)
    fitted_line, inliers_mask = sc_fit(Linear1D(), rb_magn_med, cat_magn)
    inliers_mask = np.invert(inliers_mask)

    # TODO: Round fitted parameters and plot 1\sigma uncertainty
    _plot_fitted_magn(rb_magn_med, cat_magn, flt_cname[0], out_dir, inliers_mask, fitted_line)
    print(f"{np.sum(inliers_mask)} stars, {sc_fit.fit_info['niter']} iterations, parameters {fitted_line.parameters}")

    return fitted_line, inliers_mask


def _select_inliers(flux, magn, merr, cat, inliers_mask):
    if flux is not None:
        flux = flux[:, inliers_mask]
    if magn is not None:
        magn = magn[:, inliers_mask]
    if merr is not None:
        merr = merr[:, inliers_mask]
    cat = cat[inliers_mask]

    return flux, magn, merr, cat


def _calc_magn_std(clc_magn, flt_cname, out_dir):
    clc_magn_med = np.nanmedian(clc_magn, axis=0)
    clc_magn_std = np.nanstd(clc_magn, axis=0)

    otw_data = np.loadtxt(f"{out_dir}1.2M_errors_{flt_cname[0]}.txt")
    otw_magn = otw_data[:, 0]
    otw_std = otw_data[:, 1]

    _plot_std_magn(clc_magn_med, clc_magn_std, otw_magn, otw_std, flt_cname[0], out_dir)


def _calc_merr(clc_magn, rb_merr, flt_cname, out_dir):
    clc_magn_med = np.nanmedian(clc_magn, axis=0)
    rb_merr_med = np.nanmedian(rb_merr, axis=0)

    _plot_merr(clc_magn_med, rb_merr_med, flt_cname[0], out_dir)


def _calc_total_throughput(clc_magn, rb_flux, flt_cname, out_dir):
    clc_magn_med = np.nanmedian(clc_magn, axis=0)
    rb_flux_med = np.nanmedian(rb_flux, axis=0)

    # TODO: Calc tp [%]

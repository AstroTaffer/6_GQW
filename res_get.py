import numpy as np
from astropy.stats import sigma_clip
from astropy.modeling.fitting import LinearLSQFitter, FittingWithOutlierRemoval
from astropy.modeling.models import Linear1D
import matplotlib.pyplot as plt

from res_plot import (_plot_raw_magn, _plot_fitted_magn, _plot_std_magn, _plot_merr, _plot_delta_magn, _plot_check_magn,
                      _plot_total_throughput, _plot_color_term, _plot_alt_color_term)


def rm_get_results(flux, magn, merr, cat, cfg, flt_cname):
    print('')
    out_dir = cfg['OUT_DIR']

    flux, magn, merr, cat = _select_any_finite_ba_slice(flux, magn, merr, cat, cfg['BEST_APER_ID'])

    fitted_line, inliers_mask = _get_sigma_clip_mask(magn, cat, flt_cname, out_dir)

    flux, magn, merr, cat = _select_inliers(flux, magn, merr, cat, inliers_mask)
    clc_magn = fitted_line(magn)
    clc_magn_med = np.nanmedian(clc_magn, axis=0)

    # _calc_magn_delta(clc_magn_med, cat, flt_cname, out_dir)
    # _calc_magn_std(magn, cat, flt_cname, out_dir)

    # if merr is not None:
    #     _calc_merr(clc_magn, merr, flt_cname, out_dir)

    # if flux is None:
    #     flux = 80 * np.power(10.0, -0.4 * magn)
    # _calc_total_throughput(flux, cat, flt_cname, out_dir)

    return clc_magn_med, cat


def rm_get_color_terms(r_clc_magn_med, i_clc_magn_med, r_cat, i_cat, cfg):
    print('')
    out_dir = cfg['OUT_DIR']

    inter, r_ids, i_ids = np.intersect1d(r_cat['ID'], i_cat['ID'], return_indices=True)
    r_clc_magn_med = r_clc_magn_med[r_ids]
    i_clc_magn_med = i_clc_magn_med[i_ids]
    cat = r_cat[r_ids]
    print(f"{len(cat)} stars in color calculation")

    cat_color = cat['rmag'] - cat['imag']
    r_mag_delta = cat['rmag'] - r_clc_magn_med
    i_mag_delta = cat['imag'] - i_clc_magn_med

    lin_fit_r = LinearLSQFitter()
    lin_fit_i = LinearLSQFitter()
    r_fitted_line = lin_fit_r(Linear1D(), cat_color, r_mag_delta)
    i_fitted_line = lin_fit_i(Linear1D(), cat_color, i_mag_delta)
    print(f"rmag color parameters {r_fitted_line.parameters}, {np.std(r_mag_delta - r_fitted_line(cat_color))}")
    print(f"imag color parameters {i_fitted_line.parameters}, {np.std(i_mag_delta - i_fitted_line(cat_color))}")

    _plot_color_term(cat_color, r_mag_delta, 'r', out_dir, r_fitted_line)
    _plot_color_term(cat_color, i_mag_delta, 'i', out_dir, i_fitted_line)

    # r_fin_magn_med = r_clc_magn_med + r_fitted_line(cat_color)
    # i_fin_magn_med = i_clc_magn_med + i_fitted_line(cat_color)

    clc_color = r_clc_magn_med - i_clc_magn_med

    lin_fit = LinearLSQFitter()
    alt_fitted_line = lin_fit(Linear1D(), clc_color, cat_color)
    print(f"alt color parameters {alt_fitted_line.parameters}, {np.std(cat_color - alt_fitted_line(clc_color))}")
    _plot_alt_color_term(clc_color, cat_color, out_dir, alt_fitted_line)


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

    # _plot_raw_magn(rb_magn_med, cat_magn, flt_cname[0], out_dir)

    # noinspection PyTypeChecker
    sc_fit = FittingWithOutlierRemoval(LinearLSQFitter(), sigma_clip, niter=5, sigma=3.0)
    fitted_line, outliers_mask = sc_fit(Linear1D(), rb_magn_med, cat_magn)
    inliers_mask = np.invert(outliers_mask)

    # TODO: Plot 1\sigma uncertainty
    # _plot_fitted_magn(rb_magn_med, cat_magn, flt_cname[0], out_dir, inliers_mask, fitted_line)
    print(f"{np.sum(inliers_mask)} stars, {sc_fit.fit_info['niter']} iterations, parameters {fitted_line.parameters},"
          f"std(residuals) {np.std(rb_magn_med - fitted_line(rb_magn_med))}")

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


def _calc_magn_delta(clc_magn_med, cat, flt_cname, out_dir):
    cat_magn = cat[flt_cname]
    magn_delta = cat_magn - clc_magn_med

    _plot_delta_magn(cat_magn, magn_delta, flt_cname[0], out_dir)
    _plot_check_magn(clc_magn_med, cat_magn, flt_cname[0], out_dir)


def _calc_magn_std(rb_magn, cat, flt_cname, out_dir):
    cat_magn = cat[flt_cname]

    rb_magn_std = np.nanstd(rb_magn, axis=0)

    otw_data = np.loadtxt(f"{out_dir}1.2M_errors_{flt_cname[0]}.txt")
    otw_magn = otw_data[:, 0]
    otw_std = otw_data[:, 1]

    _plot_std_magn(cat_magn, rb_magn_std, otw_magn, otw_std, flt_cname[0], out_dir)


def _calc_merr(clc_magn, rb_merr, flt_cname, out_dir):
    clc_magn_med = np.nanmedian(clc_magn, axis=0)
    rb_merr_med = np.nanmedian(rb_merr, axis=0)

    _plot_merr(clc_magn_med, rb_merr_med, flt_cname[0], out_dir)


def _calc_total_throughput(rb_flux, cat, flt_cname, out_dir):
    cat_magn = cat[flt_cname]
    rb_flux_med = np.nanmedian(rb_flux, axis=0)

    # flux [phot] = F_nu0 * 10^(-0.4 * m) / h * (Dlambda/lambda_eff) * S * T * P / G
    # P(RAper = 1 * FWHM) = 0.65
    match flt_cname:
        case 'rmag':
            cat_flux = (3631 * np.power(10, -0.4 * cat_magn) *
                        1.51e7 * (1253.71 / 6201.71) *
                        np.pi * 0.09 * 80 *
                        0.65 / 1.4)
        case 'imag':
            cat_flux = (3631 * np.power(10, -0.4 * cat_magn) *
                        1.51e7 * (1478.93 / 7672.59) *
                        np.pi * 0.09 * 80 *
                        0.65 / 1.4)
        case _:
            cat_flux = np.nan

    total_tp = rb_flux_med / cat_flux * 100
    _plot_total_throughput(cat_magn, total_tp, flt_cname[0], out_dir)
    print(f"TP in {flt_cname} {np.mean(total_tp)} +/- {np.std(total_tp)} %")

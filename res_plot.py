import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from filesys_io import check_out_dir


def _plot_raw_magn(rb_magn, cat_magn, flt, out_dir):
    fig, ax = plt.subplots(dpi=300)

    ax.plot(rb_magn, cat_magn, 'k.', markersize=3)

    ax.set_xlabel(f"$\\langle m_{{{flt}'RB}} \\rangle$ [mag]")
    ax.set_ylabel(f"$m_{{{flt}'CAT}}$ [mag]")
    ax.grid()

    fig.savefig(f'{out_dir}{flt}\\{flt}_magn_raw.png')


def _plot_fitted_magn(rb_magn, cat_magn, flt, out_dir, inliers_mask, fitted_line):
    fig, ax = plt.subplots(dpi=300)

    ax.plot(rb_magn[inliers_mask], cat_magn[inliers_mask], 'k.', markersize=3, label="Fitted data")
    ax.plot(rb_magn[~inliers_mask], cat_magn[~inliers_mask], 'r.', markersize=3, label="Clipped data")

    rb_magn_lims = [np.min(rb_magn), np.max(rb_magn)]
    ax.plot(rb_magn_lims, fitted_line(rb_magn_lims), 'g-', label="Fitted model")

    ax.set_xlabel(f"$\\langle m_{{{flt}'RB}} \\rangle$ [mag]")
    ax.set_ylabel(f"$m_{{{flt}'CAT}}$ [mag]")
    ax.grid()
    ax.legend()

    fig.savefig(f'{out_dir}{flt}\\{flt}_magn_fit.png')


def _plot_std_magn(cat_magn, rb_magn_std, otw_magn, otw_std, flt, out_dir):
    fig, ax = plt.subplots(dpi=300)

    ax.plot(cat_magn, rb_magn_std, 'k.', markersize=3, label="RoboPhot data, EXPTIME = 80 s")
    ax.plot(otw_magn, otw_std, 'r.', markersize=3, label="1.2m data, EXPTIME = 90 s")

    ax.set_xlabel(f"$m_{{{flt}'CAT}}$ [mag]")
    ax.set_ylabel(f"$\\sigma( \\langle m_{{{flt}'}} \\rangle )$ [mag]")
    ax.grid()
    ax.legend()

    fig.savefig(f'{out_dir}{flt}\\{flt}_std.png')


def _plot_merr(cat_magn, rb_merr, flt, out_dir):
    fig, ax = plt.subplots(dpi=300)

    ax.plot(cat_magn, rb_merr, 'k.', markersize=3)

    ax.set_xlabel(f"$m_{{{flt}'CAT}}$ [mag]")
    ax.set_ylabel(f"$\\langle merr_{{{flt}'RB}} \\rangle$ [mag]")
    ax.grid()

    fig.savefig(f'{out_dir}{flt}\\{flt}_merr.png')


def _plot_total_throughput(cat_magn, total_tp, flt, out_dir):
    fig, ax = plt.subplots(dpi=300)

    ax.plot(cat_magn, total_tp, 'k.', markersize=3)

    ax.set_xlabel(f"$m_{{{flt}'CAT}}$ [mag]")
    ax.set_ylabel(f"TP [%]")
    ax.grid()

    fig.savefig(f'{out_dir}{flt}\\{flt}_tp.png')


def _plot_color_term(cat_color, delta_magn, flt, out_dir, fitted_line):
    fig, ax = plt.subplots(dpi=300)

    ax.plot(cat_color, delta_magn, 'k.', markersize=3, label="Data")

    rb_magn_lims = [np.min(cat_color), np.max(cat_color)]
    ax.plot(rb_magn_lims, fitted_line(rb_magn_lims), 'r-', label="Fitted model")

    ax.set_xlabel(r"$m_{r'CAT} - m_{i'CAT}$ [mag]")
    ax.set_ylabel(f"$m_{{{flt}'CAT}} - \\langle m_{{{flt}'RB}} \\rangle$ [mag]")
    ax.grid()
    ax.legend()

    fig.savefig(f'{out_dir}{flt}\\{flt}_color_term.png')


def _plot_alt_color_term(cat_color, rb_color, out_dir, fitted_line):
    fig, ax = plt.subplots(dpi=300)

    ax.plot(cat_color, rb_color, 'k.', markersize=3, label="Data")

    rb_magn_lims = [np.min(cat_color), np.max(cat_color)]
    ax.plot(rb_magn_lims, fitted_line(rb_magn_lims), 'r-', label="Fitted model")

    ax.set_xlabel(r"$m_{r'CAT} - m_{i'CAT}$ [mag]")
    ax.set_ylabel(r"$m_{r'RB} - m_{i'RB}$ [mag]")
    ax.grid()
    ax.legend()

    fig.savefig(f'{out_dir}alt_color_term.png')


def _draw_sky_map(img_header, img_data, img_wcs, apertures, img_edge, out_dir):
    # Step 1: Prepare data
    img_data = np.log10(img_data)
    scs_mean, scs_median, scs_std = sigma_clipped_stats(img_data[img_edge: -1 - img_edge][img_edge: -1 - img_edge])
    dynr_max = scs_median + 2 * scs_std
    dynr_min = scs_median - 1 * scs_std

    # Step 2: Prepare figure
    fig, ax = plt.subplots(dpi=450, subplot_kw=dict(projection=img_wcs))

    ax.set_xlim(img_edge, img_header['NAXIS2'] - img_edge)
    ax.set_ylim(img_edge, img_header['NAXIS1'] - img_edge)

    if img_header['CD2_2'] < 0:
        ax.invert_xaxis()
    if img_header['CD1_1'] > 0:
        ax.invert_yaxis()

    # By default in WCSAxes (part of astropy.visualisation),
    # the tick and axis labels for the first coordinate (RA) are shown on the x-axis,
    # and the tick and axis labels for the second coordinate (DEC) are shown on the y-axis.
    # RoboPhot images are rotated by almost exactly 90 degrees, so we need to redefine this behaviour
    # For more info see https://github.com/astropy/astropy/issues/8521
    ax.coords[0].set_ticks_position('lr')
    ax.coords[0].set_ticklabel_position('lr')
    ax.coords[1].set_ticks_position('bt')
    ax.coords[1].set_ticklabel_position('bt')

    ax.coords[0].set_major_formatter('dd:mm:ss')
    ax.coords[1].set_major_formatter('dd:mm:ss')

    ax.coords[0].set_axislabel('DEJ2000, deg')
    ax.coords[1].set_axislabel('RAJ2000, deg')

    ax.coords.grid(color='blue', ls='--', alpha=0.7)
    # Title overlaps the tick labels
    # ax.set_title(f"{img_header['OBJNAME']}, F={img_header['FILTER']}, E={img_header['EXPTIME']}")

    # Step 3: Plot data
    ax.imshow(img_data, vmin=dynr_min, vmax=dynr_max, origin='lower', interpolation='nearest', cmap='gray_r')
    apertures.plot(ax=ax, color='red', lw=0.3, alpha=0.8)

    fig.savefig(f"{out_dir}map_{img_header['OBJNAME']}.png")


def _debug_extreme_ploting(flux, magn, merr, cat, cfg):
    check_out_dir(f".\\OUT\\DRFlux\\")
    check_out_dir(f".\\OUT\\DRMagn\\")
    check_out_dir(f".\\OUT\\DRMerr\\")

    flux_fig, flux_ax = plt.subplots(dpi=300)
    magn_fig, magn_ax = plt.subplots(dpi=300)
    merr_fig, merr_ax = plt.subplots(dpi=300)

    for _ in range(len(cat)):
        for __ in range(len(cfg['APER_RADII'])):
            flux_ax.plot(flux[__, :, _], label=f"R = {cfg['APER_RADII'][__]}")
            magn_ax.plot(magn[__, :, _], label=f"R = {cfg['APER_RADII'][__]}")
            merr_ax.plot(merr[__, :, _], label=f"R = {cfg['APER_RADII'][__]}")

        flux_ax.set_ylabel("Flux, ADU")
        magn_ax.set_ylabel("Magn, m")
        merr_ax.set_ylabel("Merr, m")
        flux_ax.grid()
        flux_ax.legend()
        magn_ax.grid()
        magn_ax.legend()
        merr_ax.grid()
        merr_ax.legend()

        flux_fig.savefig(f".\\OUT\\DRFlux\\{cat['ID'][_]}.png")
        magn_fig.savefig(f".\\OUT\\DRMagn\\{cat['ID'][_]}.png")
        merr_fig.savefig(f".\\OUT\\DRMerr\\{cat['ID'][_]}.png")
        flux_ax.cla()
        magn_ax.cla()
        merr_ax.cla()


def _debug_bool_nans_and_negs(flux, cat, cfg, prefix='unknw'):
    has_negatives = (flux < 0).any(axis=1)
    has_nans = np.isnan(flux).any(axis=1)

    check_out_dir(f".\\OUT\\{prefix}_debug\\")

    neg_fig, neg_ax = plt.subplots(dpi=300)
    nan_fig, nan_ax = plt.subplots(dpi=300)

    for _ in range(len(cfg['APER_RADII'])):
        neg_ax.plot(cat[f'{prefix}mag'], has_negatives[_], "r.", markersize=2)
        nan_ax.plot(cat[f'{prefix}mag'], has_nans[_], "r.", markersize=2)

        neg_ax.set_ylabel("Has negatives")
        neg_ax.set_xlabel(f'{prefix}mag, m')
        neg_ax.grid()
        nan_ax.set_ylabel("Has NaNs")
        nan_ax.set_xlabel(f'{prefix}mag, m')
        nan_ax.grid()

        neg_fig.savefig(f".\\OUT\\{prefix}_debug\\neg_{cfg['APER_RADII'][_]}.png")
        neg_ax.cla()
        nan_fig.savefig(f".\\OUT\\{prefix}_debug\\nan_{cfg['APER_RADII'][_]}.png")
        nan_ax.cla()

import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from filesys_io import check_out_dir


def _draw_sky_map(img_header, img_data, img_wcs, apertures, img_edge, out_dir):
    # Step 1: Prepare data
    img_data = np.log10(img_data)
    scs_mean, scs_median, scs_std = sigma_clipped_stats(img_data[img_edge: -1 - img_edge][img_edge: -1 - img_edge])
    dynr_max = scs_median + 2 * scs_std
    dynr_min = scs_median - 1 * scs_std

    # Step 2: Prepare figure
    fig, ax = plt.subplots(dpi=300, subplot_kw=dict(projection=img_wcs))

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


def _extreme_debug_plotting(flux, magn, merr, cat, cfg):
    check_out_dir(f".\\OUT\\DRFlux\\")
    check_out_dir(f".\\OUT\\DRMagn\\")
    check_out_dir(f".\\OUT\\DRMerr\\")

    flux_fig, flux_ax = plt.subplots(dpi=150)
    magn_fig, magn_ax = plt.subplots(dpi=150)
    merr_fig, merr_ax = plt.subplots(dpi=150)

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


def _bool__flux_debug(flux, cat, cfg):
    has_negatives = (flux < 0).any(axis=1)
    has_nans = np.isnan(flux).any(axis=1)

    check_out_dir(f".\\OUT\\DRBool\\")

    neg_fig, neg_ax = plt.subplots(dpi=150)
    nan_fig, nan_ax = plt.subplots(dpi=150)

    for _ in range(len(cfg['APER_RADII'])):
        neg_ax.plot(cat['imag'], has_negatives[_], "r.", markersize=2)
        nan_ax.plot(cat['imag'], has_nans[_], "r.", markersize=2)

        neg_ax.set_ylabel("Has negatives")
        neg_ax.set_xlabel("imag, m")
        neg_ax.grid()
        nan_ax.set_ylabel("Has NaNs")
        nan_ax.set_xlabel("imag, m")
        nan_ax.grid()

        neg_fig.savefig(f".\\OUT\\DRBool\\neg_{cfg['APER_RADII'][_]}.png")
        neg_ax.cla()
        nan_fig.savefig(f".\\OUT\\DRBool\\nan_{cfg['APER_RADII'][_]}.png")
        nan_ax.cla()

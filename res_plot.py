import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats


def _draw_sky_map(img_header, img_data, img_wcs, apertures, img_edge):
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

    fig.savefig(f"OUT\\map_{img_header['OBJNAME']}.png")


def plot_merr_graph(res_magn, res_merr):
    fig, ax = plt.subplots(dpi=150)

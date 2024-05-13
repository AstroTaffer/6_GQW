from filesys_io import *
from sel_images import rm_select_fits_files
from sel_sources import rm_select_sources
from phot_aper import rm_sources_aperture_photometry
from phot_diff import rm_sources_ensemble_photometry, _diff_phot_alt_flux_core
from res_get import rm_get_results, rm_get_color_terms


# Step 1: Read JSON config file
config = read_json_config('gqw_config.json')


# Step 2: FITS files selection
# r_ff_list = rm_select_fits_files(f"{config['IN_FITS_DIR']}r\\")
# i_ff_list = rm_select_fits_files(f"{config['IN_FITS_DIR']}i\\")
# write_ff_list(r_ff_list, f"{config['OUT_DIR']}r\\", prefix='r')
# write_ff_list(i_ff_list, f"{config['OUT_DIR']}i\\", prefix='i')
r_ff_list = read_ff_list(f"{config['OUT_DIR']}r\\", prefix='r')
i_ff_list = read_ff_list(f"{config['OUT_DIR']}i\\", prefix='i')


# Step 3: Sources selection
# r_src_cat = rm_select_sources(r_ff_list, config, 'rmag', config['CAT_RMAG_LIMIT'])
# i_src_cat = rm_select_sources(i_ff_list, config, 'imag', config['CAT_IMAG_LIMIT'])
# write_sel_src_cat(r_src_cat, f"{config['OUT_DIR']}r\\", prefix='r')
# write_sel_src_cat(i_src_cat, f"{config['OUT_DIR']}i\\", prefix='i')
r_src_cat = read_sel_src_cat(f"{config['OUT_DIR']}r\\", prefix='r')
i_src_cat = read_sel_src_cat(f"{config['OUT_DIR']}i\\", prefix='i')


# Step 4: Aperture photometry
# r_raw_flux, r_raw_magn, r_raw_merr = rm_sources_aperture_photometry(r_ff_list, r_src_cat, config, prefix='r')
# i_raw_flux, i_raw_magn, i_raw_merr = rm_sources_aperture_photometry(i_ff_list, i_src_cat, config, prefix='i')
# write_phot_res(r_raw_flux, r_raw_magn, r_raw_merr, f"{config['OUT_DIR']}r\\", prefix='r_raw')
# write_phot_res(i_raw_flux, i_raw_magn, i_raw_merr, f"{config['OUT_DIR']}i\\", prefix='i_raw')
# r_raw_flux, r_raw_magn, r_raw_merr = read_phot_res(f"{config['OUT_DIR']}r\\", prefix='r_raw')
# i_raw_flux, i_raw_magn, i_raw_merr = read_phot_res(f"{config['OUT_DIR']}i\\", prefix='i_raw')


# Step 5: Ensemble photometry
# Step 5.1: Simple flux-based equalizer
# r_alt_flux, r_alt_magn = _diff_phot_alt_flux_core(r_raw_flux)
# i_alt_flux, i_alt_magn = _diff_phot_alt_flux_core(i_raw_flux)
# write_phot_res(r_alt_flux, r_alt_magn, None, f"{config['OUT_DIR']}r\\", prefix='r_alt')
# write_phot_res(i_alt_flux, i_alt_magn, None, f"{config['OUT_DIR']}i\\", prefix='i_alt')
# r_alt_flux, r_alt_magn, r_alt_merr = read_phot_res(f"{config['OUT_DIR']}r\\", prefix='r_alt')
# i_alt_flux, i_alt_magn, i_alt_merr = read_phot_res(f"{config['OUT_DIR']}i\\", prefix='i_alt')

# Step 5.2 Simplified Astrokit
# r_smp_magn, r_smp_merr = rm_sources_ensemble_photometry(r_raw_magn, r_raw_merr, r_src_cat, config,
#                                                         'rmag', 'smp')
# i_smp_magn, i_smp_merr = rm_sources_ensemble_photometry(i_raw_magn, i_raw_merr, i_src_cat, config,
#                                                         'imag', 'smp')
# write_phot_res(None, r_smp_magn, r_smp_merr, f"{config['OUT_DIR']}r\\", prefix='r_smp')
# write_phot_res(None, i_smp_magn, i_smp_merr, f"{config['OUT_DIR']}i\\", prefix='i_smp')
# r_smp_flux, r_smp_magn, r_smp_merr = read_phot_res(f"{config['OUT_DIR']}r\\", prefix='r_smp')
# i_smp_flux, i_smp_magn, i_smp_merr = read_phot_res(f"{config['OUT_DIR']}i\\", prefix='i_smp')

# Step 5.2 Full Astrokit
# r_full_magn, r_full_merr = rm_sources_ensemble_photometry(r_raw_magn, r_raw_merr, r_src_cat, config,
#                                                           'rmag', 'full')
# i_full_magn, i_full_merr = rm_sources_ensemble_photometry(i_raw_magn, i_raw_merr, i_src_cat, config,
#                                                           'imag', 'full')
# write_phot_res(None, r_full_magn, r_full_merr, f"{config['OUT_DIR']}r\\", prefix='r_full')
# write_phot_res(None, i_full_magn, i_full_merr, f"{config['OUT_DIR']}i\\", prefix='i_full')
r_full_flux, r_full_magn, r_full_merr = read_phot_res(f"{config['OUT_DIR']}r\\", prefix='r_full')
i_full_flux, i_full_magn, i_full_merr = read_phot_res(f"{config['OUT_DIR']}i\\", prefix='i_full')


# Step 6: Fit and plot results
finite_mask_r = np.isfinite(r_full_magn[-1]).any(axis=0)
finite_mask_i = np.isfinite(i_full_magn[-1]).any(axis=0)
r_full_magn = r_full_magn[:, :, finite_mask_r]
i_full_magn = i_full_magn[:, :, finite_mask_i]
r_src_cat = r_src_cat[finite_mask_r]
i_src_cat = i_src_cat[finite_mask_i]

inter, r_ids, i_ids = np.intersect1d(r_src_cat['ID'], i_src_cat['ID'], return_indices=True)
r_full_magn = r_full_magn[:, :, r_ids]
i_full_magn = i_full_magn[:, :, i_ids]
cat = r_src_cat[r_ids]

r_clc_magn_med = np.nanmedian(r_full_magn[-1], axis=0)
i_clc_magn_med = np.nanmedian(i_full_magn[-1], axis=0)
rb_color = r_clc_magn_med - i_clc_magn_med

a = 0.035 * (0.21 - rb_color)
b = 0.041 * (0.21 - rb_color)
cat['rmag'] = cat['rmag'] + 0.035 * (0.21 - rb_color)
cat['imag'] = cat['imag'] + 0.041 * (0.21 - rb_color)


r_sel_magn, r_sel_cat = rm_get_results(None, r_full_magn, None, cat, config, 'rmag')
i_sel_magn, i_sel_cat = rm_get_results(None, i_full_magn, None, cat, config, 'imag')
rm_get_color_terms(r_sel_magn, i_sel_magn, r_sel_cat, i_sel_cat, config)


# Step -1: Debug
# 804 finite stars in aperture slice [-1]
# 743 stars, 3 iterations, parameters [ 0.98434301 20.97460636]
#
# 804 finite stars in aperture slice [-1]
# 730 stars, 3 iterations, parameters [ 0.98843148 20.65572853]
#
# rmag color parameters [-0.01570409  0.00314885]
# imag color parameters [-0.10649912  0.02881444]
# alt color parameters [0.90920497 0.02566559]
pass

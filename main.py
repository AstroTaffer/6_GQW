from filesys_io import *
from sel_images import rm_select_fits_files
from sel_sources import rm_select_sources
from phot_aper import rm_sources_aperture_photometry
from phot_diff import rm_sources_ensemble_photometry, _diff_phot_alt_flux_core
from res_get import rm_get_results


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
rm_get_results(r_full_flux, r_full_magn, r_full_merr, r_src_cat, config, 'rmag')
rm_get_results(i_full_flux, i_full_magn, i_full_merr, i_src_cat, config, 'imag')


# r_mrb_flux, r_mrb_magn, r_mrb_merr, r_fin_cat = rm_get_final_results(r_alt_flux, r_alt_magn, r_alt_merr, r_src_cat,
#                                                                      config['APER_RADII'], 2.5, 'rmag')
# i_mrb_flux, i_mrb_magn, i_mrb_merr, i_fin_cat = rm_get_final_results(i_alt_flux, i_alt_magn, i_alt_merr, i_src_cat,
#                                                                      config['APER_RADII'], 2.5, 'imag')
# write_phot_res(r_mrb_flux, r_mrb_magn, r_mrb_merr, f"{config['OUT_DIR']}r\\", prefix='r_fin')
# write_phot_res(i_mrb_flux, i_mrb_magn, i_mrb_merr, f"{config['OUT_DIR']}i\\", prefix='i_fin')
# write_sel_src_cat(r_src_cat, f"{config['OUT_DIR']}r\\", prefix='r_fin')
# write_sel_src_cat(i_src_cat, f"{config['OUT_DIR']}i\\", prefix='i_fin')
# r_mrb_flux, r_mrb_magn, r_mrb_merr = read_phot_res(f"{config['OUT_DIR']}r\\", prefix='r_fin')
# i_mrb_flux, i_mrb_magn, i_mrb_merr = read_phot_res(f"{config['OUT_DIR']}i\\", prefix='i_fin')
# r_src_cat = read_sel_src_cat(f"{config['OUT_DIR']}r\\", prefix='r_fin')
# i_src_cat = read_sel_src_cat(f"{config['OUT_DIR']}i\\", prefix='i_fin')


# Step 7: Results plotting
# rm_plot_results(r_alt_magn, None, r_src_cat, config, 'rmag', prefix='r_alt')
# rm_plot_results(i_alt_magn, None, i_src_cat, config, 'imag', prefix='i_alt')
# rm_plot_results(r_smp_magn, r_smp_merr, r_src_cat, config, 'rmag', prefix='r_smp')
# rm_plot_results(i_smp_magn, i_smp_merr, i_src_cat, config, 'imag', prefix='i_smp')
# rm_plot_results(r_full_magn, r_full_merr, r_src_cat, config, 'rmag', prefix='r_full')
# rm_plot_results(i_full_magn, i_full_merr, i_src_cat, config, 'imag', prefix='i_full')


# Step -1: Debug
pass
# _plot_magn(r_mrb_magn, r_fin_cat['rmag'], f"{config['OUT_DIR']}r_sg\\")
# _plot_magn(np.nanmedian(r_alt_magn[2], axis=0), r_src_cat['rmag'], f"{config['OUT_DIR']}r_sga\\")
# _plot_magn(i_mrb_magn, i_fin_cat['imag'], f"{config['OUT_DIR']}i_sg\\")
# _plot_magn(np.nanmedian(i_alt_magn[2], axis=0), i_src_cat['imag'], f"{config['OUT_DIR']}i_sga\\")
# _debug_bool_nans_and_negs(r_raw_flux, r_src_cat, config, prefix='r')
# _debug_bool_nans_and_negs(i_raw_flux, i_src_cat, config, prefix='i')
# _extreme_debug_plotting(raw_flux, raw_magn, raw_merr, sources_catalog, config)

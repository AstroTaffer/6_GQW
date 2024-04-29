from filesys_io import *
from sel_images import rm_select_fits_files
from sel_sources import rm_select_sources
from phot_aper import rm_sources_aperture_photometry, _aper_phot_altfull
from phot_diff import rm_sources_ensemble_photometry, _diff_phot_fluxcore
from res_plot import _debug_extreme_ploting, _debug_bool_nans_and_negs


# Step 1: Read JSON config file
config = read_json_config('gqw_config.json')


# Step 2: List FITS files
# fits_files_list = rm_select_fits_files(config)
# write_ff_list(fits_files_list, config['OUT_DIR'])
fits_files_list = read_ff_list(config['OUT_DIR'])


# Step 3: Select sources
# sources_catalog = rm_select_sources(fits_files_list, config)
# write_sel_src_cat(sources_catalog, config['OUT_DIR'])
sources_catalog = read_sel_src_cat(config['OUT_DIR'])


# Step 4: Aperture photometry
# raw_flux, raw_magn, raw_merr = rm_sources_aperture_photometry(fits_files_list, sources_catalog, config)
# write_phot_res(raw_flux, raw_magn, raw_merr, config['OUT_DIR'], prefix='raw')
# write_sel_src_cat(sources_catalog, config['OUT_DIR'])
raw_flux, raw_magn, raw_merr = read_phot_res(config['OUT_DIR'], prefix='raw')


# Step 5: Ensemble photometry
# alt_clr_flux, alt_clt_magn = _diff_phot_fluxcore(raw_flux)
# write_phot_res(alt_clr_flux, alt_clt_magn, None, config['OUT_DIR'], prefix='alt_clr')
clr_magn, clr_merr = rm_sources_ensemble_photometry(raw_magn, raw_merr, sources_catalog, config)
write_phot_res(None, clr_magn, clr_merr, config['OUT_DIR'], prefix='clr')
# clr_flux, clr_magn, clr_merr = read_phot_res(config['OUT_DIR'], prefix='clr')
# clr_flux, clr_magn, clr_merr = read_phot_res(config['OUT_DIR'], prefix='alt_clr')


# Step -1: Debug
# _bool__flux_debug(raw_flux, sources_catalog, config)
# _extreme_debug_plotting(raw_flux, raw_magn, raw_merr, sources_catalog, config)
# _aper_phot_altfull(sources_catalog)

from filesys_io import read_json_config
from sel_images import rm_select_fits_files
from sel_sources import rm_select_sources
from phot_aper import rm_sources_aperture_photometry
from phot_diff import rm_sources_ensemble_photometry
from res_plot import _extreme_debug_plotting

config = read_json_config('gqw_config.json')
fits_files_list = rm_select_fits_files(config)
sources_catalog = rm_select_sources(fits_files_list, config)
raw_flux, raw_magn, raw_merr = rm_sources_aperture_photometry(fits_files_list, sources_catalog, config)



clr_magn, clr_merr = rm_sources_ensemble_photometry(raw_magn, raw_merr, sources_catalog, config)


# _extreme_debug_plotting(raw_flux, raw_magn, raw_merr, sources_catalog, config)



# from diff_phot import rm_sources_ensemble_photometry
# from get_results import rm_get_results
#
#
# # clr_magn, clr_merr = rm_sources_ensemble_photometry(raw_magn, raw_merr, sources_cat, config)
# sources_cat, res_magn, res_merr = rm_get_results(raw_magn, raw_merr, sources_cat, config)

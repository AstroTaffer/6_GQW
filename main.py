from filesys_io import read_json_config
from fits_sel import rm_select_fits_files
from source_sel import rm_select_sources
from aper_phot import rm_sources_aperture_photometry
from diff_phot import rm_sources_ensemble_photometry


config = read_json_config('gqw_config.json')
fits_files_list = rm_select_fits_files(config)
sources_cat = rm_select_sources(fits_files_list, config)
raw_flux, raw_magn, raw_merr = rm_sources_aperture_photometry(fits_files_list, sources_cat, config)
rm_sources_ensemble_photometry(raw_magn, raw_merr, sources_cat, config)
from astropy.io import fits
import numpy as np
fields_list = ['gdn','gds','egs']#,
catalog_path = '/Users/lpr/Data/lirg_project/intake/CANDELS/catalog/JFang_CANDELS_Data/'
expath = '/Users/lpr/Data/lirg_project/output/catalog/'
for field in fields_list:
	catalog = fits.open(catalog_path+field+'_all.fits')[1].data
	band_counts = 0
	for column in catalog.columns:
		if 'FLUX_AUTO_' in column.name and not column.name.endswith('_08'): 
			band_counts += 1
	array_offset = np.full(band_counts,-999.)
	array_band = np.full(band_counts,'-999',dtype='U10')
	
	array_star = np.full([len(catalog),band_counts],-999.)
	lpr = 0
	col_counts = 0
	for column in catalog.columns:
		if 'FLUX_AUTO_' in column.name and not column.name.endswith('_08'):
			band = column.name[10:]
			print(column.name)
			for column2 in catalog.columns:
				if column2.name == 'ACS_' + band + '_FLUX' or column2.name == 'WFC3_' + band + '_FLUX':#egs field flux is lowercase
					print(column2.name)
					zp_band_original = input('Please input ' +field+' ' +band + ' zero-point:')
					zp_band_original = float(zp_band_original)
					ind = np.where((catalog['CLASS_STAR'] > 0.8) & ((zp_band_original-2.5*np.log10(catalog[column.name]) < 19)) & ((zp_band_original-2.5*np.log10(catalog[column.name]) > 15)))
					offset = np.nanmedian((zp_band_original - 2.5 * np.log10(catalog[ind][column.name])) - (23.9 - 2.5 * np.log10(catalog[ind][column2.name])))
					print(offset)
					array_offset[lpr] = offset
					array_band[lpr] = band
					corr = (zp_band_original - 2.5 * np.log10(catalog[ind][column.name])) - (23.9 - 2.5 * np.log10(catalog[ind][column2.name]))
					array_star[0:corr.shape[0],col_counts] = corr
					print(array_star.shape)
					lpr += 1
					col_counts += 1
	col0 = fits.Column(name='Band',array=array_band,format='A10')
	col1 = fits.Column(name='APER_TPHOT_ZP_OFFSET',array=array_offset,format='D')
	hdu = fits.BinTableHDU.from_columns([col0,col1])
	hdu.writeto(expath+field+'_zp_offset.fits',overwrite=True)
	# expath_name = expath+field+'_zp_corr.fits'
	# col = fits.Column(name='delta_'+array_band[0],array=array_star[:,0],format='D')
	# hdu = fits.BinTableHDU.from_columns([col])
	# hdu.writeto(expath_name,overwrite=True)
	# for num1 in range(1,band_counts):
	# 	col_name = 'delta_'+array_band[num1]
	# 	col = fits.Column(name=col_name,array=array_star[:,num1],format='D')
	# 	hdu_ori = fits.open(expath_name)[1]
	# 	hdu = fits.BinTableHDU.from_columns(hdu_ori.data.columns+col)
	# 	hdu.writeto(expath_name,overwrite=True)
	# col0 = fits.Column(name='delta_'+array_band[0],array=array_star[:,0],format='D')
	# col1 = fits.Column(name='delta_'+array_band[1],array=array_star[:,1],format='D')
	# col2 = fits.Column(name='delta_'+array_band[2],array=array_star[:,2],format='D')
	# col3 = fits.Column(name='delta_'+array_band[3],array=array_star[:,3],format='D')
	# col4 = fits.Column(name='delta_'+array_band[4],array=array_star[:,4],format='D')
	# col5 = fits.Column(name='delta_'+array_band[5],array=array_star[:,5],format='D')
	# col6 = fits.Column(name='delta_'+array_band[6],array=array_star[:,6],format='D')
	# col7 = fits.Column(name='delta_'+array_band[7],array=array_star[:,7],format='D')
	# col8 = fits.Column(name='delta_'+array_band[8],array=array_star[:,8],format='D')
	# hdu = fits.BinTableHDU.from_columns([col0,col1,col2,col3,col4,col5,col6,col7,col8])
	# expath_name = expath+field+'_zp_corr.fits'
	# hdu.writeto(expath_name,overwrite=True)	
# # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# # # ---------------------------------------------------------------------
# # # Use CANDELS aperture to auto photometry fraction multiply TPHOT
# # # photometry, to derive aperture flux.
# # # innerpart: aper_4 in egs and goodss, aper_5 in goodsn. because there is
# # # an offset between barroGOODSN photometry aperture and other fileds
# # # ---------------------------------------------------------------------
# # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# from astropy.io import fits
# import numpy as np
# fields_list = ['goodsn','goodss','egs']
# cata_path = '/Users/lpr/Data/lirg_project/output/catalog/'
# cata_suffix = '_Huangall_candels_radec_van_modifyz.fits'
# candels_path = '/Users/lpr/Data/lirg_project/intake/CANDELS/catalog/JFang_CANDELS_Data/'
# expath = '/Users/lpr/Data/lirg_project/output/catalog/'
# band_used = {'left1':['F606W','F125W','F814W','F160W'],'right1':['F606W','F125W','F125W','F160W']}
# inner = {'goodsn':['_APER_3','_APER_10'],'goodss':['_APER_4','_APER_11'],'egs':['_APER_4','_APER_11']}
# tphot = {'left1':['ACS','WFC3','ACS','WFC3'],'right1':['ACS','WFC3','WFC3','WFC3']}
# zp = 23.9
# for field in fields_list:
# 	cata = fits.open(cata_path+field+cata_suffix)[1].data
# 	if field == 'goodsn':
# 		cata_candels = fits.open(candels_path+'gdn_all.fits')[1].data
# 	elif field == 'goodss':
# 		cata_candels = fits.open(candels_path+'gds_all.fits')[1].data
# 	else:
# 		cata_candels = fits.open(candels_path+'egs_all.fits')[1].data
# 	array = np.full([len(cata),14],-999.)
# 	# array = np.full([len(cata),7],-999.)
# 	for num1 in range(0,len(cata)):
# 		idx = cata[num1]['ID']
# 		idx_candels = cata[num1]['ID_CANDELS']
# 		z = cata[num1]['z_used']
# 		if z > 0.8 and z < 1.:
# 			uv_filter_obs = band_used['left1'][0:2]
# 			uv_tphot_obs = tphot['left1'][0:2]
# 			vi_filter_obs = band_used['left1'][2:4]
# 			vi_tphot_obs = tphot['left1'][2:4]
# 		elif z > 1. and z < 1.3:
# 			uv_filter_obs = band_used['right1'][0:2]
# 			uv_tphot_obs = tphot['right1'][0:2]
# 			vi_filter_obs = band_used['right1'][2:4]
# 			vi_tphot_obs = tphot['right1'][2:4]
# 		# ========================= uv_u =========================
# 		# uv_u_img = fits.open(image_path+field+'/'+field+'_'+ uv_filter_obs[0] + '/'+field+'_' + uv_filter_obs[0] + '_597.fits')[0]
# 		# zp_uv_u_total = -2.5*np.log10(uv_u_img.header['PHOTFLAM'])-21.1-5*np.log10(uv_u_img.header['PHOTPLAM'])+18.692
# 		# for num2 in range(0,len(zp_corr_catalog)):
# 		# 	if zp_corr_catalog[num2][0] == uv_filter_obs[0].upper():
# 		# 		zp_corr_uv_u = zp_corr_catalog[num2][1]
# 		# uv_u_fraction = cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+uv_filter_obs[0]]/cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][1]+'_'+uv_filter_obs[0]]
# 		uv_u_fraction = cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+uv_filter_obs[0]]/cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX_AUTO_'+uv_filter_obs[0]] # use auto instead of aper_11
# 		print(uv_u_fraction)
# 		if field != 'egs':
# 			uv_u_bulge  = zp - 2.5*np.log10(uv_u_fraction*cata_candels[np.where(cata_candels['ID']==idx_candels)][uv_tphot_obs[0]+'_'+uv_filter_obs[0]+'_FLUX'])
# 			uv_u_disk = zp - 2.5*np.log10((1-uv_u_fraction)*cata_candels[np.where(cata_candels['ID']==idx_candels)][uv_tphot_obs[0]+'_'+uv_filter_obs[0]+'_FLUX'])
# 		else:
# 			uv_u_bulge  = zp - 2.5*np.log10(uv_u_fraction*cata_candels[np.where(cata_candels['ID']==idx_candels)][uv_tphot_obs[0]+'_'+uv_filter_obs[0]+'_flux'])
# 			uv_u_disk = zp - 2.5*np.log10((1-uv_u_fraction)*cata_candels[np.where(cata_candels['ID']==idx_candels)][uv_tphot_obs[0]+'_'+uv_filter_obs[0]+'_flux'])
# 		# uv_u_total = zp_uv_u_total - zp_corr_uv_u - 2.5*np.log10(cata[num1]['FLUX_AUTO_'+uv_filter_obs[0].upper()])
# 		# uv_u_bulge = zp_uv_u_total - zp_corr_uv_u - 2.5*np.log10(cata[num1]['FLUX_APER_4_'+uv_filter_obs[0].upper()])
# 		# uv_u_disk = zp_uv_u_total - zp_corr_uv_u - 2.5*np.log10(cata[num1]['FLUX_AUTO_'+uv_filter_obs[0].upper()] - cata[num1]['FLUX_APER_4_'+uv_filter_obs[0].upper()])
# 		# ========================= uv_v =========================
# 		# uv_v_img = fits.open(image_path+field+'/'+field+'_'+ uv_filter_obs[1] + '/'+field+'_' + uv_filter_obs[1] + '_597.fits')[0]
# 		# zp_uv_v_total = -2.5*np.log10(uv_v_img.header['PHOTFLAM'])-21.1-5*np.log10(uv_v_img.header['PHOTPLAM'])+18.692
# 		# for num2 in range(0,len(zp_corr_catalog)):
# 		# 	if zp_corr_catalog[num2][0] == uv_filter_obs[1].upper():
# 		# 		zp_corr_uv_v = zp_corr_catalog[num2][1]
# 		# uv_v_fraction = cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+uv_filter_obs[1]]/cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][1]+'_'+uv_filter_obs[1]]
# 		uv_v_fraction = cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+uv_filter_obs[1]]/cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX_AUTO_'+uv_filter_obs[1]] # use auto instead of aper_11
# 		print(uv_v_fraction)
# 		if field != 'egs':
# 			uv_v_bulge  = zp - 2.5*np.log10(uv_v_fraction*cata_candels[np.where(cata_candels['ID']==idx_candels)][uv_tphot_obs[1]+'_'+uv_filter_obs[1]+'_FLUX'])
# 			uv_v_disk = zp - 2.5*np.log10((1-uv_v_fraction)*cata_candels[np.where(cata_candels['ID']==idx_candels)][uv_tphot_obs[1]+'_'+uv_filter_obs[1]+'_FLUX'])
# 		else:
# 			uv_v_bulge  = zp - 2.5*np.log10(uv_v_fraction*cata_candels[np.where(cata_candels['ID']==idx_candels)][uv_tphot_obs[1]+'_'+uv_filter_obs[1]+'_flux'])
# 			uv_v_disk = zp - 2.5*np.log10((1-uv_v_fraction)*cata_candels[np.where(cata_candels['ID']==idx_candels)][uv_tphot_obs[1]+'_'+uv_filter_obs[1]+'_flux'])
# 		# uv_v_total = zp_uv_v_total - zp_corr_uv_v - 2.5*np.log10(cata[num1]['FLUX_AUTO_'+uv_filter_obs[1].upper()])
# 		# uv_v_bulge = zp_uv_v_total - zp_corr_uv_v - 2.5*np.log10(cata[num1]['FLUX_APER_4_'+uv_filter_obs[1].upper()])
# 		# uv_v_disk = zp_uv_v_total - zp_corr_uv_v - 2.5*np.log10(cata[num1]['FLUX_AUTO_'+uv_filter_obs[1].upper()] - cata[num1]['FLUX_APER_4_'+uv_filter_obs[1].upper()])
# 		# ========================= vi_v =========================
# 		# vi_v_img = fits.open(image_path+field+'/'+field+'_'+ vi_filter_obs[0] + '/'+field+'_' + vi_filter_obs[0] + '_597.fits')[0]
# 		# zp_vi_v_total = -2.5*np.log10(vi_v_img.header['PHOTFLAM'])-21.1-5*np.log10(vi_v_img.header['PHOTPLAM'])+18.692
# 		# for num2 in range(0,len(zp_corr_catalog)):
# 		# 	if zp_corr_catalog[num2][0] == vi_filter_obs[0].upper():
# 		# 		zp_corr_vi_v = zp_corr_catalog[num2][1]
# 		# vi_v_fraction = cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+vi_filter_obs[0]]/cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][1]+'_'+vi_filter_obs[0]]
# 		vi_v_fraction = cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+vi_filter_obs[0]]/cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX_AUTO_'+vi_filter_obs[0]] # use auto instead of aper_11
# 		print(vi_v_fraction)
# 		if field != 'egs':
# 			vi_v_bulge  = zp - 2.5*np.log10(vi_v_fraction*cata_candels[np.where(cata_candels['ID']==idx_candels)][vi_tphot_obs[0]+'_'+vi_filter_obs[0]+'_FLUX'])
# 			vi_v_disk = zp - 2.5*np.log10((1-vi_v_fraction)*cata_candels[np.where(cata_candels['ID']==idx_candels)][vi_tphot_obs[0]+'_'+vi_filter_obs[0]+'_FLUX'])
# 		else:
# 			vi_v_bulge  = zp - 2.5*np.log10(vi_v_fraction*cata_candels[np.where(cata_candels['ID']==idx_candels)][vi_tphot_obs[0]+'_'+vi_filter_obs[0]+'_flux'])
# 			vi_v_disk = zp - 2.5*np.log10((1-vi_v_fraction)*cata_candels[np.where(cata_candels['ID']==idx_candels)][vi_tphot_obs[0]+'_'+vi_filter_obs[0]+'_flux'])
# 		# vi_v_total = zp_vi_v_total - zp_corr_vi_v - 2.5*np.log10(cata[num1]['FLUX_AUTO_'+vi_filter_obs_copy.upper()])
# 		# vi_v_bulge = zp_vi_v_total - zp_corr_vi_v - 2.5*np.log10(cata[num1]['FLUX_APER_4_'+vi_filter_obs_copy.upper()])
# 		# vi_v_disk = zp_vi_v_total - zp_corr_vi_v - 2.5*np.log10(cata[num1]['FLUX_AUTO_'+vi_filter_obs_copy.upper()] - cata[num1]['FLUX_APER_4_'+vi_filter_obs_copy.upper()])
# 		# ========================= vi_i =========================
# 		# vi_i_img = fits.open(image_path+field+'/'+field+'_'+ vi_filter_obs[1] + '/'+field+'_' + vi_filter_obs[1] + '_597.fits')[0]
# 		# zp_vi_i_total = -2.5*np.log10(vi_i_img.header['PHOTFLAM'])-21.1-5*np.log10(vi_i_img.header['PHOTPLAM'])+18.692
# 		# for num2 in range(0,len(zp_corr_catalog)):
# 		# 	if zp_corr_catalog[num2][0] == vi_filter_obs[1].upper():
# 		# 		zp_corr_vi_i = zp_corr_catalog[num2][1]
# 		# vi_i_fraction = cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+vi_filter_obs[1]]/cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][1]+'_'+vi_filter_obs[1]]
# 		vi_i_fraction = cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+vi_filter_obs[1]]/cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX_AUTO_'+vi_filter_obs[1]] # use auto instead of aper_11
# 		print(vi_i_fraction)
# 		if field != 'egs':
# 			vi_i_bulge  = zp - 2.5*np.log10(vi_i_fraction*cata_candels[np.where(cata_candels['ID']==idx_candels)][vi_tphot_obs[1]+'_'+vi_filter_obs[1]+'_FLUX'])
# 			vi_i_disk = zp - 2.5*np.log10((1-vi_i_fraction)*cata_candels[np.where(cata_candels['ID']==idx_candels)][vi_tphot_obs[1]+'_'+vi_filter_obs[1]+'_FLUX'])
# 		else:
# 			vi_i_bulge  = zp - 2.5*np.log10(vi_i_fraction*cata_candels[np.where(cata_candels['ID']==idx_candels)][vi_tphot_obs[1]+'_'+vi_filter_obs[1]+'_flux'])
# 			vi_i_disk = zp - 2.5*np.log10((1-vi_i_fraction)*cata_candels[np.where(cata_candels['ID']==idx_candels)][vi_tphot_obs[1]+'_'+vi_filter_obs[1]+'_flux'])
# 		# vi_i_total = zp_vi_i_total - zp_corr_vi_i - 2.5*np.log10(cata[num1]['FLUX_AUTO_'+vi_filter_obs[1].upper()])
# 		# vi_i_bulge = zp_vi_i_total - zp_corr_vi_i - 2.5*np.log10(cata[num1]['FLUX_APER_4_'+vi_filter_obs[1].upper()])
# 		# vi_i_disk = zp_vi_i_total - zp_corr_vi_i - 2.5*np.log10(cata[num1]['FLUX_AUTO_'+vi_filter_obs[1].upper()] - cata[num1]['FLUX_APER_4_'+vi_filter_obs[1].upper()])
# 		# print(idx,uv_u_bulge,uv_v_bulge,uv_u_disk,uv_v_disk,vi_v_bulge,vi_i_bulge,vi_v_disk,vi_i_disk)
# 		array[num1] = [idx,uv_u_bulge-uv_v_bulge,uv_u_disk-uv_v_disk,vi_v_bulge-vi_i_bulge,vi_v_disk-vi_i_disk,uv_u_bulge,uv_v_bulge,uv_u_disk,uv_v_disk,vi_v_bulge,vi_i_bulge,vi_v_disk,vi_i_disk,z]
# 		# array[num1] = [idx,u_bulge,v_bulge,i_bulge,u_total,v_total,i_total]
# 		# print(array[num1])
# 	col0 = fits.Column(name='id_huang',array=array[:,0],format='K')
# 	col1 = fits.Column(name='uv_bulge_CANDELS_obs',array=array[:,1],format='D')
# 	col2 = fits.Column(name='uv_disk_CANDELS_obs',array=array[:,2],format='D')
# 	col3 = fits.Column(name='vi_bulge_CANDELS_obs',array=array[:,3],format='D')
# 	col4 = fits.Column(name='vi_disk_CANDELS_obs',array=array[:,4],format='D')
# 	col5 = fits.Column(name='uv_u_bulge_obs',array=array[:,5],format='D')
# 	col6 = fits.Column(name='uv_v_bulge_obs',array=array[:,6],format='D')
# 	col7 = fits.Column(name='uv_u_disk_obs',array=array[:,7],format='D')
# 	col8 = fits.Column(name='uv_v_disk_obs',array=array[:,8],format='D')
# 	col9 = fits.Column(name='vi_v_bulge_obs',array=array[:,9],format='D')
# 	col10 = fits.Column(name='vi_i_bulge_obs',array=array[:,10],format='D')
# 	col11 = fits.Column(name='vi_v_disk_obs',array=array[:,11],format='D')
# 	col12 = fits.Column(name='vi_i_disk_obs',array=array[:,12],format='D')
# 	col13 = fits.Column(name='z_used',array=array[:,13],format='D')
# 	hdu = fits.BinTableHDU.from_columns([col0,col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13])
# 	hdu.writeto(expath+field+'_uvi_obs_CANDELS_useauto.fits',overwrite=True)
# 	print('======================== '+field+' done =======================')




# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# # ---------------------------------------------------------------------
# # Use CANDELS aperture to auto photometry fraction multiply TPHOT
# # photometry, to derive aperture flux. This procedure only applied to
# # GOODS-N field, because only GOODS-N field have zero-point problem.
# # innerpart: aper_4 in egs and goodss, aper_5 in goodsn. because there is
# # an offset between barroGOODSN photometry aperture and other fileds
# # ---------------------------------------------------------------------
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
from astropy.io import fits
import numpy as np
fields_list = ['goodsn','goodss','egs']
cata_path = '/Users/lpr/Data/lirg_project/output/catalog/'
cata_suffix = '_Huangall_candels_radec_van_modifyz.fits'
candels_path = '/Users/lpr/Data/lirg_project/intake/CANDELS/catalog/JFang_CANDELS_Data/'
expath = '/Users/lpr/Data/lirg_project/output/catalog/'
band_used = {'left1':['F606W','F125W','F814W','F160W'],'right1':['F606W','F125W','F125W','F160W']}
inner = {'goodsn':['_APER_3','_APER_10'],'goodss':['_APER_4','_APER_11'],'egs':['_APER_4','_APER_11']}
tphot = {'left1':['ACS','WFC3','ACS','WFC3'],'right1':['ACS','WFC3','WFC3','WFC3']}
zp = 23.9
for field in fields_list:
	cata = fits.open(cata_path+field+cata_suffix)[1].data
	if field == 'goodsn':
		cata_candels = fits.open(candels_path+'gdn_all.fits')[1].data
	elif field == 'goodss':
		cata_candels = fits.open(candels_path+'gds_all.fits')[1].data
	else:
		cata_candels = fits.open(candels_path+'egs_all.fits')[1].data
	array = np.full([len(cata),14],-999.)
	# array = np.full([len(cata),7],-999.)
	for num1 in range(0,len(cata)):
		idx = cata[num1]['ID']
		idx_candels = cata[num1]['ID_CANDELS']
		z = cata[num1]['z_used']
		if z > 0.8 and z < 1.:
			uv_filter_obs = band_used['left1'][0:2]
			uv_tphot_obs = tphot['left1'][0:2]
			vi_filter_obs = band_used['left1'][2:4]
			vi_tphot_obs = tphot['left1'][2:4]
		elif z > 1. and z < 1.3:
			uv_filter_obs = band_used['right1'][0:2]
			uv_tphot_obs = tphot['right1'][0:2]
			vi_filter_obs = band_used['right1'][2:4]
			vi_tphot_obs = tphot['right1'][2:4]
		# ========================= uv_u =========================
		uv_u_fraction = cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+uv_filter_obs[0]]/cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][1]+'_'+uv_filter_obs[0]]
		print(uv_u_fraction)
		if field == 'goodsn':
			uv_u_bulge  = zp - 2.5*np.log10(uv_u_fraction*cata_candels[np.where(cata_candels['ID']==idx_candels)][uv_tphot_obs[0]+'_'+uv_filter_obs[0]+'_FLUX'])
			uv_u_disk = zp - 2.5*np.log10((1-uv_u_fraction)*cata_candels[np.where(cata_candels['ID']==idx_candels)][uv_tphot_obs[0]+'_'+uv_filter_obs[0]+'_FLUX'])
		else:
			uv_u_bulge  = zp - 2.5*np.log10(cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+uv_filter_obs[0]])
			uv_u_disk = zp - 2.5*np.log10(cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][1]+'_'+uv_filter_obs[0]]-cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+uv_filter_obs[0]])
		# ========================= uv_v =========================
		uv_v_fraction = cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+uv_filter_obs[1]]/cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][1]+'_'+uv_filter_obs[1]]
		print(uv_v_fraction)
		if field == 'goodsn':
			uv_v_bulge  = zp - 2.5*np.log10(uv_v_fraction*cata_candels[np.where(cata_candels['ID']==idx_candels)][uv_tphot_obs[1]+'_'+uv_filter_obs[1]+'_FLUX'])
			uv_v_disk = zp - 2.5*np.log10((1-uv_v_fraction)*cata_candels[np.where(cata_candels['ID']==idx_candels)][uv_tphot_obs[1]+'_'+uv_filter_obs[1]+'_FLUX'])
		else:
			uv_v_bulge  = zp - 2.5*np.log10(cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+uv_filter_obs[1]])
			uv_v_disk = zp - 2.5*np.log10(cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][1]+'_'+uv_filter_obs[1]]-cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+uv_filter_obs[1]])
		# ========================= vi_v =========================
		vi_v_fraction = cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+vi_filter_obs[0]]/cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][1]+'_'+vi_filter_obs[0]]
		print(vi_v_fraction)
		if field == 'goodsn':
			vi_v_bulge  = zp - 2.5*np.log10(vi_v_fraction*cata_candels[np.where(cata_candels['ID']==idx_candels)][vi_tphot_obs[0]+'_'+vi_filter_obs[0]+'_FLUX'])
			vi_v_disk = zp - 2.5*np.log10((1-vi_v_fraction)*cata_candels[np.where(cata_candels['ID']==idx_candels)][vi_tphot_obs[0]+'_'+vi_filter_obs[0]+'_FLUX'])
		else:
			vi_v_bulge  = zp - 2.5*np.log10(cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+vi_filter_obs[0]])
			vi_v_disk = zp - 2.5*np.log10(cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][1]+'_'+vi_filter_obs[0]]-cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+vi_filter_obs[0]])
		# ========================= vi_i =========================
		vi_i_fraction = cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+vi_filter_obs[1]]/cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][1]+'_'+vi_filter_obs[1]]
		print(vi_i_fraction)
		if field == 'goodsn':
			vi_i_bulge  = zp - 2.5*np.log10(vi_i_fraction*cata_candels[np.where(cata_candels['ID']==idx_candels)][vi_tphot_obs[1]+'_'+vi_filter_obs[1]+'_FLUX'])
			vi_i_disk = zp - 2.5*np.log10((1-vi_i_fraction)*cata_candels[np.where(cata_candels['ID']==idx_candels)][vi_tphot_obs[1]+'_'+vi_filter_obs[1]+'_FLUX'])
		else:
			vi_i_bulge  = zp - 2.5*np.log10(cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+vi_filter_obs[1]])
			vi_i_disk = zp - 2.5*np.log10(cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][1]+'_'+vi_filter_obs[1]]-cata_candels[np.where(cata_candels['ID']==idx_candels)]['FLUX'+inner[field][0]+'_'+vi_filter_obs[1]])
		array[num1] = [idx,uv_u_bulge-uv_v_bulge,uv_u_disk-uv_v_disk,vi_v_bulge-vi_i_bulge,vi_v_disk-vi_i_disk,uv_u_bulge,uv_v_bulge,uv_u_disk,uv_v_disk,vi_v_bulge,vi_i_bulge,vi_v_disk,vi_i_disk,z]
	col0 = fits.Column(name='id_huang',array=array[:,0],format='K')
	col1 = fits.Column(name='uv_bulge_CANDELS_obs',array=array[:,1],format='D')
	col2 = fits.Column(name='uv_disk_CANDELS_obs',array=array[:,2],format='D')
	col3 = fits.Column(name='vi_bulge_CANDELS_obs',array=array[:,3],format='D')
	col4 = fits.Column(name='vi_disk_CANDELS_obs',array=array[:,4],format='D')
	col5 = fits.Column(name='uv_u_bulge_obs',array=array[:,5],format='D')
	col6 = fits.Column(name='uv_v_bulge_obs',array=array[:,6],format='D')
	col7 = fits.Column(name='uv_u_disk_obs',array=array[:,7],format='D')
	col8 = fits.Column(name='uv_v_disk_obs',array=array[:,8],format='D')
	col9 = fits.Column(name='vi_v_bulge_obs',array=array[:,9],format='D')
	col10 = fits.Column(name='vi_i_bulge_obs',array=array[:,10],format='D')
	col11 = fits.Column(name='vi_v_disk_obs',array=array[:,11],format='D')
	col12 = fits.Column(name='vi_i_disk_obs',array=array[:,12],format='D')
	col13 = fits.Column(name='z_used',array=array[:,13],format='D')
	hdu = fits.BinTableHDU.from_columns([col0,col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13])
	hdu.writeto(expath+field+'_uvi_obs_CANDELS_onlygdnfrac.fits',overwrite=True)
	print('======================== '+field+' done =======================')



# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# ---------------------------------------------------------------------
# Do not use CANDELS aperture photometry and auto photometry.
# ---------------------------------------------------------------------
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# import matplotlib.pyplot as plt
# from lpr.image import photometry as pm
# from astropy.io import fits
# import numpy as np
# import os
# from scipy.integrate import quad
# import matplotlib.gridspec as gridspec
# from scipy.optimize import curve_fit
# from astropy.cosmology import FlatLambdaCDM as flcdm
# import matplotlib.colors as colors
# from matplotlib.patches import Ellipse, Circle
# import numpy as np
# from scipy.optimize import curve_fit
# from photutils.utils import calc_total_error as photoerr
# from photutils.aperture import CircularAnnulus as can
# from photutils.aperture import CircularAperture as cap
# from photutils.aperture import aperture_photometry as ap
# def kpc_per_arcsec(z):
#     angular_distance = flcdm(H0=70,Om0=0.3).angular_diameter_distance(z).to_value() # to_value transform the value with unit quantity to the value with unit float
#     arc_scale = angular_distance*np.pi*1000/(180*3600)
#     return(arc_scale)
# def rescale(arr):
#     arr_min = arr.min()
#     arr_max = arr.max()
#     return (arr - arr_min) / (arr_max - arr_min)
# def gaussian(x,mu,sigma):
#     return np.exp(-(x-mu)**2/(2*sigma**2))/(sigma*np.sqrt(2*np.pi))
# def r_determine(img,cx,cy,totalerr):
# 	for r in range(1,int(img.shape[0]/2)-1):
# 		aperture = can([cx,cy],r,r+1)
# 		aper_flux = ap(img,aperture,totalerr,method='exact')
# 		if aper_flux['aperture_sum'][0]/aper_flux['aperture_sum_err'][0]<10:
# 			break;
# 	return r
# print('------------------ HERE BEGINS THE CODE -----------------')
# cata_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/goodsn_Huangall_radecmatch_modifyz_kpc.fits'
# band_used = ['f606w','f850l','f125w','f160w']
# img_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/convolveGALAXYimage_gaussian/'
# short_wave_path = '/goodsn_f606w/'
# long_wave_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/convolveGALAXYimage_gaussian/goodsn_f125w/'
# expath = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/color_difference/'
# mask_path = '/Users/lpr/Data/SExtractor/goodsn_f125w/'
# wht_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/singleGALAXYimage_wht/'
# cata = fits.open(cata_path)[1].data
# # k-correction is valid when 160w cover the i band, which means zmax<=1.98-1
# cata = cata[np.where((cata['z_used']>0.8)&(cata['z_used']<1.3)&(cata['re_f160w']>0))]
# u = 365
# v = 551
# i = 806
# array = np.zeros([len(cata),12])
# wavelength_list = []
# for band in band_used:
# 	wavelength = int(band[1:-1])
# 	if wavelength < 200:
# 		wavelength *= 10
# 	wavelength_list.append(wavelength)
# wavelength_list = np.array(wavelength_list)
# for num1 in range(0,len(cata)):
# 	idx = cata[num1]['ID']
# 	print(idx)
# 	boo = True
# 	for num2 in band_used:
# 		if 'goodsn_'+num2+'_'+str(idx)+'.fits' not in os.listdir(img_path+'goodsn_'+num2+'_photutils/'):
# 			boo = False
# 	if boo:
# 		print('--------- every band has data of '+str(idx)+' ---------')
# 		z = cata[num1]['z_used']
# 		re = cata[num1]['re_f160w']
# 		img_160 = fits.open(img_path+'goodsn_f160w_photutils/'+'goodsn_f160w_'+str(idx)+'.fits')[0]
# 		img_wht_160 = fits.open(wht_path+'goodsn_f160w/'+'goodsn_'+str(idx)+'_f160w_wht.fits')[0]
# 		seg_cat = fits.open(mask_path+str(idx)+'_cat.fits')[1].data
# 		seg_img = fits.open(mask_path+str(idx)+'_seg.fits')[0].data
# 		distance = np.sqrt((seg_cat['X_IMAGE']-50)**2 + (seg_cat['Y_IMAGE']-50)**2)
# 		ind = np.where(distance==np.min(distance))[0]
# 		totalerr_160 = photoerr(img_160.data,np.sqrt(1/np.median(img_wht_160.data.flatten())),img_160.header['EXPTIME'])
# 		r_bulge = re * 0.5 / 0.06
# 		r_disk = r_determine(img_160.data,cx=50,cy=50,totalerr=totalerr_160)
# 		print('r_bulge: '+str(r_bulge)+'\n'+'r_disk: '+str(r_disk))
# 		# --------------------------------------------------
# 		if r_bulge < r_disk:
# 			u_z_hst_mag = []
# 			u_z_hst_magerr = []
# 			v_z_hst_mag = []
# 			v_z_hst_magerr = []
# 			i_z_hst_mag = []
# 			i_z_hst_magerr = []
# 			u_z = u * (1+z)
# 			v_z = v * (1+z)
# 			i_z = i * (1+z)
# 			band_used_copy = np.delete(band_used,3) # in case, v band and i band both use f160w
# 			wavelength_list_copy = np.delete(wavelength_list,3) # in case, v band and i band both use f160w
# 			print(u_z,v_z,i_z)
# 			# -------------------------------------- u band determined ------------------------------------------
# 			u_z_hst = band_used[np.argmin(abs(wavelength_list - u_z))]
# 			print('In observe-frame, I use '+u_z_hst+' for rest-frame U band.')
# 			# -------------------------------------- v band determined ------------------------------------------
# 			v_z_hst = band_used_copy[np.argmin(abs(wavelength_list_copy - v_z))]# v_z_hst = band_used[np.argmin(abs(wavelength_list - v_z))]
# 			print('In observe-frame, I use '+v_z_hst+' for rest-frame V band.')
# 			# -------------------------------------- i band determined ------------------------------------------
# 			i_z_hst = band_used[np.argmin(abs(wavelength_list - i_z))]
# 			if i_z_hst == v_z_hst:
# 				i_z_hst = 'f160w'
# 			if v_z_hst == 'f850l':
# 				u_z_hst = 'f606w'
# 			print('In observe-frame, I use '+i_z_hst+' for rest-frame I band.')
# 			img_u_z_hst = fits.open(img_path+'goodsn_'+u_z_hst+'_photutils/'+'goodsn_'+u_z_hst+'_'+str(idx)+'.fits')[0]
# 			img_v_z_hst = fits.open(img_path+'goodsn_'+v_z_hst+'_photutils/'+'goodsn_'+v_z_hst+'_'+str(idx)+'.fits')[0]
# 			img_i_z_hst = fits.open(img_path+'goodsn_'+i_z_hst+'_photutils/'+'goodsn_'+i_z_hst+'_'+str(idx)+'.fits')[0]
# 			img_u_z_hst_data = img_u_z_hst.data
# 			img_v_z_hst_data = img_v_z_hst.data
# 			img_i_z_hst_data = img_i_z_hst.data
# 			img_u_z_hst_wht = fits.open(wht_path+'goodsn_'+u_z_hst+'/goodsn_'+str(idx)+'_'+u_z_hst+'_wht.fits')[0]
# 			img_v_z_hst_wht = fits.open(wht_path+'goodsn_'+v_z_hst+'/goodsn_'+str(idx)+'_'+v_z_hst+'_wht.fits')[0]
# 			img_i_z_hst_wht = fits.open(wht_path+'goodsn_'+i_z_hst+'/goodsn_'+str(idx)+'_'+i_z_hst+'_wht.fits')[0]
# 			whether_do_photometry = True
# 			if img_u_z_hst_data.any() == 0 or img_v_z_hst_data.any() == 0 or img_i_z_hst_data.any() == 0 or img_u_z_hst_wht.data.any() == 0 or img_v_z_hst_wht.data.any() == 0 or img_i_z_hst_wht.data.any() == 0:
# 				whether_do_photometry = False
# 			if whether_do_photometry:
# 				u_z_hst_totalerr = photoerr(img_u_z_hst.data,np.sqrt(1/np.median(np.delete(img_u_z_hst_wht.data.flatten(),img_u_z_hst_wht.data.flatten() == 0))),img_u_z_hst.header['EXPTIME'])
# 				v_z_hst_totalerr = photoerr(img_v_z_hst.data,np.sqrt(1/np.median(np.delete(img_v_z_hst_wht.data.flatten(),img_v_z_hst_wht.data.flatten() == 0))),img_v_z_hst.header['EXPTIME'])
# 				i_z_hst_totalerr = photoerr(img_i_z_hst.data,np.sqrt(1/np.median(np.delete(img_i_z_hst_wht.data.flatten(),img_i_z_hst_wht.data.flatten() == 0))),img_i_z_hst.header['EXPTIME'])
# 				img_u_z_hst_data[np.where((seg_img.data!=ind+1)&(seg_img.data!=0))] = 0
# 				img_v_z_hst_data[np.where((seg_img.data!=ind+1)&(seg_img.data!=0))] = 0
# 				img_i_z_hst_data[np.where((seg_img.data!=ind+1)&(seg_img.data!=0))] = 0
# 				print('Here comes the photometry codes.')
# 				# ----------------------- u band photometry -----------------------
# 				aper_flux = ap(img_u_z_hst_data,cap([50,50],r_bulge),u_z_hst_totalerr,method='exact')
# 				u_z_hst_mag.append(aper_flux['aperture_sum'][0])
# 				u_z_hst_magerr.append(aper_flux['aperture_sum_err'][0])
# 				aper_flux = ap(img_u_z_hst_data,can([50,50],r_bulge,r_disk),u_z_hst_totalerr,method='exact')
# 				u_z_hst_mag.append(aper_flux['aperture_sum'][0])
# 				u_z_hst_magerr.append(aper_flux['aperture_sum_err'][0])
# 				u_zp = -2.5*np.log10(img_u_z_hst.header['PHOTFLAM'])-21.1-5*np.log10(img_u_z_hst.header['PHOTPLAM'])+18.692
# 				u_z_hst_magerr = abs(-2.5*np.array(u_z_hst_magerr)/(np.array(u_z_hst_mag)*np.log(10)))
# 				u_z_hst_mag = u_zp-2.5*np.log10(np.array(u_z_hst_mag))
# 				# ----------------------- v band photometry -----------------------
# 				aper_flux = ap(img_v_z_hst_data,cap([50,50],r_bulge),v_z_hst_totalerr,method='exact')
# 				v_z_hst_mag.append(aper_flux['aperture_sum'][0])
# 				v_z_hst_magerr.append(aper_flux['aperture_sum_err'][0])
# 				aper_flux = ap(img_v_z_hst_data,can([50,50],r_bulge,r_disk),v_z_hst_totalerr,method='exact')
# 				v_z_hst_mag.append(aper_flux['aperture_sum'][0])
# 				v_z_hst_magerr.append(aper_flux['aperture_sum_err'][0])
# 				v_zp = -2.5*np.log10(img_v_z_hst.header['PHOTFLAM'])-21.1-5*np.log10(img_v_z_hst.header['PHOTPLAM'])+18.692
# 				v_z_hst_magerr = abs(-2.5*np.array(v_z_hst_magerr)/(np.array(v_z_hst_mag)*np.log(10)))
# 				v_z_hst_mag = v_zp-2.5*np.log10(np.array(v_z_hst_mag))
# 				# ----------------------- i band photometry -----------------------
# 				aper_flux = ap(img_i_z_hst_data,cap([50,50],r_bulge),i_z_hst_totalerr,method='exact')
# 				i_z_hst_mag.append(aper_flux['aperture_sum'][0])
# 				i_z_hst_magerr.append(aper_flux['aperture_sum_err'][0])
# 				aper_flux = ap(img_i_z_hst_data,can([50,50],r_bulge,r_disk),i_z_hst_totalerr,method='exact')
# 				i_z_hst_mag.append(aper_flux['aperture_sum'][0])
# 				i_z_hst_magerr.append(aper_flux['aperture_sum_err'][0])
# 				i_zp = -2.5*np.log10(img_i_z_hst.header['PHOTFLAM'])-21.1-5*np.log10(img_i_z_hst.header['PHOTPLAM'])+18.692
# 				i_z_hst_magerr = abs(-2.5*np.array(i_z_hst_magerr)/(np.array(i_z_hst_mag)*np.log(10)))
# 				i_z_hst_mag = i_zp-2.5*np.log10(np.array(i_z_hst_mag))
# 				# ----------------------- write into array -----------------------
# 				uv_b,uv_berr = u_z_hst_mag[0]-v_z_hst_mag[0],np.sqrt(u_z_hst_magerr[0]**2+v_z_hst_magerr[0]**2)
# 				uv_d,uv_derr = u_z_hst_mag[1]-v_z_hst_mag[1],np.sqrt(u_z_hst_magerr[1]**2+v_z_hst_magerr[1]**2)
# 				vi_b,vi_berr = v_z_hst_mag[0]-i_z_hst_mag[0],np.sqrt(v_z_hst_magerr[0]**2+i_z_hst_magerr[0]**2)
# 				vi_d,vi_derr = v_z_hst_mag[1]-i_z_hst_mag[1],np.sqrt(v_z_hst_magerr[1]**2+i_z_hst_magerr[1]**2)
# 				array[num1] = [idx,uv_b,uv_berr,uv_d,uv_derr,vi_b,vi_berr,vi_d,vi_derr,int(u_z_hst[1:-1]),int(v_z_hst[1:-1]),int(i_z_hst[1:-1])]
# 		print('=========================== '+str(idx)+' is done ===========================')
# col1 = fits.Column(name='ID_Huang',array=array[:,0],format='K')
# col2 = fits.Column(name='UV_bulge',array=array[:,1],format='D')
# col3 = fits.Column(name='UV_bulgeerr',array=array[:,2],format='D')
# col4 = fits.Column(name='UV_disk',array=array[:,3],format='D')
# col5 = fits.Column(name='UV_diskerr',array=array[:,4],format='D')
# col6 = fits.Column(name='VI_bulge',array=array[:,5],format='D')
# col7 = fits.Column(name='VI_bulgeerr',array=array[:,6],format='D')
# col8 = fits.Column(name='VI_disk',array=array[:,7],format='D')
# col9 = fits.Column(name='VI_diskerr',array=array[:,8],format='D')
# col10 = fits.Column(name='U_HST_filter',array=array[:,9],format='K')
# col11 = fits.Column(name='V_HST_filter',array=array[:,10],format='K')
# col12 = fits.Column(name='I_HST_filter',array=array[:,11],format='K')
# hdu = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12])
# hdu.writeto(expath+'uvi_color_obs.fits',overwrite=True)
# print('=========================== GOODSN DONE ===========================')
# ||==|| The final color-color diagram showed two groups of galaxies.
# ||==|| Which may due to the filters used are too complicated.
# ||==|| Then, I tried the method used in wang+17. They used same filters for same
# ||==|| redshift bin.
# import matplotlib.pyplot as plt
# from lpr.image import photometry as pm
# from astropy.io import fits
# import numpy as np
# import os
# from scipy.integrate import quad
# import matplotlib.gridspec as gridspec
# from scipy.optimize import curve_fit
# from astropy.cosmology import FlatLambdaCDM as flcdm
# import matplotlib.colors as colors
# from matplotlib.patches import Ellipse, Circle
# import numpy as np
# from scipy.optimize import curve_fit
# from photutils.utils import calc_total_error as photoerr
# from photutils.aperture import CircularAnnulus as can
# from photutils.aperture import CircularAperture as cap
# from photutils.aperture import aperture_photometry as ap
# def kpc_per_arcsec(z):
#     angular_distance = flcdm(H0=70,Om0=0.3).angular_diameter_distance(z).to_value() # to_value transform the value with unit quantity to the value with unit float
#     arc_scale = angular_distance*np.pi*1000/(180*3600)
#     return(arc_scale)
# def rescale(arr):
#     arr_min = arr.min()
#     arr_max = arr.max()
#     return (arr - arr_min) / (arr_max - arr_min)
# def gaussian(x,mu,sigma):
#     return np.exp(-(x-mu)**2/(2*sigma**2))/(sigma*np.sqrt(2*np.pi))
# def r_determine(img,cx,cy,totalerr):
# 	for r in range(1,int(img.shape[0]/2)-1):
# 		aperture = can([cx,cy],r,r+1)
# 		aper_flux = ap(img,aperture,totalerr,method='exact')
# 		if aper_flux['aperture_sum'][0]/aper_flux['aperture_sum_err'][0]<10:
# 			break;
# 	return r

# print('------------------ HERE BEGINS THE CODE -----------------')
# fields_list = ['goodss','egs']#'goodsn',
# cata_path = '/Users/lpr/Data/lirg_project/output/catalog/'
# cata_suffix = '_Huangall_candels_radec_van_modifyz.fits'
# img_path = '/Users/lpr/Data/lirg_project/output/'
# # short_wave_path = '/goodsn_f606w/'
# # long_wave_path = '/Users/lpr/Data/lirg_project/output/'
# expath = '/Users/lpr/Data/lirg_project/output/catalog/'
# mask_path = '/Users/lpr/Data/SExtractor/'
# band_used = {'left1':['f606w','f125w','f814w','f160w'],'right1':['f606w','f125w','f125w','f160w']}
# for field in fields_list:
# 	zp_corr_catalog = fits.open('/Users/lpr/Data/lirg_project/output/catalog/'+field+'_zp_offset.fits')[1].data
# 	print(field)
# 	cata = fits.open(cata_path+field+cata_suffix)[1].data
# 	# k-correction is valid when 160w cover the i band, which means zmax<=1.98-1
# 	cata = cata[np.where((cata['z_used']>0.8)&(cata['z_used']<1.3)&(cata['re_f160w']>0.))]
# 	array = np.full([len(cata),26],-999.)
# 	for num1 in range(0,len(cata)):
# 		idx = cata[num1]['ID']
# 		boo = True
# 		z = cata[num1]['z_used']
# 		if z < 1.:
# 			for num2 in band_used['left1']:
# 				if field+'_'+num2+'_'+str(idx)+'_photutils.fits' not in os.listdir(img_path+field+'/'+field+'_'+num2) or field+'_'+num2+'_'+str(idx)+'_wht.fits' not in os.listdir(img_path+field+'/'+field+'_'+num2):
# 					boo = False
# 		elif z > 1.:
# 			for num2 in band_used['right1']:
# 				if field+'_'+num2+'_'+str(idx)+'_photutils.fits' not in os.listdir(img_path+field+'/'+field+'_'+num2) or field+'_'+num2+'_'+str(idx)+'_wht.fits' not in os.listdir(img_path+field+'/'+field+'_'+num2):
# 					boo = False
# 		if boo:
# 			print('--------- every band has data of '+str(idx)+' ---------')
# 			re = cata[num1]['re_f160w']
# 			img_160 = fits.open(img_path+field+'/'+field+'_f160w/'+field+'_f160w_'+str(idx)+'.fits')[0]
# 			img_wht_160 = fits.open(img_path+field+'/'+field+'_f160w/'+field+'_f160w_'+str(idx)+'_wht.fits')[0]
# 			seg_cat = fits.open(mask_path+field+'_f160w/'+str(idx)+'_cat.fits')[1].data
# 			seg_img = fits.open(mask_path+field+'_f160w/'+str(idx)+'_seg.fits')[0].data
# 			distance = np.sqrt((seg_cat['X_IMAGE']-50)**2 + (seg_cat['Y_IMAGE']-50)**2)
# 			ind = np.where(distance==np.min(distance))[0]
# 			totalerr_160 = photoerr(img_160.data,np.sqrt(1/np.median(img_wht_160.data.flatten())),img_160.header['EXPTIME'])
# 			if re < 0.19*2: # The resolution of HST in F160W is 0.19"
# 				r_bulge = 0.19
# 				r_disk = r_determine(img_160.data,cx=50,cy=50,totalerr=totalerr_160)
# 			else:
# 				r_bulge = re * 0.5 / 0.06
# 				r_disk = r_determine(img_160.data,cx=50,cy=50,totalerr=totalerr_160)
# 			print('r_bulge: '+str(r_bulge)+'\n'+'r_disk: '+str(r_disk))
# 			# --------------------------------------------------
# 			if r_bulge < r_disk :
# 				uv_u_mag = []
# 				uv_u_magerr = []
# 				uv_v_mag = []
# 				uv_v_magerr = []
# 				vi_v_mag = []
# 				vi_v_magerr = []
# 				vi_i_mag = []
# 				vi_i_magerr = []
# 				if z > 0.8 and z < 1.:
# 					uv_filter_obs = band_used['left1'][0:2]
# 					vi_filter_obs = band_used['left1'][2:4]
# 				elif z > 1. and z < 1.3:
# 					uv_filter_obs = band_used['right1'][0:2]
# 					vi_filter_obs = band_used['right1'][2:4]
# 				img_uv_u = fits.open(img_path+field+'/'+field+'_'+uv_filter_obs[0]+'/'+field+'_'+uv_filter_obs[0]+'_'+str(idx)+'_photutils.fits')[0]
# 				img_uv_v = fits.open(img_path+field+'/'+field+'_'+uv_filter_obs[1]+'/'+field+'_'+uv_filter_obs[1]+'_'+str(idx)+'_photutils.fits')[0]
# 				img_vi_v = fits.open(img_path+field+'/'+field+'_'+vi_filter_obs[0]+'/'+field+'_'+vi_filter_obs[0]+'_'+str(idx)+'_photutils.fits')[0]
# 				img_vi_i = fits.open(img_path+field+'/'+field+'_'+vi_filter_obs[1]+'/'+field+'_'+vi_filter_obs[1]+'_'+str(idx)+'.fits')[0]
# 				img_uv_u_data = img_uv_u.data
# 				img_uv_v_data = img_uv_v.data
# 				img_vi_v_data = img_vi_v.data
# 				img_vi_i_data = img_vi_i.data
# 				img_uv_u_wht = fits.open(img_path+field+'/'+field+'_'+uv_filter_obs[0]+'/'+field+'_'+uv_filter_obs[0]+'_'+str(idx)+'_wht.fits')[0]
# 				img_uv_v_wht = fits.open(img_path+field+'/'+field+'_'+uv_filter_obs[1]+'/'+field+'_'+uv_filter_obs[1]+'_'+str(idx)+'_wht.fits')[0]
# 				img_vi_v_wht = fits.open(img_path+field+'/'+field+'_'+vi_filter_obs[0]+'/'+field+'_'+vi_filter_obs[0]+'_'+str(idx)+'_wht.fits')[0]
# 				img_vi_i_wht = fits.open(img_path+field+'/'+field+'_'+vi_filter_obs[1]+'/'+field+'_'+vi_filter_obs[1]+'_'+str(idx)+'_wht.fits')[0]
# 				whether_do_photometry = True
# 				if img_uv_u.data.any() == 0 or img_uv_v.data.any() == 0 or img_vi_v.data.any() == 0 or img_vi_i.data.any() == 0 or img_uv_u_wht.data.any() == 0 or img_uv_v_wht.data.any() == 0 or img_vi_v_wht.data.any() == 0 or img_vi_i_wht.data.any() == 0:
# 					whether_do_photometry = False
# 				if whether_do_photometry:
# 					# ----------------------- uv color -----------------------
# 					uv_u_totalerr = photoerr(img_uv_u.data,np.sqrt(1/np.median(np.delete(img_uv_u_wht.data.flatten(),img_uv_u_wht.data.flatten() == 0))),img_uv_u.header['EXPTIME'])
# 					uv_v_totalerr = photoerr(img_uv_v.data,np.sqrt(1/np.median(np.delete(img_uv_v_wht.data.flatten(),img_uv_v_wht.data.flatten() == 0))),img_uv_v.header['EXPTIME'])
# 					vi_v_totalerr = photoerr(img_vi_v.data,np.sqrt(1/np.median(np.delete(img_vi_v_wht.data.flatten(),img_vi_v_wht.data.flatten() == 0))),img_vi_v.header['EXPTIME'])
# 					vi_i_totalerr = photoerr(img_vi_i.data,np.sqrt(1/np.median(np.delete(img_vi_i_wht.data.flatten(),img_vi_i_wht.data.flatten() == 0))),img_vi_i.header['EXPTIME'])
# 					img_uv_u_data[np.where((seg_img.data!=ind+1)&(seg_img.data!=0))] = 0
# 					img_uv_v_data[np.where((seg_img.data!=ind+1)&(seg_img.data!=0))] = 0
# 					img_vi_v_data[np.where((seg_img.data!=ind+1)&(seg_img.data!=0))] = 0
# 					img_vi_i_data[np.where((seg_img.data!=ind+1)&(seg_img.data!=0))] = 0
# 					print('Here comes the photometry codes.')
# 					# ----------------------- uv_u photometry -----------------------
# 					for num3 in range(0,len(zp_corr_catalog)):
# 						if zp_corr_catalog[num3][0] == uv_filter_obs[0].upper():
# 							zp_corr_uv_u = zp_corr_catalog[num3][1]
# 					aper_flux = ap(img_uv_u_data,cap([50,50],r_bulge),uv_u_totalerr,method='exact')
# 					uv_u_mag.append(aper_flux['aperture_sum'][0])
# 					uv_u_magerr.append(aper_flux['aperture_sum_err'][0])
# 					aper_flux = ap(img_uv_u_data,can([50,50],r_bulge,r_disk),uv_u_totalerr,method='exact')
# 					uv_u_mag.append(aper_flux['aperture_sum'][0])
# 					uv_u_magerr.append(aper_flux['aperture_sum_err'][0])
# 					uv_u_zp = -2.5*np.log10(img_uv_u.header['PHOTFLAM'])-21.1-5*np.log10(img_uv_u.header['PHOTPLAM'])+18.692
# 					print('uv_u_zp='+str(uv_u_zp))
# 					uv_u_magerr = abs(-2.5*np.array(uv_u_magerr)/(np.array(uv_u_mag)*np.log(10)))
# 					uv_u_mag = uv_u_zp-zp_corr_uv_u-2.5*np.log10(np.array(uv_u_mag))
# 					# ----------------------- uv_v photometry -----------------------
# 					for num3 in range(0,len(zp_corr_catalog)):
# 						if zp_corr_catalog[num3][0] == uv_filter_obs[1].upper():
# 							zp_corr_uv_v = zp_corr_catalog[num3][1]
# 					aper_flux = ap(img_uv_v_data,cap([50,50],r_bulge),uv_v_totalerr,method='exact')
# 					uv_v_mag.append(aper_flux['aperture_sum'][0])
# 					uv_v_magerr.append(aper_flux['aperture_sum_err'][0])
# 					aper_flux = ap(img_uv_v_data,can([50,50],r_bulge,r_disk),uv_v_totalerr,method='exact')
# 					uv_v_mag.append(aper_flux['aperture_sum'][0])
# 					uv_v_magerr.append(aper_flux['aperture_sum_err'][0])
# 					uv_v_zp = -2.5*np.log10(img_uv_v.header['PHOTFLAM'])-21.1-5*np.log10(img_uv_v.header['PHOTPLAM'])+18.692
# 					print('uv_v_zp='+str(uv_v_zp))
# 					uv_v_magerr = abs(-2.5*np.array(uv_v_magerr)/(np.array(uv_v_mag)*np.log(10)))
# 					uv_v_mag = uv_v_zp-zp_corr_uv_v-2.5*np.log10(np.array(uv_v_mag))
# 					# ----------------------- vi_v photometry -----------------------
# 					for num3 in range(0,len(zp_corr_catalog)):
# 						if zp_corr_catalog[num3][0] == vi_filter_obs[0].upper():
# 							zp_corr_vi_v= zp_corr_catalog[num3][1]
# 					aper_flux = ap(img_vi_v_data,cap([50,50],r_bulge),vi_v_totalerr,method='exact')
# 					vi_v_mag.append(aper_flux['aperture_sum'][0])
# 					vi_v_magerr.append(aper_flux['aperture_sum_err'][0])
# 					aper_flux = ap(img_vi_v_data,can([50,50],r_bulge,r_disk),vi_v_totalerr,method='exact')
# 					vi_v_mag.append(aper_flux['aperture_sum'][0])
# 					vi_v_magerr.append(aper_flux['aperture_sum_err'][0])
# 					vi_v_zp = -2.5*np.log10(img_vi_v.header['PHOTFLAM'])-21.1-5*np.log10(img_vi_v.header['PHOTPLAM'])+18.692
# 					print('vi_v_zp='+str(vi_v_zp))
# 					vi_v_magerr = abs(-2.5*np.array(vi_v_magerr)/(np.array(vi_v_mag)*np.log(10)))
# 					vi_v_mag = vi_v_zp-zp_corr_vi_v-2.5*np.log10(np.array(vi_v_mag))
# 					# ----------------------- vi_i photometry -----------------------
# 					for num3 in range(0,len(zp_corr_catalog)):
# 						if zp_corr_catalog[num3][0] == vi_filter_obs[1].upper():
# 							zp_corr_vi_i= zp_corr_catalog[num3][1]
# 					aper_flux = ap(img_vi_i_data,cap([50,50],r_bulge),vi_i_totalerr,method='exact')
# 					vi_i_mag.append(aper_flux['aperture_sum'][0])
# 					vi_i_magerr.append(aper_flux['aperture_sum_err'][0])
# 					aper_flux = ap(img_vi_i_data,can([50,50],r_bulge,r_disk),vi_i_totalerr,method='exact')
# 					vi_i_mag.append(aper_flux['aperture_sum'][0])
# 					vi_i_magerr.append(aper_flux['aperture_sum_err'][0])
# 					vi_i_zp = -2.5*np.log10(img_vi_i.header['PHOTFLAM'])-21.1-5*np.log10(img_vi_i.header['PHOTPLAM'])+18.692
# 					print('vi_i_zp='+str(vi_i_zp))
# 					vi_i_magerr = abs(-2.5*np.array(vi_i_magerr)/(np.array(vi_i_mag)*np.log(10)))
# 					vi_i_mag = vi_i_zp-zp_corr_vi_i-2.5*np.log10(np.array(vi_i_mag))
# 					# ----------------------- write into array -----------------------
# 					uv_b,uv_berr = uv_u_mag[0]-uv_v_mag[0],np.sqrt(uv_u_magerr[0]**2+uv_v_magerr[0]**2)
# 					uv_d,uv_derr = uv_u_mag[1]-uv_v_mag[1],np.sqrt(uv_u_magerr[1]**2+uv_v_magerr[1]**2)
# 					vi_b,vi_berr = vi_v_mag[0]-vi_i_mag[0],np.sqrt(vi_v_magerr[0]**2+vi_i_magerr[0]**2)
# 					vi_d,vi_derr = vi_v_mag[1]-vi_i_mag[1],np.sqrt(vi_v_magerr[1]**2+vi_i_magerr[1]**2)
# 					array[num1] = [idx,uv_b,uv_berr,uv_d,uv_derr,vi_b,vi_berr,vi_d,vi_derr,z,uv_u_mag[0],uv_u_magerr[0],uv_v_mag[0],uv_v_magerr[0],vi_v_mag[0],vi_v_magerr[0],vi_i_mag[0],vi_i_magerr[0],uv_u_mag[1],uv_u_magerr[1],uv_v_mag[1],uv_v_magerr[1],vi_v_mag[1],vi_v_magerr[1],vi_i_mag[1],vi_i_magerr[1]]
# 			print('=========================== '+str(idx)+' is done ===========================')
# 	col1 = fits.Column(name='ID_Huang',array=array[:,0],format='K')
# 	col2 = fits.Column(name='UV_bulge',array=array[:,1],format='D')
# 	col3 = fits.Column(name='UV_bulgeerr',array=array[:,2],format='D')
# 	col4 = fits.Column(name='UV_disk',array=array[:,3],format='D')
# 	col5 = fits.Column(name='UV_diskerr',array=array[:,4],format='D')
# 	col6 = fits.Column(name='VI_bulge',array=array[:,5],format='D')
# 	col7 = fits.Column(name='VI_bulgeerr',array=array[:,6],format='D')
# 	col8 = fits.Column(name='VI_disk',array=array[:,7],format='D')
# 	col9 = fits.Column(name='VI_diskerr',array=array[:,8],format='D')
# 	col10 = fits.Column(name='z_used',array=array[:,9],format='D')
# 	col11 = fits.Column(name='uv_u_bulge_mag',array=array[:,10],format='D')
# 	col12 = fits.Column(name='uv_u_bulge_magerr',array=array[:,11],format='D')
# 	col13 = fits.Column(name='uv_v_bulge_mag',array=array[:,12],format='D')
# 	col14 = fits.Column(name='uv_v_bulge_magerr',array=array[:,13],format='D')
# 	col15 = fits.Column(name='vi_v_bulge_mag',array=array[:,14],format='D')
# 	col16 = fits.Column(name='vi_v_bulge_magerr',array=array[:,15],format='D')
# 	col17 = fits.Column(name='vi_i_bulge_mag',array=array[:,16],format='D')
# 	col18 = fits.Column(name='vi_i_bulge_magerr',array=array[:,17],format='D')
# 	col19 = fits.Column(name='uv_u_disk_mag',array=array[:,18],format='D')
# 	col20 = fits.Column(name='uv_u_disk_magerr',array=array[:,19],format='D')
# 	col21 = fits.Column(name='uv_v_disk_mag',array=array[:,20],format='D')
# 	col22 = fits.Column(name='uv_v_disk_magerr',array=array[:,21],format='D')
# 	col23 = fits.Column(name='vi_v_disk_mag',array=array[:,22],format='D')
# 	col24 = fits.Column(name='vi_v_disk_magerr',array=array[:,23],format='D')
# 	col25 = fits.Column(name='vi_i_disk_mag',array=array[:,24],format='D')
# 	col26 = fits.Column(name='vi_i_disk_magerr',array=array[:,25],format='D')
# 	hdu = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25,col26])
# 	hdu.writeto(expath+field+'_uvi_obs.fits',overwrite=True)
# 	print('=========================== '+field+' DONE ===========================')
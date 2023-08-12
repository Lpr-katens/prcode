from astropy.io import fits
import numpy as np
from tqdm import tqdm

# catalog = ['goodsn_Huangall_candels_radec_van.fits','goodss_Huangall_candels_radec_van.fits','egs_Huangall_candels_radec_van.fits']#
candels_catalog = ['gdn_all.fits','gds_all.fits','egs_all.fits']
catalog = ['goodsn_Huangall_candels_van_params.fits','goodss_Huangall_candels_van_params.fits','egs_Huangall_candels_van_params.fits']
# # -|-|-|-|-|-|Firstly, add z_used to catalog.
# for num1 in tqdm(range(len(catalog))):
# 	hdu = fits.open('/Users/lpr/Data/lirg_project/output/catalog_radec/'+catalog[num1])[1].data
# 	array_column = np.full(len(hdu),-999.,dtype='float32')
# 	ind = np.where((hdu['ZSPEC1']>0.8)&(hdu['ZSPEC1']<1.3))
# 	array_column[ind] = hdu[ind]['ZSPEC1']
# 	ind = np.where(array_column==-999.)[0]
# 	ind = ind[np.where((hdu[ind]['col1']>0.8)&(hdu[ind]['col1']<1.3))]
# 	array_column[ind] = hdu[ind]['col1']
# 	col = fits.Column(name='z_used',array=array_column,format='D')
# 	hdu = fits.BinTableHDU.from_columns(hdu.columns+col)
# 	hdu.writeto('/Users/lpr/Data/lirg_project/output/catalog_radec/'+catalog[num1][:-5]+'_modifyz.fits',overwrite=True)
# # -|-|-|-|-|-|Secondly, add kpc of re to catalog.
# from astropy.cosmology import FlatLambdaCDM as flcdm
# # Define a function to calculate kpc/"
# def kpc_per_arcsec(z):
#     angular_distance = flcdm(H0=70,Om0=0.3).angular_diameter_distance(z).to_value()
#     arc_scale = angular_distance*np.pi*1000/(180*3600)
#     return(arc_scale)
# for num1 in range(0,len(catalog)):
# 	hdu = fits.open('/Users/lpr/Data/lirg_project/output/catalog_radec/'+catalog[num1][:-5]+'_modifyz.fits')[1].data
# 	kpc = np.full(len(hdu),-999.,dtype='float32')
# 	for num2 in range(0,len(hdu)):
# 		if hdu[num2]['z_used'] > 0.8 and hdu[num2]['z_used'] < 1.3 and hdu[num2]['re_f160w'] > 0:
# 			kpc[num2] = hdu[num2]['re_f160w']*kpc_per_arcsec(hdu[num2]['z_used'])
# 	col = fits.Column(name='re_kpc',array=kpc,format='D')
# 	hdu = hdu = fits.BinTableHDU.from_columns(hdu.columns+col)
# 	hdu.writeto('/Users/lpr/Data/lirg_project/output/catalog_radec/'+catalog[num1][:-5]+'_kpc.fits',overwrite=True)
# -|-|-|-|-|-|Thirdly, add SFR with unit Msun/yr.
# ================== modify mass and mass error of sample catalog
# for num1 in tqdm(range(len(catalog))):
# 	hdu = fits.open('/Users/lpr/Data/lirg_project/output/catalog_radec/'+catalog[num1][:-5]+'_modifyz.fits')[1].data
# 	candels = fits.open('/Users/lpr/Data/lirg_project/intake/CANDELS/catalog/JFang_CANDELS_Data/'+candels_catalog[num1])[1].data
# 	candels = candels[np.isin(candels['id'],hdu['id_candels'])]
# 	lmass_candels = np.full(len(hdu),-999.,dtype='float32')
# 	lmass_candels_err = np.full(len(hdu),-999.,dtype='float32')
# 	for num2 in range(len(hdu)):
# 		if hdu[num2]['z_used'] > 0:
# 			idx_candels = hdu[num2]['ID_CANDELS']
# 			lmass_candels[num2] = candels[np.where(candels['ID']==idx_candels)]['M_med']
# 			lmass_candels_err[num2] = candels[np.where(candels['ID']==idx_candels)]['s_med']
# 	col0 = fits.Column(name='lmass_candels',array=lmass_candels,format='D')
# 	col1 = fits.Column(name='lmass_candelserr',array=lmass_candels_err,format='D')
# 	hdu = fits.BinTableHDU.from_columns(hdu.columns+col0+col1)
# 	hdu.writeto('/Users/lpr/Data/lirg_project/output/catalog_radec/'+catalog[num1][:-5]+'_zmass.fits',overwrite=True)
# ================= use whole 8um luminosity to calculate sfr
# for num1 in range(0,len(catalog)):
# 	hdu = fits.open('/Users/lpr/Data/lirg_project/output/catalog_radec/'+catalog[num1][:-5]+'_zmass.fits')[1].data
# 	sfr = np.full(len(hdu),-999.,dtype='float32')
# 	sfr_err = np.full(len(hdu),-999.,dtype='float32')
# 	for num2 in range(0,len(hdu)):
# 		if hdu[num2]['z_used'] > 0:
# 			sfr[num2] = ((10**hdu[num2]['L8'][0])*4*1.73*1e-10)*0.63 # multiply 0.63 is for transformimg from salpeter IMF to chabrier IMF
# 			sfr_err[num2] = (10**hdu[num2]['L8'][0] * np.log(10) * hdu[num2]['E8'][0] * 4*1.73*1e-10)*0.63
# 	col0 = fits.Column(name='SFR_8um',array=sfr,format='D')
# 	col1 = fits.Column(name='SFR_8um_err',array=sfr_err,format='D')
# 	hdu = fits.BinTableHDU.from_columns(hdu.columns+col0+col1)
# 	hdu.writeto('/Users/lpr/Data/lirg_project/output/catalog_radec/'+catalog[num1][:-5]+'_zmasssfr.fits',overwrite=True)
# ================= use 8um star forming component to calculate sfr
# for num1 in range(0,len(catalog)):
# 	hdu = fits.open('/Users/lpr/Data/lirg_project/output/catalog_radec/'+catalog[num1][:-5]+'_zmasssfr.fits')[1].data
# 	sfr = np.full(len(hdu),-999.,dtype='float32')
# 	sfr_err = np.full(len(hdu),-999.,dtype='float32')
# 	for num2 in range(0,len(hdu)):
# 		if hdu[num2]['z_used'] > 0:
# 			sfr[num2] = (hdu[num2]['L8SFR']*4*1.73*1e-10)*0.63 # multiply 0.63 is for transformimg from salpeter IMF to chabrier IMF
# 			sfr_err[num2] = (hdu[num2]['E8SFR'] * 4*1.73*1e-10)*0.63
# 	col0 = fits.Column(name='SFR_8umsfr',array=sfr,format='D')
# 	col1 = fits.Column(name='SFR_8umsfr_err',array=sfr_err,format='D')
# 	hdu = fits.BinTableHDU.from_columns(hdu.columns+col0+col1)
# 	hdu.writeto('/Users/lpr/Data/lirg_project/output/catalog_radec/'+catalog[num1][:-5]+'_zmasssfr8sfr.fits',overwrite=True)
# # ================= use 8um star forming component for l8sf/l8 >0.85 to calculate sfr and use lfir for l8sf/l8 <0.85 to calculate sfr
for num1 in tqdm(range(len(catalog))):
	hdu = fits.open('/Users/lpr/Data/lirg_project/output/catalog_radec/'+catalog[num1][:-5]+'_zmass.fits')[1].data
	sfr = np.full(len(hdu),-999.,dtype='float32')
	sfr_err = np.full(len(hdu),-999.,dtype='float32')
	flag = np.full(len(hdu),'8umSFR',dtype='U6')
	for num2 in range(len(hdu)):
		if hdu[num2]['z_used'] > 0:
			if hdu[num2]['L8SFR']/10.**hdu[num2]['L8'][0] >= 0.5: # 0.5的分界线是从Huang21 sec4.2中来的
				sfr[num2] = (hdu[num2]['L8SFR']*4*1.73*1e-10)*0.63 # 0.63 is the correction from salpeter IMF to chabrier IMF
				sfr_err[num2] = (hdu[num2]['E8SFR'] * 4*1.73*1e-10)*0.63
				flag[num2] = '8umSFR'
			elif hdu[num2]['L8SFR']/10**hdu[num2]['L8'][0] < 0.5:
				sfr[num2] = (10**hdu[num2]['LIR_FULL']*4.5*3.826*1e-11)*0.63 # 4.5*1e-44 come form kennicuut+98
				sfr_err[num2] = (10**hdu[num2]['LIR_FULL']*np.log(10)*hdu[num2]['LERR_FULL']*4.5*3.826*1e-11)*0.63
				flag[num2] = 'FIRSFR'
	col0 = fits.Column(name='SFR_huang',array=sfr,format='D')
	col1 = fits.Column(name='SFR_huangerr',array=sfr_err,format='D')
	col2 = fits.Column(name='SFR_Flag',array=flag,format='6A')
	hdu = fits.BinTableHDU.from_columns(hdu.columns+col0+col1+col2)
	hdu.writeto('/Users/lpr/Data/lirg_project/output/catalog_radec/'+catalog[num1][:-5]+'_zmasshuangsfr.fits',overwrite=True)
# ================== match morphological classification catalog to sample catalog
# from astropy.io import fits
# import numpy as np
# catalog = ['goodss_all/goodss_Huangall_radecmatch_modifyz.fits','goodsn_all/goodsn_Huangall_radecmatch_modifyz.fits','aegis_all/aegis_Huangall_radecmatch_modifyz.fits']
# for num1 in range(0,len(catalog)):
# 	hdu = fits.open('/Users/lpr/Data/fits/expdata/HST/'+catalog[num1])[1].data
# 	morph = np.full(len(hdu),-1,dtype='h')
# 	col = fits.Column(name='Morph_class',array=morph,format='K')
# 	hdu = fits.BinTableHDU.from_columns(hdu.columns+col)
# 	hdu.writeto('/Users/lpr/Data/fits/expdata/HST/'+catalog[num1][:-5]+'_morph.fits',overwrite=True)
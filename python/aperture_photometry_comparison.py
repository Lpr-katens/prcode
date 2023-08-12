import matplotlib.pyplot as plt
# from lpr.image import photometry as pm
from astropy.io import fits
import numpy as np
import os
# from scipy.integrate import quad
# import matplotlib.gridspec as gridspec
# from scipy.optimize import curve_fit
# from astropy.cosmology import FlatLambdaCDM as flcdm
# import matplotlib.colors as colors
# from matplotlib.patches import Ellipse, Circle
import numpy as np
from photutils.utils import calc_total_error as photoerr
from photutils.aperture import CircularAnnulus as can
from photutils.aperture import CircularAperture as cap
from photutils.aperture import aperture_photometry as ap
# def kpc_per_arcsec(z):
#     angular_distance = flcdm(H0=70,Om0=0.3).angular_diameter_distance(z).to_value()
#     # to_value transform the value with unit quantity to the value with unit float
#     arc_scale = angular_distance*np.pi*1000/(180*3600)
#     return(arc_scale)
# def rescale(arr):
#     arr_min = arr.min()
#     arr_max = arr.max()
#     return (arr - arr_min) / (arr_max - arr_min)
# def gaussian(x,mu,sigma):
#     return np.exp(-(x-mu)**2/(2*sigma**2))/(sigma*np.sqrt(2*np.pi))
fields_list = ['goodsn','goodss','egs']    
image_path = '/Users/lpr/Data/lirg_project/output/'
cata_path = '/Users/lpr/Data/lirg_project/output/catalog/'
candels_cata_path = '/Users/lpr/Data/lirg_project/intake/CANDELS/catalog/JFang_CANDELS_Data/'
cata_suffix = '_Huangall_radec_candels.fits'
aper_list = [0.088,0.125,0.176,0.25,0.35,0.5,0.71,1,1.414,2,2.828] # arcsec"
x,y = 50,50
for field in fields_list:
	cata = fits.open(cata_path+field+cata_suffix)[1].data
	band_list = []
	if field == 'goodsn':
		candels_cata = fits.open(candels_cata_path+'gdn_all.fits')[1].data
		for column in candels_cata.columns:
			if 'FLUX_APER_1_' in column.name:
				if field+'_'+column.name[column.name.rfind('_')+1:][0:5].lower() in os.listdir(image_path+field):
					band_list.append(column.name[column.name.rfind('_',)+1:])
	elif field == 'goodss':
		candels_cata = fits.open(candels_cata_path+'gds_all.fits')[1].data
		for column in candels_cata.columns:
			if 'FLUX_APER_1_' in column.name:
				if field+'_'+column.name[column.name.rfind('_')+1:][0:5].lower() in os.listdir(image_path+field):
					band_list.append(column.name[column.name.rfind('_',)+1:])
	else:
		candels_cata = fits.open(candels_cata_path+'egs_all.fits')[1].data
		for column in candels_cata.columns:
			if 'FLUX_APER_1_' in column.name:
				if field+'_'+column.name[column.name.rfind('_')+1:][0:5].lower() in os.listdir(image_path+field):
					band_list.append(column.name[column.name.rfind('_',)+1:])
	array = np.full([len(cata),len(band_list)*len(aper_list)],-999.,dtype='d')
	idx_col = np.full([len(cata),1],-1)
	for num1 in range(0,len(cata)):
		idx = cata[num1]['ID']
		temp_row = []
		for band in band_list:
			if band!='f160w'.upper():
				if field+'_'+band[0:5].lower()+'_'+str(idx)+'_photutils.fits' in os.listdir(image_path+field+'/'+field+'_'+band[0:5].lower()):
					img = fits.open(image_path+field+'/'+field+'_'+band[0:5].lower()+'/'+field+'_'+band[0:5].lower()+'_'+str(idx)+'_photutils.fits')[0]
					zp = -2.5*np.log10(img.header['PHOTFLAM'])-21.1-5*np.log10(img.header['PHOTPLAM'])+18.692
					for aper in aper_list:
						aper_pixels = aper/(img.header['CD2_2']*3600)
						aperture = cap([x,y],aper_pixels)
						aper_flux = ap(img.data,aperture,method='exact')
						mag_band_aper = zp-2.5*np.log10(aper_flux['aperture_sum'][0])
						temp_row.append(mag_band_aper)
				elif field+'_'+band[0:5].lower()+'_'+str(idx)+'.fits' in os.listdir(image_path+field+'/'+field+'_'+band[0:5].lower()):
					img = fits.open(image_path+field+'/'+field+'_'+band[0:5].lower()+'/'+field+'_'+band[0:5].lower()+'_'+str(idx)+'.fits')[0]
					zp = -2.5*np.log10(img.header['PHOTFLAM'])-21.1-5*np.log10(img.header['PHOTPLAM'])+18.692
					for aper in aper_list:
						aper_pixels = aper/(img.header['CD2_2']*3600)
						aperture = cap([x,y],aper_pixels)
						aper_flux = ap(img.data,aperture,method='exact')
						mag_band_aper = zp-2.5*np.log10(aper_flux['aperture_sum'][0])
						temp_row.append(mag_band_aper)
						# print(mag_band_aper)
			elif band=='f160w'.upper():
				if field+'_'+band[0:5].lower()+'_'+str(idx)+'.fits' in os.listdir(image_path+field+'/'+field+'_'+band[0:5].lower()):
					img = fits.open(image_path+field+'/'+field+'_'+band[0:5].lower()+'/'+field+'_'+band[0:5].lower()+'_'+str(idx)+'.fits')[0]
					zp = -2.5*np.log10(img.header['PHOTFLAM'])-21.1-5*np.log10(img.header['PHOTPLAM'])+18.692
					for aper in aper_list:
						aper_pixels = aper/(img.header['CD2_2']*3600)
						aperture = cap([x,y],aper_pixels)
						aper_flux = ap(img.data,aperture,method='exact')
						mag_band_aper = zp-2.5*np.log10(aper_flux['aperture_sum'][0])
						temp_row.append(mag_band_aper)
		if np.array(temp_row).shape[0] == array.shape[1]:
			print('-------------------------- '+str(idx)+' is done --------------------------')
			idx_col[num1] = idx
			array[num1] = np.array(temp_row)
	col = fits.Column(name='ID_Huang',array=idx_col,format='K')
	hdu = fits.BinTableHDU.from_columns([col])
	hdu.writeto(cata_path+field+'_aper_mag_lpr.fits',overwrite=True)
	for num2 in range(0,len(band_list)):
		for num3 in range(0,len(aper_list)):
			col_name = 'flux_aper_'+str(num3+1)+'_'+band_list[num2].lower()
			col = fits.Column(name=col_name,array=array[:,num2*len(aper_list)+num3],format='D')
			hdu_ori = fits.open(cata_path+field+'_aper_mag_lpr.fits')[1]
			hdu = fits.BinTableHDU.from_columns(hdu_ori.data.columns+col)
			hdu.writeto(cata_path+field+'_aper_mag_lpr.fits',overwrite=True)
	print('=========================== '+field+' DONE ===========================')	
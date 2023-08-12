# The first idea to evaluate the fitting result is to calculate the average residual inside re_f160w.
# The second idea to evaluate the fitting result is to calculate residual/original_image.
# Then recieve the fitting result whose residual/original_image<30%.
# ---------------------------------------------------------------------------------------------------
# --------------------------------------realize first idea-------------------------------------------
# ---------------------------------------------------------------------------------------------------
# from astropy.io import fits
# import numpy as np
# import os
# galfit_path = '/Users/lpr/Data/galfit/LIRGproject/goodss/'
# gala_list = os.listdir(galfit_path)
# gala = []
# for file in gala_list:
# 	if file.endswith('.fits') or file.endswith('.DS_Store'):
# 		print(file)
# 	else:
# 		gala.append(file)
# catalog_path = '/Users/lpr/Data/fits/expdata/HST/goodss_all/goodss_Huang_van.fits'
# catalog = fits.open(catalog_path)[1].data
# evaluate_catalog = np.zeros((len(gala),5))
# evaluate_catalog_count = 0
# for num1 in range(0,len(gala)):
# 	temp = os.listdir(galfit_path+gala[num1])
# 	if 'residual.fits' in temp:
# 		galaid = int(gala[num1])
# 		re = catalog[np.where(catalog['ID_Huang']==galaid)]['re_f160w']/0.06
# 		residual_image = fits.open(galfit_path+str(galaid)+'/residual.fits')[3].data
# 		flat_image = residual_image.flatten()
# 		background = np.median(flat_image[np.where((flat_image>-0.05)&(flat_image<0.05))])
# 		residual_image = residual_image-background
# 		core_residual = 0
# 		core_count = 0
# 		for x in range(0,residual_image.shape[0]):
# 			x0 = (residual_image.shape[0]-1)/2
# 			for y in range(0,residual_image.shape[1]):
# 				y0 = (residual_image.shape[1]-1)/2
# 				if np.sqrt((x-x0)**2+(y-y0)**2)<re:
# 					core_residual += residual_image[x,y]
# 					core_count += 1
# 		core_residual_sb = core_residual/core_count
# 		evaluate_catalog[num1] = [galaid,re,background,core_residual,core_residual_sb]
# 		evaluate_catalog_count += 1
# 		print(galaid)
# 	else:
# 		print(galaid[num1]+' is not fitted')
# col1 = fits.Column(name='ID_Huang', array=evaluate_catalog[:,0], format='K')
# col2 = fits.Column(name='re_f160w_pixels',array=evaluate_catalog[:,1], format='D')
# col3 = fits.Column(name='background',array=evaluate_catalog[:,2], format='D')
# col4 = fits.Column(name='residual_core',array=evaluate_catalog[:,3], format='D')
# col5 = fits.Column(name='residual_core_sb',array=evaluate_catalog[:,4], format='D')
# table = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5])
# table.writeto(galfit_path+'eval_galfit_result.fits',overwrite=True)
# print('CATALOG done')
# ---------------------------------------------------------------------------------------------------
# --------------------------------------realize second idea------------------------------------------
# ---------------------------------------------------------------------------------------------------
from astropy.io import fits
import numpy as np
import os
galfit_path = '/Users/lpr/Data/galfit/LIRGproject/goodss/'
gala_list = os.listdir(galfit_path)
gala = []
for file in gala_list:
	if file.endswith('.fits') or file.endswith('.DS_Store'):
		print(file)
	else:
		gala.append(file)
catalog_path = '/Users/lpr/Data/fits/expdata/HST/goodss_all/goodss_Huang_van.fits'
catalog = fits.open(catalog_path)[1].data
evaluate_catalog = np.zeros((len(gala),5))
for num1 in range(0,len(gala)):
	temp = os.listdir(galfit_path+gala[num1])
	if 'residual.fits' in temp:
		galaid = int(gala[num1])
		re = catalog[np.where(catalog['ID_Huang']==galaid)]['re_f160w']/0.06
		residual_image = fits.open(galfit_path+str(galaid)+'/residual.fits')[3].data
		original_image = fits.open(galfit_path+str(galaid)+'/residual.fits')[1].data
		deviation_image = abs(residual_image/original_image)
		core_deviation = 0
		core_count = 0
		for x in range(0,deviation_image.shape[0]):
			x0 = (deviation_image.shape[0]-1)/2
			for y in range(0,deviation_image.shape[1]):
				y0 = (deviation_image.shape[1]-1)/2
				if np.sqrt((x-x0)**2+(y-y0)**2)<re:
					core_deviation += deviation_image[x,y]
					core_count += 1
		core_deviation_sb = core_deviation/core_count
		total_deviation = np.sum(deviation_image)
		evaluate_catalog[num1] = [galaid,re,core_deviation_sb,core_count,total_deviation]
		print(galaid)
	else:
		print(str(galaid)+' is not fitted')
col1 = fits.Column(name='ID_Huang', array=evaluate_catalog[:,0], format='K')
col2 = fits.Column(name='re_f160w_pixels',array=evaluate_catalog[:,1], format='D')
col3 = fits.Column(name='core_deviation_sb',array=evaluate_catalog[:,2], format='D')
col4 = fits.Column(name='core_count',array=evaluate_catalog[:,3], format='D')
col5 = fits.Column(name='total_deviation',array=evaluate_catalog[:,4], format='D')
table = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5])
table.writeto(galfit_path+'eval_galfit_result2.fits',overwrite=True)
print('CATALOG done')

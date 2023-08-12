from astropy.io import fits
import os
import numpy as np
import matplotlib.pyplot as plt
path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/spectral/'
cata = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/goodsn_Huang_van_3dhstFITemissionline.fits')[1].data
temp = os.listdir(path)
spectral = []
for file in temp:
	if file.endswith('.fits'):
		spectral.append(file)
for file in spectral:
	hdu = fits.open(path+file)[1].data
	ID_van = file[file.index('_',9)+1:file.rindex('_')]
	number = file[file.rindex('_')+1:file.index('.fits')]
	redshift = cata[np.where(cata['ID_van']==int(ID_van))]['REDSHIFT']
	wave_obs = hdu['wave']
	flux_obs = hdu['flux']
	err_obs = hdu['error']
	wave_res = np.around(wave_obs/(1+redshift),3)
	col1 = fits.Column(name='wave_rest',array=wave_res,format='D')
	col2 = fits.Column(name='flux',array=flux_obs,format='D')
	col3 = fits.Column(name='error',array=err_obs,format='D')
	table = fits.BinTableHDU.from_columns([col1,col2,col3])
	name = path+file[0:-5]
	table.writeto(name+'_rest.fits')
	plt.figure(figsize=[15,7])
	plt.plot(wave_res,flux_obs,linewidth=1,color='blue')
	plt.ylabel('electrons/s')
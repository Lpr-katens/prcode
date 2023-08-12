from astropy.io import fits
import numpy as np
catalog_path = '/Users/lpr/Data/fits/expdata/BTratio/equation/ori_gal_match.fits'
galaxy_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/goodsn_wellFIT_van.fits'
catalog = fits.open(catalog_path)[1].data
galaxy = fits.open(galaxy_path)[1].data
array = np.empty((0,6))
for num1 in range(0,len(galaxy)):
	ID_Huang = int(galaxy[num1]['ID_Huang'])
	re_galaxy = np.around(galaxy[num1]['re_f160w_pix'],2)
	n_galaxy = np.around(galaxy[num1]['n_f160w'],2)
	index_list = np.array([])
	distance = np.array([])
	for num2 in range(0,len(catalog)):
		re_catalog = catalog[num2]['re_total_gal']
		n_catalog = catalog[num2]['n_gal']
		if abs(re_catalog-re_galaxy)<5 and abs(n_catalog-n_galaxy)<1:
			index_list = np.append(index_list,num2)
			distance = np.append(distance,np.sqrt((re_catalog-re_galaxy)**2+(n_catalog-n_galaxy)**2))
	if index_list.size == 0:
		print(str(ID_Huang)+' have no match')
	else:
		index = int(index_list[np.where(distance==np.amin(distance))][0])
		temp = np.array([[ID_Huang,re_galaxy,n_galaxy,np.around(catalog[index]['B2T'],2),np.around(catalog[index]['re_total_gal'],2),np.around(catalog[index]['n_gal'],2)]])
		array = np.append(array,temp,axis=0)
col1 = fits.Column(name='ID_Huang',array=array[:,0],format='K')
col2 = fits.Column(name='re_galaxy',array=array[:,1],format='D')
col3 = fits.Column(name='n_galaxy',array=array[:,2],format='D')
col4 = fits.Column(name='B2T_galfit_match',array=array[:,3],format='D')
col5 = fits.Column(name='re_galfit_match',array=array[:,4],format='D')
col6 = fits.Column(name='n_galfit_match',array=array[:,5],format='D')
hdu = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6])
hdu.writeto('/Users/lpr/Data/fits/expdata/HST/goodsn_all/find_B2T_in_catalog.fits',overwrite=True)
print('done')
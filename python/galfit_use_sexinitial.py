from astropy.io import fits
import numpy as np
import os
root_path = '/Users/lpr/Data/SExtractor/goodsn_sex/'
temp = os.listdir(root_path)
galaxy_list = []
for file in temp:
	if file.endswith('.DS_Store') or file.endswith('.fits'):
		print(file)
	else:
		galaxy_list.append(file)
sextractor_catalog = fits.open(root_path+'decom_result.fits')[1].data
catalog = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/goodsn_Huang_van.fits')[1].data
for file in galaxy_list:
	galaxy_id = int(file)
	galaxy_id_sex = np.where(sextractor_catalog['ID_Huang']==galaxy_id)[0][0]
	galaxy_id_catalog = np.where(catalog['ID_Huang']==galaxy_id)[0][0]
	mag_sph_sex = sextractor_catalog[galaxy_id_sex]['MAG_SPHEROID']
	re_sph_sex = sextractor_catalog[galaxy_id_sex]['SPHEROID_REFF_IMAGE']
	mag_disk_sex = sextractor_catalog[galaxy_id_sex]['MAG_DISK']
	rs_disk_sex = sextractor_catalog[galaxy_id_sex]['DISK_SCALE_IMAGE']
	param_file = open(root_path+file+'/galfit_sex.feedme','w')
	param_file.write(
		' A) gdn_f160w_'+file+'.fits\n'
		' B) residual_sex.fits\n'
		' C) gdn_f160w_wht_'+file+'.fits\n'
		' D) psf_160_stackimage.fits diff_kernel.txt\n'
		' E) 1\n'
		' F) none\n'
		' G) galfit.CONSTRAINTS\n'
		' H) 1    101   1    101\n'
		' I) 40    40\n'
		' J) 25.96\n'
		' K) 0.06  0.06\n'
		' O) regular\n'
		' P) 0\n'
		'# object number: 1\n'
		' 0) devauc\n'
		' 1) 51	51	1	1\n'
		' 3) '+str(np.around(mag_sph_sex,3))+'   1\n'
		' 4) '+str(np.around(re_sph_sex,3))+'   1\n'
		' 9) 1   0\n'
		'10) 0 	 0\n'
		' Z) 0\n'
		'# Object number: 2\n'
		' 0) expdisk\n'
		' 1) 51	51	1	1\n'
		' 3) '+str(np.around(mag_disk_sex,3))+'   1\n'
		' 4) '+str(np.around(rs_disk_sex,3))+'   1\n'
		' 9) '+str(np.around(catalog[galaxy_id_catalog]['q_f160w'],3))+'	1\n'
		'10) '+str(np.around(catalog[galaxy_id_catalog]['pa_f160w'],3))+' 	 1\n'
		' Z) 0\n'
		'# object number: 3\n'
		' 0) sky\n'
		' 1) 0      		1\n'
		' 2) 0.0000      0\n'
		' 3) 0.0000      0\n'
		' Z) 0\n'
		)
	print('------------------ '+file+' all done ------------------')
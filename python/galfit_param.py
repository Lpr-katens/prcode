# ===========================================================
# this code edit galfit configuration parameter file in batch
# and copy configuration, weight images and CCD diffusion 
# kernel to every directory
# ===========================================================
import numpy as np
from astropy.io import fits
from scipy.integrate import quad
import os
from shutil import copy
import subprocess as sp
print('------------------ code begin now ------------------')
expath = '/Users/lpr/Data/galfit/LIRGproject/goodss/'
catalog = fits.open('/Users/lpr/Data/fits/expdata/HST/goodss_all/goodss_Huang_van.fits')[1]
# ori_config_path = '/Users/lpr/Data/galfit/1603/EXAMPLE.CONSTRAINTS'
ori_kernel_path = '/Users/lpr/Data/galfit/1849/diff_kernel.txt'
ori_psf_path = '/Users/lpr/Data/fits/pridata/CANDELS_image/goodss_3dhst.v4.0.F160W_psf.fits'
ori_ima_path = '/Users/lpr/Data/fits/expdata/HST/goodss_all/singleGALAXYimage/goodss_f160w'
ori_wht_path = '/Users/lpr/Data/fits/expdata/HST/goodss_all/singleGALAXYimage_wht/goodss_f160w_wht'
hdu = fits.BinTableHDU.copy(catalog).data
# extract sources with strange re or redshift
hdu = hdu[np.where((hdu['REDSHIFT']<1.3)&(hdu['REDSHIFT']>0.8))]
hdu = hdu[np.where((hdu['f_f160w']==1)|(hdu['f_f160w']==0))]
temp = os.listdir(ori_ima_path)
gala = []
for file in temp:
	if file.endswith('.DS_Store'):
		print(file)
	else:
		gala.append(file)
print('------------------ loop begin now ------------------')
for num in range(0,len(hdu)):
	galaid = hdu[num]['ID_Huang']
	if 'goodss_f160w_'+str(galaid)+'.fits' in gala:
		expath_gala = expath+str(galaid)
		os.system('mkdir '+expath_gala)
		idx = np.where(hdu['ID_Huang']==int(galaid))[0][0]
		file = open(expath_gala+'/galfit.feedme','w')
		file.write(
			' A) goodss_f160w_'+str(galaid)+'.fits\n' # ' A) gds_f160w_'+galaid+'.fits\n'
			' B) residual.fits\n'
			' C) goodss_f160w_wht_'+str(galaid)+'.fits\n' # ' C) gds_f160w_wht_'+galaid+'.fits\n'
			' D) goodss_3dhst.v4.0.F160W_psf.fits  diff_kernel.txt\n'
			' E) 1\n'
			' F) none\n'
			' G) galfit.CONSTRAINTS\n'
			' H) 1    101   1    101\n'
			' I) 60    60\n'
			' J) 25.95\n'
			' K) 0.06  0.06\n'
			' O) regular\n'
			' P) 0\n'
			'# object number: 1\n'
			' 0) devauc\n'
			' 1) 51	51	1	1\n'
			' 3) '+str(np.around(hdu[num]['mag_f160w'],3))+'   1\n'
			' 4) '+str(np.around(hdu[num]['re_f160w']/0.06,3))+'   1\n'
			' 9) 1   0\n'
			'10) 0 	 0\n'
			' Z) 0\n'
			'# Object number: 2\n'
			' 0) expdisk\n'
			' 1) 51	51	1	1\n'
			' 3) '+str(np.around(hdu[num]['mag_f160w'],3))+'   1\n'
			' 4) '+str(np.around(hdu[num]['re_f160w']/0.06+3,3))+'   1\n'
			' 9) '+str(np.around(hdu[num]['q_f160w'],3))+'	1\n'
			'10) '+str(np.around(hdu[num]['pa_f160w'],3))+' 	 1\n'
			' Z) 0\n'
			'# object number: 3\n'
			' 0) sky\n'
			' 1) 0      		1\n'
			' 2) 0.0000      0\n'
			' 3) 0.0000      0\n'
			' Z) 0\n'
			)
		file.close()
		file2 = open(expath_gala+'/galfit.CONSTRAINTS','w')
		file2.write(
			'1    re    1 to 40	 \n'
			'1	 x     -5   5   \n'
			'1	 y     -5   5   \n'
			'1    mag   '+str(np.around(hdu[num]['mag_f160w'],3))+' to 40\n'
			'2    rs    1 to 40	 \n'
			'2	 x     -5   5   \n'
			'2	 y     -5   5   \n'
			'2    mag   '+str(np.around(hdu[num]['mag_f160w'],3))+' to 40\n'
			'2    q		-0.5	0.5\n'
			'2    pa	 -30		30\n'
			)
		file2.close()
		copy(ori_psf_path,expath_gala+'/goodss_3dhst.v4.0.F160W_psf.fits') # copy psf file
		copy(ori_kernel_path,expath_gala+'/diff_kernel.txt') # copy kernel file
		copy(ori_ima_path+'/goodss_f160w_'+str(galaid)+'.fits',expath_gala+'/goodss_f160w_'+str(galaid)+'.fits')
		copy(ori_wht_path+'/goodss_f160w_wht_'+str(galaid)+'.fits',expath_gala+'/goodss_f160w_wht_'+str(galaid)+'.fits') # copy weight images file
		print('------------------ '+str(galaid)+' all done ------------------')
# # import os
# root_path = '/Users/lpr/Data/galfit/LIRGproject/aegis/'
# temp = os.listdir(root_path)
# galaxy_list = []
# for file in temp:
# 	if file.endswith('.DS_Store') or file.endswith('fits'):
# 		print(file)
# 	else:
# 		galaxy_list.append(file)
# for file in galaxy_list:
# 	f = open(root_path+file+'/galfit.01','r')
# 	lines = f.readlines()
# 	lines[8] = 'B) subcomps.fits       # Output data image block\n'
# 	lines[19] = 'P) 3                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n'
# 	f.close()
# 	f = open(root_path+file+'/test.01','w')
# 	f.writelines(lines)
# 	f.close()
# 	print('------------------- '+file+' test.01 is done ------------------')
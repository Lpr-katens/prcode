from astropy.io import fits
import numpy as np
import os
from astropy.modeling.models import Sersic2D
from scipy.special import gammaincinv
from scipy.optimize import least_squares
import subprocess as sp
import time
def sersic_profile(Ire,re,n,r):
	# calculate surface brightness at radius r
	# input:
	# 	Ire: surface brightness in re
	# 	re: effective radius
	# 	n: sersic index
	# 	r: radius
	bn = gammaincinv(2.*n,0.5)
	Intensity_r = Ire*np.exp(-bn*(r/re)**(1/n)-1)
	# r = np.sqrt((x-x0)**2+(y-y0)**2)
	return Intensity_r
root_path = '/Users/lpr/Data/galfit/test/goodsn/'
catalog = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/goodsn_Huang_van.fits')[1].data
temp = os.listdir(root_path)
galaxy_list = []
for file in temp:
	if file.endswith('.DS_Store') or file.endswith('.fits'):
		print(file)
	else:
		galaxy_list.append(file)
for iteration in range(0,10):
	if iteration == 0:
		for file in galaxy_list:
			galaxy_id = int(file)
			idx = np.where(catalog['ID_Huang']==galaxy_id)[0][0]
			os.chdir(root_path)
			os.chdir(file)
			original = fits.open('subcomps.fits')
			f = open('fit.log')
			data = f.readlines()
			for num2 in range(0,len(data)): 
				if 'devauc' in data[num2]:
					x_b = data[num2].split()[3].replace('*',' ')[:-1]
					y_b = data[num2].split()[4].replace('*',' ')[:-1]
					MAG_B = data[num2].split()[5].replace('*',' ')
					re = data[num2].split()[6].replace('*',' ')
				elif 'disk' in data[num2]: 
					x_d = data[num2].split()[3].replace('*',' ')[:-1]
					y_d = data[num2].split()[4].replace('*',' ')[:-1]
					MAG_D = data[num2].split()[5].replace('*',' ')
					rs = data[num2].split()[6].replace('*',' ')
					ba = data[num2].split()[7].replace('*',' ')
					pa = data[num2].split()[8].replace('*',' ')
			feedme = open('devau.feedme','w')
			feedme.write(
				' A) devau_'+str(iteration)+'.fits\n'
				' B) res_devau_'+str(iteration)+'.fits\n'
				' C) gdn_f160w_wht_'+file+'.fits\n'
				' D) psf_160_stackimage.fits diff_kernel.txt\n'
				' E) 1\n'
				' F) none\n'
				' G) devau.CONSTRAINTS\n'
				' H) 1    101   1    101\n'
				' I) 101    101\n'
				' J) 25.96\n'
				' K) 0.06  0.06\n'
				' O) regular\n'
				' P) 0\n'
				'# object number: 1\n'
				' 0) devauc\n'
				' 1) '+x_b+'  '+y_b+'  1 1\n'
				' 3) '+MAG_B+'     1\n'
				' 4) '+re+'   1\n'
				' 9) 1   0\n'
				'10) 0 	 0\n'
				' Z) 0\n'
				)
			feedme.close()
			const = open('devau.CONSTRAINTS','w')
			const.write(
				'1    re    1 to 40	\n'
				'1	 x     -5   5   \n'
				'1	 y     -5   5   \n'
				'1    mag   '+str(catalog[idx]['mag_f160w'])+' to 40\n'
				)
			const.close()
			feedme2 = open('exp.feedme','w')
			feedme2.write(
				' A) exp_'+str(iteration)+'.fits\n'
				' B) res_exp_'+str(iteration)+'.fits\n'
				' C) gdn_f160w_wht_'+file+'.fits\n'
				' D) psf_160_stackimage.fits diff_kernel.txt\n'
				' E) 1\n'
				' F) none\n'
				' G) exp.CONSTRAINTS\n'
				' H) 1    101   1    101\n'
				' I) 101    101\n'
				' J) 25.96\n'
				' K) 0.06  0.06\n'
				' O) regular\n'
				' P) 0\n'
				'# Object number: 1\n'
				' 0) expdisk\n'
				' 1) '+x_d+'    '+y_d+'  1 1\n'
				' 3) '+MAG_D+'     1\n'
				' 4) '+rs+'      1\n'
				' 9) '+ba+'      1\n'
				'10) '+pa+'     1\n'
				' Z) 0\n'
				)
			feedme2.close()
			const2 = open('exp.CONSTRAINTS','w')
			const2.write(
				'1    rs    1 to 40	\n'
				'1	 x     -5   5   \n'
				'1	 y     -5   5   \n'
				'1    mag   '+str(catalog[idx]['mag_f160w'])+' to 40\n'
				'1    q		-0.5	0.5\n'
				'1    pa	-30		30\n'
				)
			const2.close()
		for file in galaxy_list:
			os.chdir(root_path)
			os.chdir(file)
			original = fits.open('subcomps.fits')
			devau = original[0].data-original[2].data
			devau = fits.PrimaryHDU(devau)
			devau.writeto('devau_'+str(iteration)+'.fits',overwrite=True)
			print('devau_'+str(iteration)+'.fits has been written')
			print(file) 
			sp.run('galfit devau.feedme',shell=True,check=True)
			print('devau '+str(iteration))
			exp = original[0].data-fits.open('res_devau_'+str(iteration)+'.fits')[2].data
			exp = fits.PrimaryHDU(exp)
			exp.writeto('exp_'+str(iteration)+'.fits',overwrite=True)
			print('exp_'+str(iteration)+'.fits has been written')
		for file in galaxy_list:
			os.chdir(root_path)
			os.chdir(file)
			print(file)
			sp.run('galfit exp.feedme',shell=True,check=True)
			print('exp '+str(iteration))
			print('------------'+str(iteration)+' '+file+' is done------------')
	else:
		for file in galaxy_list:
			galaxy_id = int(file)
			idx = np.where(catalog['ID_Huang']==galaxy_id)[0][0]
			os.chdir(root_path)
			os.chdir(file)
			original = fits.open('subcomps.fits')
			f = open('fit.log')
			data = f.readlines()
			for num2 in range(0,len(data)): 
				if 'devauc' in data[num2]:
					x_b = data[num2].split()[3].replace('*',' ')[:-1]
					y_b = data[num2].split()[4].replace('*',' ')[:-1]
					MAG_B = data[num2].split()[5].replace('*',' ')
					re = data[num2].split()[6].replace('*',' ')
				elif 'disk' in data[num2]: 
					x_d = data[num2].split()[3].replace('*',' ')[:-1]
					y_d = data[num2].split()[4].replace('*',' ')[:-1]
					MAG_D = data[num2].split()[5].replace('*',' ')
					rs = data[num2].split()[6].replace('*',' ')
					ba = data[num2].split()[7].replace('*',' ')
					pa = data[num2].split()[8].replace('*',' ')
			feedme = open('devau.feedme','w')
			feedme.write(
				' A) devau_'+str(iteration)+'.fits\n'
				' B) res_devau_'+str(iteration)+'.fits\n'
				' C) gdn_f160w_wht_'+file+'.fits\n'
				' D) psf_160_stackimage.fits diff_kernel.txt\n'
				' E) 1\n'
				' F) none\n'
				' G) devau.CONSTRAINTS\n'
				' H) 1    101   1    101\n'
				' I) 101    101\n'
				' J) 25.96\n'
				' K) 0.06  0.06\n'
				' O) regular\n'
				' P) 0\n'
				'# object number: 1\n'
				' 0) devauc\n'
				' 1) '+x_b+'  '+y_b+'  1 1\n'
				' 3) '+MAG_B+'     1\n'
				' 4) '+re+'   1\n'
				' 9) 1   0\n'
				'10) 0 	 0\n'
				' Z) 0\n'
				)
			feedme.close()
			feedme2 = open('exp.feedme','w')
			feedme2.write(
				' A) exp_'+str(iteration)+'.fits\n'
				' B) res_exp_'+str(iteration)+'.fits\n'
				' C) gdn_f160w_wht_'+file+'.fits\n'
				' D) psf_160_stackimage.fits diff_kernel.txt\n'
				' E) 1\n'
				' F) none\n'
				' G) exp.CONSTRAINTS\n'
				' H) 1    101   1    101\n'
				' I) 101    101\n'
				' J) 25.96\n'
				' K) 0.06  0.06\n'
				' O) regular\n'
				' P) 0\n'
				'# Object number: 1\n'
				' 0) expdisk\n'
				' 1) '+x_d+'    '+y_d+'  1 1\n'
				' 3) '+MAG_D+'     1\n'
				' 4) '+rs+'      1\n'
				' 9) '+ba+'      1\n'
				'10) '+pa+'     1\n'
				' Z) 0\n'
				)
			feedme2.close()
		for file in galaxy_list:
			os.chdir(root_path)
			os.chdir(file)
			original = fits.open('subcomps.fits')
			devau = original[0].data-fits.open('res_exp_'+str(iteration-1)+'.fits')[2].data
			devau = fits.PrimaryHDU(devau)
			devau.writeto('devau_'+str(iteration)+'.fits',overwrite=True)
			print('devau_'+str(iteration)+'.fits has been written')
			print(file)
			sp.run('galfit devau.feedme',shell=True,check=True)
			print('devau '+str(iteration))
			exp = original[0].data-fits.open('res_devau_'+str(iteration)+'.fits')[2].data
			exp = fits.PrimaryHDU(exp)
			exp.writeto('exp_'+str(iteration)+'.fits',overwrite=True)
			print('exp_'+str(iteration)+'.fits has been written')
		for file2 in galaxy_list:
			os.chdir(root_path)
			os.chdir(file2)
			print(file2)
			sp.run('galfit exp.feedme',shell=True,check=True)
			print('exp '+str(iteration))
			print('------------'+str(iteration)+' '+file2+' is done------------')
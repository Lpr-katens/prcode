import numpy as np
from astropy.io import fits
from scipy.integrate import quad
import os
from shutil import copy
import subprocess as sp

rootpath = '/Users/lpr/Data/galfit/LIRGproject/BT_catalog/'
# ddd = os.listdir(expath)
# galaxy = []
# for files in ddd:
# 	if files.endswith('.DS_Store') or files.endswith('.fits'):
# 		print(files)
# 	elif files.startswith('0.9') or files.startswith('0.95'):
# 		temp = os.listdir(expath+'/'+files)
# 		if 'fit.log' not in temp:
# 			galaxy.append(files)
# for num in range(0,len(galaxy)):
# 	print('     *              *\n'
# 		'     ***            ***\n'
# 		'    *****          *****\n'
# 		'   **********************\n'
# 		' **************************\n'
# 		'*********  galfit  *********\n'
# 		'*********   begin  *********\n'
# 		' ********   now    ********\n'
# 		'   **********************\n'
# 		'     ******************\n'
# 		'       **************\n'
# 		'         **********\n'
# 		)
# 	os.chdir(expath)
# 	os.chdir(galaxy[num])
# 	sp.run('galfit galfit.feedme',shell=True,check=True)
# 	print('------------------ '+str(galaxy[num])+' is done ------------------')

re_devau_size = [51,2]
re_exp_size = [61,2]
for BT_ratio in np.arange(0.05,1.01,0.05):
	for re_devau in np.arange(1,5,re_devau_size[1]):
		for re_exp in np.arange(re_devau,re_exp_size[0],re_exp_size[1]):
			expath = rootpath+str(np.around(BT_ratio,2))+'-'+str(np.around(re_devau,2))+'-'+str(np.around(re_exp,2))
			temp = os.listdir(expath)
			if 'fit.log' not in temp:
				print('     *              *\n'
					'     ***            ***\n'
					'    *****          *****\n'
					'   **********************\n'
					' **************************\n'
					'*********  galfit  *********\n'
					'*********   begin  *********\n'
					' ********   now    ********\n'
					'   **********************\n'
					'     ******************\n'
					'       **************\n'
					'         **********\n'
					)
				os.chdir(expath)
				sp.run('galfit galfit.feedme',shell=True,check=True)
				print('------------------ '+str(np.around(BT_ratio,2))+'-'+str(np.around(re_devau,2))+'-'+str(np.around(re_exp,2))+' is done ------------------')
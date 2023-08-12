import numpy as np
from astropy.io import fits
from scipy.integrate import quad
import os
from shutil import copy
import subprocess as sp


expath = '/Users/lpr/Data/galfit/LIRGproject/goodss/'# expath = '/Users/lpr/Data/SExtractor/goodsn_sex'
ddd = os.listdir(expath)
galaxy = []
for files in ddd:
	if files.endswith('.DS_Store') or files.endswith('.fits'):
		print(files)
	else:
		galaxy.append(files)
for num in range(0,len(galaxy)):
	temp = os.listdir(expath+'/'+galaxy[num])
	if 'residual.fits' not in temp:
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
		os.chdir(galaxy[num])
		sp.run('galfit galfit.feedme',shell=True,check=True)
		print('------------------ '+str(galaxy[num])+' is done ------------------')
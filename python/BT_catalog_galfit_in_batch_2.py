import numpy as np
from astropy.io import fits
from scipy.integrate import quad
import os
from shutil import copy
import subprocess as sp

expath = '/Users/lpr/Data/galfit/LIRGproject/BT_catalog/'
ddd = os.listdir(expath)
galaxy = []
for files in ddd:
	if files.endswith('.DS_Store') or files.endswith('.fits'):
		print(files)
	elif files.startswith('0.8') or files.startswith('0.85'):
		temp = os.listdir(expath+'/'+files)
		if 'fit.log' not in temp:
			galaxy.append(files)
for num in range(0,len(galaxy)):
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
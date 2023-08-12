import numpy as np
from astropy.io import fits
from scipy.integrate import quad
import os
from shutil import move
import subprocess as sp

expath = '/Users/lpr/Data/galfit/LIRGproject/BT_catalog/'
ori_ima_path = '/Users/lpr/Data/fits/expdata/BTratio/equation/'
image_list = []
for file in os.listdir(ori_ima_path):
	if file.endswith('.fits'):
		image_list.append(file)
# for file in os.listdir(expath):
# 	if file.endswith('.fits') or file.endswith('DS_Store'):
# 		print(file)
# 	else:
# 		image_list.append(file)
for num1 in range(0,len(image_list)):
	galaname = image_list[num1][0:image_list[num1].index('.fits')] # [0:image_list[num1].index('.fits')]
	bt_ratio = image_list[num1][0:image_list[num1].index('-')]
	re_devau = image_list[num1][image_list[num1].index('-')+1:image_list[num1].index('-',5)]
	re_exp = image_list[num1][image_list[num1].index('-',5)+1:image_list[num1].index('.fits')] # image_list[num1].index('.fits')
	print(bt_ratio,re_devau,re_exp)
	expath_gala = expath+galaname
	os.system('mkdir '+expath_gala)
	file = open(expath_gala+'/galfit.feedme','w')
	file.write(
		' A) '+galaname+'.fits\n'
		' B) residual.fits\n'
		' C) none\n'
		' D) none\n'
		' E) 1\n'
		' F) none\n'
		' G) galfit.CONSTRAINTS\n'
		' H) 1	2001	1	2001\n'
		' I) 100	100\n'
		' J) 25\n'
		' K) 0.06  0.06\n'
		' O) regular\n'
		' P) 0\n'
		'# object number: 1\n'
		' 0) sersic\n'
		' 1) 1001	1001	1	1\n'
		' 3) 20   1\n'
		' 4) '+str(np.around(float(re_devau)+(1-float(bt_ratio))*(float(re_exp)-float(re_devau)),2))+'   1\n'
		' 5) '+str(np.around(float(bt_ratio)*3+1,2))+'		1\n'
		' 9) 1	0\n'
		'10) 0 	 0\n'
		' Z) 0\n'
		)
	file2 = open(expath_gala+'/galfit.CONSTRAINTS','w')
	file2.write(
		'1    re    '+re_devau+' to '+re_exp+'	\n'
		'1	 x     -5   5   \n'
		'1	 y     -5   5   \n'
		'1	 n     0.2	to   8   \n'
		)
	move(ori_ima_path+image_list[num1],expath_gala+'/'+image_list[num1])
	print('------------------ '+image_list[num1]+' is done ------------------')
from astropy.nddata import CCDData
import numpy as np
from petrofit.segmentation import make_catalog
from petrofit.photometry import order_cat
import os
import matplotlib.pyplot as plt
from astropy.io import fits
import statmorph
from astropy.wcs import WCS
import photutils
from lpr.image import morphology as mp
from astropy.convolution import Gaussian2DKernel, convolve
from photutils.aperture import EllipticalAnnulus as ea
from photutils.isophote import Ellipse
from photutils.morphology import data_properties
import time
from scipy import ndimage as ndi
from shutil import copy
import subprocess as sp
import astropy.units as u

for eta in [0.2,0.25]:#np.arange(0.05,0.5,0.05)
	txt = open('/Users/lpr/Data/lirg_project/output/catalog_radec/goodsn_gini_m20_115petro_'+str(np.around(eta,2))+'eta_test.txt','w')
	fields = ['goodsn']#,'goodss','egs'
	psf_path = {'goodsn':'goodsn','goodss':'goodss','egs':'aegis'}
	image_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
	ctg_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
	segm_path = '/Users/lpr/Data/lirg_project/output/gini_m20_segmap/'
	band_used = 'f160w'
	even_check = {False:0,True:1}
	aper = 17 #设置计算gini和m20的孔径为1个角秒，大概是17个像素,0.5"=8.3pixels
	# 满足E/S0/Sa条件的星系
	regular_list = {'goodsn':[15112,13404,15995,15722,16115,13918,10534,24051,9988,11762,10685,13054,5383,14258,9457,22595,6981,5788,5271,2929,5075,20279,29082,10738,16441,14130,17749,11003,6248,24758,26797,13814,5297,29349,14780,30540,12762,16957,18654,27912,6513,3668,5858,1323,17308,1209,26506,1111,28294,32419,32851,27324,32014,34293,1114,12538,33568,12288,22758,3631,35739,29451,34043,23486,30995,13012,26118,13130,10645,26024,2142,28448,14887,25856,36636,32849,12336,31462,32411,9720,24141,14480,10776,23076,15742,14500,34371,29304,31417,38001,15758,25587,23844,19436,6796,25645,35985,27302,14565,23316,24082,28688,30235,23971,24705,30322,26926,33021,18711,32963,33322,15477,29835,24080,23226,16767,17565,33272,34501,24035,30265,13906],'goodss':[1928,2767,3452,4490,3919,26981,3631,6721,6825,8519,8326,8049,7223,10680,11793,11657,11957,12040,12993,12780,12353,13193,12997,13326,13342,14707],'egs':[13010519,13101979,13017707,13033913,13026220,13018656,13003614,13026194,13042381,12024267,13010765,12101109,13004254,13034520,13034076,13011067,13042029,13004240,13026626,12028644,12028644,12028094,12028092,13011748,13026391,13034272,12028563,13034019,13026825,13019229,13026703,13026797,12028490,13035131,13026418,13034617,13011723]}
	for field in fields:
		psf = fits.getdata('/Users/lpr/Data/lirg_project/intake/CANDELS/'+psf_path[field]+'_3dhst_v4.0_wfc3_psf/'+psf_path[field]+'_3dhst.v4.0.F160W_psf.fits',0)
		path_1 = image_path + field + '_10/' + field + '_' + band_used
		id_ctg = fits.getdata(ctg_path+field+'_id.fits',1)
		idx_list = np.full(len(id_ctg),-99)
		gini_m20_statmorph = np.full((len(id_ctg),2),-999.)
		gini_m20_lprnoaper = np.full((len(id_ctg),2),-999.) #自己写的程序，不用aper
		gini_m20_lpraper = np.full((len(id_ctg),2),-999.) #自己写的程序，设置aper
		h_ctg = fits.getdata(ctg_path+field+'_Huangall_candels_radec_van_zmasssfr8sfrhuangsfr.fits',1)
		for num1 in range(0,len(id_ctg)):#
			idx = id_ctg[num1]['id']
			st = time.process_time()
			if idx in regular_list[field]:#idx==12773:#
				print('ID:'+str(idx))
				print('------------------------------------------')
				# print('         *              *\n'
				#   	'        %$@            !@#\n'
				#   	'       &%*&%          ^$#!@\n'
				#   	'      <:{!+@"|&^*&%$#!@%"}#)\n'
				#   	'     #!@%$**************&%$@~\n'
				#   	'    &++&@%***  code  ***%&%*&%\n'
				#   	'    @%*#!&**   begin  **"?)^&%\n'
				#   	'     %*#!&**   now    **?%*@#\n'
				#   	'      ++&@**************@~#^\n'
				#   	'          #!*&&^$#!@%*()\n'
				#   	'            $%(_+(%@?^\n'
				#   	'         )!(#*&++&@%%!*_#+\n'
				#   	'      !@#$)(*&^$#!@%*()*&%$@\n'
				#   	'    (*&^$#!@%*()*&%$@?#^%*(*&*\n'
				#   	'  #$)(*&^$#<{}"?)^&%$@?#^%*(*&*)\n'
				#   	' !@#$)(*&^$#!@%*&%*&%$@~#^%*(*&*)\n'
				#   	'  )(*&^$#!@%*$#*&^$#!@%*()*&%$@~\n'
				#   	'   $)(*&^$#!@%*!@*&%$@|#^%&%$@~\n'
				#   	'    (*&^$#!@%*@|_>"}#)$(*()*&%          **   **\n'
				#   	'     (*&^$#!@%*#!*&&^$#!@%*()       **          **\n'
				#   	'      )*&%$@|%*@%*&&^$#!@%*(      **      **     **\n'
				#   	'       @|#^%&%$>:}?%*@#!@%*     **       **      **\n'
				#   	'        @%*@#!@%*!@%*$#*&@*     **         **   **\n'
				#   	'        !@%*()*&:>%#!~>:"& ** **\n'
				#   	'        <:{!+@"|&^*&%$#!@%   *\n'
				#   	'       *$#!+@"*&  &%$@~<{}"')
				ra,dec = h_ctg[h_ctg['id']==idx]['ra_candels'],h_ctg[h_ctg['id']==idx]['dec_candels']
				image = fits.getdata(path_1 + '/' + field + '_f160w_' + str(idx) + '.fits',0)
				hdr = fits.getheader(path_1 + '/' + field + '_f160w_' + str(idx) + '.fits',0)
				hdr['CTYPE1'] = 'RA---TAN-SIP'
				hdr['CTYPE2'] = 'DEC--TAN-SIP'
				x,y = WCS(hdr).wcs_world2pix(ra,dec,0)
				print('coordinate:'+str([x,y]))
				idx_list[num1] = idx
				# |先计算背景|
				noise = np.median(image[(image>-0.01)&(image<0.01)].flatten())
				noise_sigma = np.std(image[(image>-0.01)&(image<0.01)].flatten())
				noise_map = np.random.normal(noise,noise_sigma,image.shape)
				image = image - noise
				npixels = 21
				try:
					threshold = photutils.detect_threshold(image,3)
					segm = photutils.detect_sources(image,threshold,npixels)
					label = np.argmax(segm.areas) + 1
					segmap = segm.data == label
					segmap = np.where(segmap==True,1,0)
					source_morphs = statmorph.source_morphology(image,segmap,gain=1,psf=None,eta=eta)
					morph = source_morphs[segmap[int(x),int(y)]-1]
					gini_stat = morph.gini
					m20_stat = morph.m20
					gini_m20_statmorph[num1] = [m20_stat,gini_stat]
					print(gini_m20_statmorph[num1])
				except:
					continue
				et = time.process_time()
				print(field+' '+str(idx)+'time consume: ' + str(et - st))
		col0 = fits.Column(name='id', array = idx_list, format = 'K')
		col1 = fits.Column(name='Moment_20_statmorph', array=gini_m20_statmorph[:,0], format='D')
		col2 = fits.Column(name='Gini_coeff_statmorph', array = gini_m20_statmorph[:,1], format='D')
		hdu = fits.BinTableHDU.from_columns([col0,col1,col2])
		print(image_path + field + '_GiniM20_regulars_115petro_3detectsigma_'+str(np.around(eta,2))+'eta.fits')
		hdu.writeto(image_path + field + '_GiniM20_regulars_115petro_3detectsigma_'+str(np.around(eta,2))+'eta.fits',overwrite=True)
	
	# 开始画图
	ctg_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
	img_path = '/Users/lpr/Data/lirg_project/output/images/'
	merger_path = '/Users/lpr/Data/lirg_project/output/gini_m20_segmap/'
	fields = ['goodsn']#,'goodss','egs'
	pair_flag = {True:'red',False:'blue'}
	x = np.linspace(-3,0,100)
	y1 = -0.14*x+0.33
	y2 = 0.14*x+0.8
	count_specpair = []
	regular_list = {'goodsn':[15112,13404,15995,15722,16115,13918,10534,24051,9988,11762,10685,13054,5383,14258,9457,22595,6981,5788,5271,2929,5075,20279,29082,10738,16441,14130,17749,11003,6248,24758,26797,13814,5297,29349,14780,30540,12762,16957,18654,27912,6513,3668,5858,1323,17308,1209,26506,1111,28294,32419,32851,27324,32014,34293,1114,12538,33568,12288,22758,3631,35739,29451,34043,23486,30995,13012,26118,13130,10645,26024,2142,28448,14887,25856,36636,32849,12336,31462,32411,9720,24141,14480,10776,23076,15742,14500,34371,29304,31417,38001,15758,25587,23844,19436,6796,25645,35985,27302,14565,23316,24082,28688,30235,23971,24705,30322,26926,33021,18711,32963,33322,15477,29835,24080,23226,16767,17565,33272,34501,24035,30265,13906],'goodss':[1928,2767,3452,4490,3919,26981,3631,6721,6825,8519,8326,8049,7223,10680,11793,11657,11957,12040,12993,12780,12353,13193,12997,13326,13342,14707],'egs':[13010519,13101979,13017707,13033913,13026220,13018656,13003614,13026194,13042381,12024267,13010765,12101109,13004254,13034520,13034076,13011067,13042029,13004240,13026626,12028644,12028644,12028094,12028092,13011748,13026391,13034272,12028563,13034019,13026825,13019229,13026703,13026797,12028490,13035131,13026418,13034617,13011723]}
	for field in fields:
		ginim20_ctg = fits.getdata(ctg_path+field+'_GiniM20_regulars_115petro_3detectsigma_'+str(np.around(eta,2))+'eta.fits',1)
		pair_ctg = fits.getdata(ctg_path+field+'_pair_1000kms.fits',1)
		plt.figure(figsize=[5,5])
		plt.xlabel('$M_{20}$',fontsize=15)
		plt.ylabel('Gini',fontsize=15)
		plt.plot(x,y1,color='black')
		plt.plot(x,y2,color='black')
		for num in range(0,len(ginim20_ctg)):
			idx = ginim20_ctg[num]['id']
			if idx != -99:
				gini = ginim20_ctg[num]['gini_coeff_statmorph']
				m20 = ginim20_ctg[num]['moment_20_statmorph']
				count_spec_pair = pair_ctg[pair_ctg['id']==idx]['count_spec'][0]
				plt.scatter(m20,gini,color=pair_flag[count_spec_pair>0],s=5)
				if gini<-0.14*m20+0.33 and gini>0.14*m20+0.8:#pair_ctg[pair_ctg['id']==idx]['count_spec'][0]>0 and 
					count_specpair.append(idx)
		plt.xlim(0,-3)
		plt.ylim(0.3,0.8)
		plt.savefig(img_path+field+'_gini_m20_regulars_115petro_3detectsigma_'+str(np.around(eta,2))+'eta.png')
	print(len(count_specpair)/len(regular_list['goodsn']))
	txt.write('E/S0/Sa:'+str(len(count_specpair))+' All:'+str(len(regular_list['goodsn'])))
	txt.write('ratio:'+str(len(count_specpair)/len(regular_list['goodsn'])))
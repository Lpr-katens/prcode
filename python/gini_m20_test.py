# from astropy.nddata import CCDData
# import numpy as np
# from petrofit.segmentation import make_catalog
# from petrofit.photometry import order_cat
# import os
# import matplotlib.pyplot as plt
# from astropy.io import fits
# import statmorph
# from astropy.wcs import WCS
# import photutils
# from lpr.image import morphology as mp
# from astropy.convolution import Gaussian2DKernel, convolve
# from photutils.aperture import EllipticalAnnulus as ea
# from photutils.isophote import Ellipse
# import time
# from scipy import ndimage as ndi
# from shutil import copy
# import subprocess as sp

# fields = ['goodsn','goodss','egs']#
# psf_path = {'goodsn':'goodsn','goodss':'goodss','egs':'aegis'}
# image_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
# ctg_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
# segm_path = '/Users/lpr/Data/lirg_project/output/gini_m20_segmap/'
# band_used = 'f160w'
# even_check = {False:0,True:1}
# # detection_threshold = 3
# # mtotal_box = 5
# # npixels = 5
# aper = 17 #设置计算gini和m20的孔径为1个角秒，大概是17个像素,0.5"=8.3pixels
# for field in fields:
# 	psf = fits.open('/Users/lpr/Data/lirg_project/intake/CANDELS/'+psf_path[field]+'_3dhst_v4.0_wfc3_psf/'+psf_path[field]+'_3dhst.v4.0.F160W_psf.fits')[0].data
# 	path_1 = image_path + field + '/' + field + '_' + band_used#_10
# 	id_ctg = fits.getdata(ctg_path+field+'_id.fits',1)
# 	idx_list = np.full(len(id_ctg),-99)
# 	gini_m20_statmorph = np.full((len(id_ctg),2),-999.)
# 	gini_m20_statmorph_segm = np.full((len(id_ctg),2),-999.)
# 	gini_m20_lprnoaper = np.full((len(id_ctg),2),-999.) #自己写的程序，不用aper
# 	gini_m20_lpraper = np.full((len(id_ctg),2),-999.) #自己写的程序，设置aper
# 	h_ctg = fits.open(ctg_path+field+'_Huangall_candels_radec_van_zmasssfr8sfrhuangsfr.fits')[1].data
# 	for num1 in range(0,len(id_ctg)):#
# 		idx = id_ctg[num1]['id']
# 		st = time.process_time()
# 		if idx!= 32411:#==12773:
# 			print('ID:'+str(idx))
# 			print('         *              *\n'
# 				  '        %$@            !@#\n'
# 				  '       &%*&%          ^$#!@\n'
# 				  '      <:{!+@"|&^*&%$#!@%"}#)\n'
# 				  '     #!@%$**************&%$@~\n'
# 				  '    &++&@%***  code  ***%&%*&%\n'
# 				  '    @%*#!&**   begin  **"?)^&%\n'
# 				  '     %*#!&**   now    **?%*@#\n'
# 				  '      ++&@**************@~#^\n'
# 				  '          #!*&&^$#!@%*()\n'
# 				  '            $%(_+(%@?^\n'
# 				  '         )!(#*&++&@%%!*_#+\n'
# 				  '      !@#$)(*&^$#!@%*()*&%$@\n'
# 				  '    (*&^$#!@%*()*&%$@?#^%*(*&*\n'
# 				  '  #$)(*&^$#<{}"?)^&%$@?#^%*(*&*)\n'
# 				  ' !@#$)(*&^$#!@%*&%*&%$@~#^%*(*&*)\n'
# 				  '  )(*&^$#!@%*$#*&^$#!@%*()*&%$@~\n'
# 				  '   $)(*&^$#!@%*!@*&%$@|#^%&%$@~\n'
# 				  '    (*&^$#!@%*@|_>"}#)$(*()*&%          **   **\n'
# 				  '     (*&^$#!@%*#!*&&^$#!@%*()       **          **\n'
# 				  '      )*&%$@|%*@%*&&^$#!@%*(      **      **     **\n'
# 				  '       @|#^%&%$>:}?%*@#!@%*     **       **      **\n'
# 				  '        @%*@#!@%*!@%*$#*&@*     **         **   **\n'
# 				  '        !@%*()*&:>%#!~>:"& ** **\n'
# 				  '        <:{!+@"|&^*&%$#!@%   *\n'
# 				  '       *$#!+@"*&  &%$@~<{}"')
# 			ra,dec = h_ctg[h_ctg['id']==idx]['ra_candels'],h_ctg[h_ctg['id']==idx]['dec_candels']
# 			image = fits.getdata(path_1 + '/' + field + '_f160w_' + str(idx) + '.fits',0)
# 			hdr = fits.getheader(path_1 + '/' + field + '_f160w_' + str(idx) + '.fits',0)
# 			hdr['CTYPE1'] = 'RA---TAN-SIP'
# 			hdr['CTYPE2'] = 'DEC--TAN-SIP'
# 			x,y = WCS(hdr).wcs_world2pix(ra,dec,0)
# 			print('coordinate:'+str([x,y]))
# 			idx_list[num1] = idx
# 			# |先计算背景|
# 			noise = np.median(image[(image>-0.01)&(image<0.01)].flatten())
# 			noise_sigma = np.std(image[(image>-0.01)&(image<0.01)].flatten())
# 			noise_map = np.random.normal(noise,noise_sigma,image.shape)
# 			# |把计算半光半径的程序生成的segmentation拿来直接用|
# 			segmap = fits.getdata('/Users/lpr/Data/lirg_project/output/van_re_half_light_radius/segmentation/'+str(idx)+'.fits',0)
# 			# segm = fits.getdata('/Users/lpr/Data/lirg_project/output/gini_m20_segmap/'+field+'_'+str(idx)+'.fits',0)
# 			# |扣背景和污染源|
# 			# image[(segm != segm[int(x),int(y)])&(segm != 0)] = np.nan
# 			image = image - noise
# 			# |用statmorph直接出结果|
# 			try:
# 				source_morphs = statmorph.source_morphology(image,segmap,gain=1,psf=None)
# 				morph = source_morphs[segmap[int(x),int(y)]-1]
# 				gini_stat = morph.gini
# 				m20_stat = morph.m20
# 				gini_m20_statmorph[num1] = [m20_stat,gini_stat]
# 				print(gini_m20_statmorph[num1])
# 			except:
# 				continue
# 			et = time.process_time()
# 			print(field+' '+str(idx)+'time consume: ' + str(et - st))
# 	col0 = fits.Column(name='id', array = idx_list, format = 'K')
# 	col1 = fits.Column(name='Moment_20_statmorph', array=gini_m20_statmorph[:,0], format='D')
# 	col2 = fits.Column(name='Gini_coeff_statmorph', array = gini_m20_statmorph[:,1], format='D')
# 	hdu = fits.BinTableHDU.from_columns([col0,col1,col2])#
# 	print(image_path + field + '_GiniM20_test.fits')
# 	hdu.writeto(image_path + field + '_GiniM20_test.fits',overwrite=True)
# -------------------------------------------------------------------------------------------------
# ================================================================================================
# -------------------------------------------------------------------------------------------------
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
			print('         *              *\n'
			  	'        %$@            !@#\n'
			  	'       &%*&%          ^$#!@\n'
			  	'      <:{!+@"|&^*&%$#!@%"}#)\n'
			  	'     #!@%$**************&%$@~\n'
			  	'    &++&@%***  code  ***%&%*&%\n'
			  	'    @%*#!&**   begin  **"?)^&%\n'
			  	'     %*#!&**   now    **?%*@#\n'
			  	'      ++&@**************@~#^\n'
			  	'          #!*&&^$#!@%*()\n'
			  	'            $%(_+(%@?^\n'
			  	'         )!(#*&++&@%%!*_#+\n'
			  	'      !@#$)(*&^$#!@%*()*&%$@\n'
			  	'    (*&^$#!@%*()*&%$@?#^%*(*&*\n'
			  	'  #$)(*&^$#<{}"?)^&%$@?#^%*(*&*)\n'
			  	' !@#$)(*&^$#!@%*&%*&%$@~#^%*(*&*)\n'
			  	'  )(*&^$#!@%*$#*&^$#!@%*()*&%$@~\n'
			  	'   $)(*&^$#!@%*!@*&%$@|#^%&%$@~\n'
			  	'    (*&^$#!@%*@|_>"}#)$(*()*&%          **   **\n'
			  	'     (*&^$#!@%*#!*&&^$#!@%*()       **          **\n'
			  	'      )*&%$@|%*@%*&&^$#!@%*(      **      **     **\n'
			  	'       @|#^%&%$>:}?%*@#!@%*     **       **      **\n'
			  	'        @%*@#!@%*!@%*$#*&@*     **         **   **\n'
			  	'        !@%*()*&:>%#!~>:"& ** **\n'
			  	'        <:{!+@"|&^*&%$#!@%   *\n'
			  	'       *$#!+@"*&  &%$@~<{}"')
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
			# # |把计算半光半径的程序生成的segmentation拿来直接用|
			# segmap = fits.getdata('/Users/lpr/Data/lirg_project/output/van_re_half_light_radius/segmentation/'+str(idx)+'.fits',0)
			# segmap = fits.getdata('/Users/lpr/Data/lirg_project/output/gini_m20_segmap/'+field+'_'+str(idx)+'.fits',0)
			# |扣背景和污染源|
			# image[(segmap != segmap[int(x),int(y)])&(segmap != 0)] = 0.
			image = image - noise
			# # |整理segmap，是目标源的地方赋值label=1，其余地方赋值为0|
			# segmap_copy = np.copy(segmap)
			# segmap = np.where(segmap_copy==segmap[int(x),int(y)],1,0)
			# |用statmorph直接出结果|
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
			# # |先用photutils来拟合星系的椭圆|
			# cat = data_properties(image)
			# position = (cat.xcentroid, cat.ycentroid)
			# r = 3.0  # approximate isophotal extent
			# a = cat.semimajor_sigma.value*r
			# b = cat.semiminor_sigma.value*r
			# theta = cat.orientation.to(u.rad).value
			# # |先做petrosianr, 然后smooth 0.2*petrosian半径
			# center = mp.m20_center(image,segmap,segmap[int(x),int(y)])
			# petror = mp.petror_ellip(image,center,amax=a,bmax=b,theta=theta)
			# image_smooth = ndi.gaussian_filter(image,0.2*petror)
			# # kernel = Gaussian2DKernel(x_stddev=0.2*petror,y_stddev=0.2*petror,x_size=image.shape[0]-even_check[(image.shape[0] % 2) == 0],y_size=image.shape[1]-even_check[(image.shape[1] % 2) == 0])
			# # image_copy = convolve(image,kernel)
			# # |再拟合平滑后的星系的椭圆|
			# cat = data_properties(image_smooth)
			# position = (cat.xcentroid, cat.ycentroid)
			# r = 3.0  # approximate isophotal extent
			# a = cat.semimajor_sigma.value*r
			# b = cat.semiminor_sigma.value*r
			# theta = cat.orientation.to(u.rad).value
			# center = mp.m20_center(image_smooth,segmap,segmap[int(x),int(y)])
			# petror = mp.petror_ellip(image_smooth,center,amax=a,bmax=b,theta=theta)
			# temp_img,segmap = mp.galaxy_pixels(image_smooth,petror,center,axis_ratio=b/a,theta=theta)
			# # try:
			# # |用自己写的程序，不设置aper|
			# # gini_lprnoaper = mp.gini(image,segmap,segmap[int((image.shape[0]-1)/2),int((image.shape[1]-1)/2)])
			# # m20_lprnoaper = mp.mtwenty(mp.mtotal(image,x,y,mtotal_box,segmap,segmap[int((image.shape[0]-1)/2),int((image.shape[1]-1)/2)])[2],image,x,y,segmap,segmap[int((image.shape[0]-1)/2),int((image.shape[1]-1)/2)])
			# gini_lprnoaper = mp.gini(image,segmap,petror,center,axis_ratio=b/a,theta=theta,label=segmap[int(x),int(y)])
			# m20_lprnoaper = mp.mtwenty(image,segmap,petror,center,axis_ratio=b/a,theta=theta,label=segmap[int(x),int(y)])
			# gini_m20_lprnoaper[num1] = [m20_lprnoaper,gini_lprnoaper]
			# # |用自己写的程序，设置aper=17个像素|
			# # gini_lpraper = mp.gini(image,segmap,segmap[int((image.shape[0]-1)/2),int((image.shape[1]-1)/2)],aper=aper,x_centroid=x,y_centroid=y)
			# # m20_lpraper = mp.mtwenty(mp.mtotal(image,x,y,mtotal_box,segmap,segmap[int((image.shape[0]-1)/2),int((image.shape[1]-1)/2)],aper=aper,x_centroid=x,y_centroid=y)[2],image,x,y,segmap,segmap[int((image.shape[0]-1)/2),int((image.shape[1]-1)/2)],aper=aper)
			# gini_lpraper = mp.gini(image,segmap,petror,center,axis_ratio=b/a,theta=theta,aper=aper,label=segmap[int(x),int(y)])
			# m20_lpraper = mp.mtwenty(image,segmap,petror,center,axis_ratio=b/a,theta=theta,aper=aper,label=segmap[int(x),int(y)])
			# gini_m20_lpraper[num1] = [m20_lpraper,gini_lpraper]
			# # except:
			# # 	continue
			# print(gini_stat,m20_stat,gini_lprnoaper,m20_lprnoaper,gini_lpraper,m20_lpraper)#
	col0 = fits.Column(name='id', array = idx_list, format = 'K')
	col1 = fits.Column(name='Moment_20_statmorph', array=gini_m20_statmorph[:,0], format='D')
	col2 = fits.Column(name='Gini_coeff_statmorph', array = gini_m20_statmorph[:,1], format='D')
	# col3 = fits.Column(name='Moment_20_noaper', array=gini_m20_lprnoaper[:,0], format='D')
	# col4 = fits.Column(name='Gini_coeff_noaper', array = gini_m20_lprnoaper[:,1], format='D')
	# col5 = fits.Column(name='Moment_20_aper', array=gini_m20_lpraper[:,0], format='D')
	# col6 = fits.Column(name='Gini_coeff_aper', array = gini_m20_lpraper[:,1], format='D')
	hdu = fits.BinTableHDU.from_columns([col0,col1,col2])#,col3,col4,col5,col6
	# print(image_path + field + '_GiniM20_10arcseccutout_1petro.fits')
	# hdu.writeto(image_path + field + '_GiniM20_10arcseccutout_1petro.fits',overwrite=True)
	print(image_path + field + '_GiniM20_regulars_13petro_3detectsigma_'+str(eta)[-1]+'eta.fits')
	hdu.writeto(image_path + field + '_GiniM20_regulars_13petro_3detectsigma_'+str(eta)[-1]+'eta.fits',overwrite=True)
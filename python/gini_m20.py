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
import time
from scipy import ndimage as ndi
from shutil import copy
import subprocess as sp

fields = ['egs']#'goodsn','goodss',
psf_path = {'goodsn':'goodsn','goodss':'goodss','egs':'aegis'}
image_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
ctg_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
segm_path = '/Users/lpr/Data/lirg_project/output/gini_m20_segmap/'
band_used = 'f160w'
even_check = {False:0,True:1}
# detection_threshold = 3
# mtotal_box = 5
# npixels = 5
aper = 17 #设置计算gini和m20的孔径为1个角秒，大概是17个像素,0.5"=8.3pixels
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
		if idx !=32411:#idx==2415:#
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
			# |根据10角秒的cutout，跑一遍sextractor得到segmentation|
			# os.chdir(segm_path)
			# sp.run('sex '+path_1 + '/' + field + '_f160w_' + str(idx) + '.fits'+' -c default.sex -CATALOG_NAME '+field+'_'+str(idx)+'.cat -CHECKIMAGE_NAME '+field+'_'+str(idx)+'.fits',shell=True,check=True)# -DETECT_MINAREA 21 -DETECT_THRESH 5
			# # |把计算半光半径的程序生成的segmentation拿来直接用|
			# segmap = fits.getdata('/Users/lpr/Data/lirg_project/output/van_re_half_light_radius/segmentation/'+str(idx)+'.fits',0)
			# |把用10角秒的图像进行sextractor得到的segmentation拿来直接用|
			segmap = fits.getdata('/Users/lpr/Data/lirg_project/output/gini_m20_segmap/'+field+'_'+str(idx)+'.fits',0)
			# # |直接用photutils做segmap|
			# threshold = photutils.detect_threshold(image, 3)
			# npixels = 5
			# segm = photutils.detect_sources(image, threshold, npixels)
			# segmap = segm.data
			# label = np.argmax(segm.areas) + 1
			# segmap_1 = segm.data == label
			# segmap_float = ndi.uniform_filter(np.float64(segmap_1), size=10)
			# segmap_1 = np.int64(segmap_float > 0.5)
			# plt.imshow(segmap,cmap='binary')
			# plt.savefig('/Users/lpr/Desktop/temp.png')
			# |用statmorph直接出结果|
			try:
				image = image - noise
				source_morphs = statmorph.source_morphology(image,segmap,gain=1,psf=None)
				morph = source_morphs[segmap[int(x),int(y)]-1]
				gini_stat = morph.gini
				m20_stat = morph.m20
				gini_m20_statmorph[num1] = [m20_stat,gini_stat]
				print(gini_m20_statmorph[num1])
			except:
				continue
			# # |扣背景和污染源|
			# image[(segmap != segmap[int(x),int(y)])&(segmap != 0)] = noise_map[(segmap != segmap[int(x),int(y)])&(segmap != 0)]
			# image = image - noise
			# # plt.imshow(image,cmap='binary')
			# # plt.savefig(segm_path+field+'_'+str(idx)+'_sex.png')
			# # |整理segmap，是目标源的地方赋值label=1，其余地方赋值为0|
			# segmap_copy = np.copy(segmap)
			# segm = np.where(segmap_copy==segmap[int(x),int(y)],1,0)
			# # |先做petrosianr, 然后smooth 0.2*petrosianr, 然后再算petrosianr和petrosianr的面亮度, 再定义segmap|
			# center = mp.center(image,segm,1)
			# petror = mp.petror_ellip(image,center,amax=a,bmax=b,theta=theta)
			# kernel = Gaussian2DKernel(x_stddev=0.2*petror,y_stddev=0.2*petror,x_size=image.shape[0]-even_check[(image.shape[0] % 2) == 0],y_size=image.shape[1]-even_check[(image.shape[1] % 2) == 0])
			# image = convolve(image,kernel)
			# try:
			# annu = ea([int(x),int(y)],1.5*(morph.rpetro_ellip)-1,1.5*(morph.rpetro_ellip)+1,(1.5*(morph.rpetro_ellip)+1)*(1-morph.ellipticity_centroid),theta=morph.orientation_centroid)
			# flux_petro = photutils.aperture.aperture_photometry(image,annu,method='exact',subpixels=5)
			# mu_petro = flux_petro['aperture_sum']/annu.area
			# segmap_copy_temp = np.zeros_like(segmap_copy)
			# segmap_copy_temp[image > mu_petro] = 1
			# segmap_copy &= segmap_copy_temp
			# try:
			# # except:
			# # 	continue
			# # |用自己写的程序，不设置aper|
			# # gini_lprnoaper = mp.gini(image,segmap,segmap[int((image.shape[0]-1)/2),int((image.shape[1]-1)/2)])
			# # m20_lprnoaper = mp.mtwenty(mp.mtotal(image,x,y,mtotal_box,segmap,segmap[int((image.shape[0]-1)/2),int((image.shape[1]-1)/2)])[2],image,x,y,segmap,segmap[int((image.shape[0]-1)/2),int((image.shape[1]-1)/2)])
			# 	gini_lprnoaper = mp.gini(image,petror,center)
			# 	m20_lprnoaper = mp.mtwenty(image,petror,center)
			# 	gini_m20_lprnoaper[num1] = [m20_lprnoaper,gini_lprnoaper]
			# # |用自己写的程序，设置aper=17个像素|
			# # gini_lpraper = mp.gini(image,segmap,segmap[int((image.shape[0]-1)/2),int((image.shape[1]-1)/2)],aper=aper,x_centroid=x,y_centroid=y)
			# # m20_lpraper = mp.mtwenty(mp.mtotal(image,x,y,mtotal_box,segmap,segmap[int((image.shape[0]-1)/2),int((image.shape[1]-1)/2)],aper=aper,x_centroid=x,y_centroid=y)[2],image,x,y,segmap,segmap[int((image.shape[0]-1)/2),int((image.shape[1]-1)/2)],aper=aper)
			# 	gini_lpraper = mp.gini(image,petror,center,aper=aper)
			# 	m20_lpraper = mp.mtwenty(image,petror,center,aper=aper)
			# 	gini_m20_lpraper[num1] = [m20_lpraper,gini_lpraper]
			# except:
			# 	continue
			# print(gini_stat,m20_stat,gini_lprnoaper,m20_lprnoaper)#,gini_lpraper,m20_lpraper)
			print('---------')
		et = time.process_time()
		print(field+' '+str(idx)+'time consume: ' + str(et - st))
	col0 = fits.Column(name='id', array = idx_list, format = 'K')
	col1 = fits.Column(name='Moment_20_statmorph', array=gini_m20_statmorph[:,0], format='D')
	col2 = fits.Column(name='Gini_coeff_statmorph', array = gini_m20_statmorph[:,1], format='D')
	# col3 = fits.Column(name='Moment_20_noaper', array=gini_m20_lprnoaper[:,0], format='D')
	# col4 = fits.Column(name='Gini_coeff_noaper', array = gini_m20_lprnoaper[:,1], format='D')
	# col5 = fits.Column(name='Moment_20_aper', array=gini_m20_lpraper[:,0], format='D')
	# col6 = fits.Column(name='Gini_coeff_aper', array = gini_m20_lpraper[:,1], format='D')
	hdu = fits.BinTableHDU.from_columns([col0,col1,col2])#,col3,col4,col5,col6
	print(image_path + field + '_GiniM20_10arcseccutout_sex_nomask.fits')
	hdu.writeto(image_path + field + '_GiniM20_10arcseccutout_sex_nomask.fits',overwrite=True)
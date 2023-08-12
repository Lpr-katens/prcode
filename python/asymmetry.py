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

fields = ['goodsn','goodss','egs']#
psf_path = {'goodsn':'goodsn','goodss':'goodss','egs':'aegis'}
image_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
ctg_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
segm_path = '/Users/lpr/Data/lirg_project/output/gini_m20_segmap/'
band_used = 'f160w'
even_check = {False:0,True:1}
aper = 17 #设置计算gini和m20的孔径为1个角秒，大概是17个像素,0.5"=8.3pixels
for field in fields:
	psf = fits.open('/Users/lpr/Data/lirg_project/intake/CANDELS/'+psf_path[field]+'_3dhst_v4.0_wfc3_psf/'+psf_path[field]+'_3dhst.v4.0.F160W_psf.fits')[0].data
	path_1 = image_path + field + '/' + field + '_' + band_used
	id_ctg = fits.getdata(ctg_path+field+'_id.fits',1)
	idx_list = np.full(len(id_ctg),-99)
	asym_statmorph = np.full((len(id_ctg),1),-999.)
	h_ctg = fits.open(ctg_path+field+'_Huangall_candels_radec_van_zmasssfr8sfrhuangsfr.fits')[1].data
	for num1 in range(0,len(id_ctg)):#
		idx = id_ctg[num1]['id']
		st = time.process_time()
		if idx != 32411:#idx==2415:#
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
			# # |根据10角秒的cutout，跑一遍sextractor得到segmentation|
			# os.chdir(segm_path)
			# sp.run('sex '+path_1 + '/' + field + '_f160w_' + str(idx) + '.fits'+' -c default.sex -CATALOG_NAME '+field+'_'+str(idx)+'.cat -CHECKIMAGE_NAME '+field+'_'+str(idx)+'.fits',shell=True,check=True)# -DETECT_MINAREA 21 -DETECT_THRESH 5
			# |把计算半光半径的程序生成的segmentation拿来直接用|
			segmap = fits.getdata('/Users/lpr/Data/lirg_project/output/van_re_half_light_radius/segmentation/'+str(idx)+'.fits',0)
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
			# |整理segmap，是目标源的地方赋值label=1，其余地方赋值为0|
			segmap_copy = segmap
			segmap_copy[segmap!=segmap[int(x),int(y)]] = 0
			segm = np.where(segmap_copy==segmap[int(x),int(y)],1,0)
			# |用statmorph直接出结果|
			try:
				source_morphs = statmorph.source_morphology(image,segm,gain=1,psf=None)
				morph = source_morphs[0]
				asym_statmorph[num1] = morph.asymmetry
			except:
				continue
	# 		print('---------')
	# 	et = time.process_time()
	# 	print(field+' '+str(idx)+'time consume: ' + str(et - st))
	col0 = fits.Column(name='id', array = idx_list, format = 'K')
	col1 = fits.Column(name='asymmetry',array=asym_statmorph,format='D')
	hdu = fits.BinTableHDU.from_columns([col0,col1])#,col2,col3,col4,col5,col6
	# print(image_path + field + '_GiniM20_10arcseccutout.fits')
	hdu.writeto(image_path + field + '_asymmetry.fits',overwrite=True)
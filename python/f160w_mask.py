# =====================================================================================
# Do sextractor mask for lirg project three fields in f160w images.
# =====================================================================================
import os
from shutil import copy
import subprocess as sp
fields_list = ['goodsn','goodss','egs']
path = '/Users/lpr/Data/SExtractor/'
for field in fields_list:
	temp = os.listdir(path+field+'_f160w')
	temp_list = []
	for file in temp:
		if file.endswith('.fits'):
			temp_list.append(file)
	file = open(path+field+'_f160w/seg.param','w')
	file.write(
		'X_IMAGE\n'
		'Y_IMAGE\n'
		)
	file.close()
	os.chdir(path+field+'_f160w')
	print(os.getcwd())
	for num in range(0,len(temp_list)):
		idx = temp_list[num][temp_list[num].rindex('_')+1:temp_list[num].index('.fits')]
		file = open(path+field+'_f160w/'+idx+'.sex','w')
		file.write(
			'CATALOG_NAME     '+idx+'_cat.fits\n'
			'CATALOG_TYPE     FITS_1.0\n'
			'PARAMETERS_NAME  seg.param\n'
			'DETECT_TYPE      CCD\n'
			'DETECT_MINAREA   5\n'
			'DETECT_THRESH    2.5\n'
			'ANALYSIS_THRESH  2.5\n'
			'FILTER           N\n'
			'FILTER_NAME      default.conv\n'
			'DEBLEND_NTHRESH  32\n'
			'DEBLEND_MINCONT  0.005\n'
			'CLEAN            Y\n'
			'CLEAN_PARAM      1.0\n'
			'MASK_TYPE        NONE\n'
			'SEEING_FWHM      0.18\n'
			'STARNNW_NAME     /usr/local/astromatic/SEx/share/sextractor/default.nnw\n'
			'BACK_SIZE        90\n'
			'BACK_FILTERSIZE  3\n'
			'BACKPHOTO_TYPE   LOCAL\n'
			'CHECKIMAGE_TYPE  SEGMENTATION\n'
			'CHECKIMAGE_NAME  '+idx+'_seg.fits\n'
			)
		file.close()
		# copy(psf,path+temp_list[num]+'/test.psf')
		line = 'sex '+temp_list[num]+' -c '+idx+'.sex'
		sp.run(line,shell=True,check=True)
		print('------------------------ '+idx+' id done ------------------------')	
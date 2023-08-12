from astropy.io import fits
import numpy as np
import os
import subprocess as sp

# # load bias file
# bias_filelist = os.listdir('/Users/lpr/Data/fits/pridata/xin_jiang/BIAS')
# bias_list = []
# for file in bias_filelist:
# 	if file.endswith('.fits'):
# 		bias_list.append(file)

# load flat file
flat_filelist = os.listdir('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT')
flat_list = []
for file in flat_filelist:
	if file.endswith('.fits') and file[0:6] != 'median':
		flat_list.append(file)
flat_b = []
flat_v = []
flat_r = []
flat_i = []
for file in flat_list:
	if file[file.index('_')+1]=='B':
		flat_b.append(file)
	elif file[file.index('_')+1]=='V':
		flat_v.append(file)
	elif file[file.index('_')+1]=='R':
		flat_r.append(file)
	elif file[file.index('_')+1]=='I':
		flat_i.append(file)
print('B band have:'+str(len(flat_b))+'\nV band have:'+str(len(flat_v))+'\nR band have:'+str(len(flat_r))+'\nI band have:'+str(len(flat_i)))
print('bias and flat file load done')

# # #=====================================
# # # bias median(remove cosmic radiation)
# # #=====================================

# # open fits file as ndarray
# bias_fits = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/BIAS/' + bias_list[0])[0]
# bias_fits_array = bias_fits.data
# bias_array = np.zeros([len(bias_list),bias_fits_array.shape[0],bias_fits_array.shape[1]])
# # store every bias array in bias_array
# for num in range(0,len(bias_list)):
#  	bias_array[num] = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/BIAS/' + bias_list[num])[0].data
# print('all bias fits are in bias_array')
# bias_median = np.zeros(bias_fits_array.shape)
# # rank values in every pixel
# for x in range(0,bias_fits_array.shape[0]):
#  	for y in range(0,bias_fits_array.shape[1]):
#          value = np.zeros(len(bias_list))
#          for num in range(0,len(bias_list)):
#              value[num] = bias_array[num,x,y]
#              bias_median[x,y] = np.median(value)
# print('get the median of bias')
# # write bias fits file
# bias_fits.data = bias_median

# bias_fits.writeto('/Users/lpr/Data/fits/pridata/xin_jiang/BIAS/median_bias.fits')
# print('write bias median in bias_fits')

#=====================================
# eliminate sources of sky flat
#=====================================

# # SExtractor: remove sources
# for num in range(0,len(flat_list)):
#  	os.chdir('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_SEx')
#  	fits_file = '/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/'+flat_list[num]
#  	sex_file = 'eliminate_sou.sex'
#  	cat_name = ' -CATALOG_NAME sou_'+flat_list[num][5:flat_list[num].index('.fits')]+'.cat'
#  	chima_name = ' -CHECKIMAGE_NAME sou_'+flat_list[num][5:]
#  	SEx = 'sex ' + fits_file + ' -c ' + sex_file + cat_name + chima_name	
#  	sp.run(SEx,shell=True,check=True)
# print('SExtractor program done')

# make pixel value = 0, where segmentation pixel value of SExtractor is greater than 0
for num in range(0,len(flat_list)):
    flat_array = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/' + flat_list[num])[0]
    seg_flat = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_SEx/sou_'+flat_list[num][5:])[0]
    flat_array_float = flat_array.data.astype(float)
    flat_array_float[np.where(seg_flat.data>0)]= np.nan
    flat_array.data = flat_array_float
    flat_array.writeto('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/processed_'+flat_list[num][5:])
print('sources removement done')
#---------------------------------------------
# store every flat_b array in flat_b_array
flat_b_fits = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/processed_'+flat_b[0][5:])[0]
flat_b_fits_array = flat_b_fits.data
flat_b_array = np.zeros([len(flat_b),flat_b_fits_array.shape[0],flat_b_fits_array.shape[1]])
# store every flat_b array in flat_b_array
# falt need to divide exposure time
for num in range(0,len(flat_b)):
 	flat_temp = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/processed_' + flat_b[num][5:])[0].data
 	exptime = flat_b[num][flat_b[num].index('_',20)+1:flat_b[num].index('s_')]
 	flat_b_array[num] = flat_temp/int(exptime)

# median values in every pixel
flat_b_median = np.nanmedian(flat_b_array,axis=0)
print('get the median of flat_b')

# write flat fits file
flat_b_fits.data = flat_b_median
flat_b_fits.writeto('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/median_flat_b.fits')
print('write flat median in flat_b_fits')

# store every flat_v array in flat_v_array
flat_v_fits = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/processed_'+flat_v[0][5:])[0]
flat_v_fits_array = flat_v_fits.data
flat_v_array = np.zeros([len(flat_v),flat_v_fits_array.shape[0],flat_v_fits_array.shape[1]])
# store every flat_v array in flat_v_array
# falt need to divide exposure time
for num in range(0,len(flat_v)):
 	flat_temp = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/processed_' + flat_v[num][5:])[0].data
 	exptime = flat_v[num][flat_v[num].index('_',20)+1:flat_v[num].index('s_')]
 	flat_v_array[num] = flat_temp/int(exptime)

# value values in every pixel
flat_v_median = np.nanmedian(flat_v_array,axis=0)
print('get the median of flat_v')

# write flat fits file
flat_v_fits.data = flat_v_median
flat_v_fits.writeto('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/median_flat_v.fits')
print('write flat median in flat_v_fits')

# store every flat_r array in flat_r_array
flat_r_fits = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/processed_'+flat_r[0][5:])[0]
flat_r_fits_array = flat_r_fits.data
flat_r_array = np.zeros([len(flat_r),flat_r_fits_array.shape[0],flat_r_fits_array.shape[1]])
# store every flat_r array in flat_r_array
# falt need to divide exposure time
for num in range(0,len(flat_r)):
 	flat_temp = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/processed_' + flat_r[num][5:])[0].data
 	exptime = flat_r[num][flat_r[num].index('_',20)+1:flat_r[num].index('s_')]
 	flat_r_array[num] = flat_temp/int(exptime)

# median values in every pixel
flat_r_median = np.nanmedian(flat_r_array,axis=0)
print('get the median of flat_r')

# write flat fits file
flat_r_fits.data = flat_r_median
flat_r_fits.writeto('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/median_flat_r.fits')
print('write flat median in flat_r_fits')

# store every flat_i array in flat_i_array
flat_i_fits = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/processed_'+flat_i[0][5:])[0]
flat_i_fits_array = flat_i_fits.data
flat_i_array = np.zeros([len(flat_i),flat_i_fits_array.shape[0],flat_i_fits_array.shape[1]])
# store every flat_i array in flat_i_array
# falt need to divide exposure time
for num in range(0,len(flat_i)):
 	flat_temp = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/processed_' + flat_i[num][5:])[0].data
 	exptime = flat_i[num][flat_i[num].index('_',20)+1:flat_i[num].index('s_')]
 	flat_i_array[num] = flat_temp/int(exptime)

# median values in every pixel
flat_i_median = np.nanmedian(flat_i_array,axis=0)
print('get the median of flat_i')

# write flat fits file
flat_i_fits.data = flat_i_median
flat_i_fits.writeto('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/flat_sEx/median_flat_i.fits')
print('write flat median in flat_i_fits')
# #--------------------------------------------
# #--------------------------------------------
# #--------------------------------------------
# # median directly rather than removing the sources(wrong)
# # store every flat_b array in flat_b_array
# flat_b_fits = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/'+flat_b[0])[0]
# flat_b_fits_array = flat_b_fits.data
# flat_b_array = np.zeros([len(flat_b),flat_b_fits_array.shape[0],flat_b_fits_array.shape[1]])
# # store every flat_b array in flat_b_array
# # falt need to divide exposure time
# for num in range(0,len(flat_b)):
#  	flat_temp = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/' + flat_b[num])[0].data
#  	exptime = flat_b[num][flat_b[num].index('_',20)+1:flat_b[num].index('s_')]
#  	flat_b_array[num] = flat_temp/int(exptime)

# # median values in every pixel
# flat_b_median = np.median(flat_b_array,axis=0)
# print('get the median of flat_b')

# # write flat fits file
# flat_b_fits.data = flat_b_median
# flat_b_fits.writeto('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/median_flat_b.fits')
# print('write flat median in flat_b_fits')

# # store every flat_v array in flat_b_array
# flat_v_fits = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/'+flat_v[0])[0]
# flat_v_fits_array = flat_v_fits.data
# flat_v_array = np.zeros([len(flat_v),flat_v_fits_array.shape[0],flat_v_fits_array.shape[1]])
# # store every flat_v array in flat_v_array
# # falt need to divide exposure time
# for num in range(0,len(flat_v)):
#  	flat_temp = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/' + flat_v[num])[0].data
#  	exptime = flat_v[num][flat_v[num].index('_',20)+1:flat_v[num].index('s_')]
#  	flat_v_array[num] = flat_temp/int(exptime)

# # median values in every pixel
# flat_v_median = np.median(flat_v_array,axis=0)
# print('get the median of flat_v')

# # write flat fits file
# flat_v_fits.data = flat_v_median
# flat_v_fits.writeto('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/median_flat_v.fits')
# print('write flat median in flat_v_fits')

# # store every flat_r array in flat_r_array
# flat_r_fits = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/'+flat_r[0])[0]
# flat_r_fits_array = flat_r_fits.data
# flat_r_array = np.zeros([len(flat_r),flat_r_fits_array.shape[0],flat_r_fits_array.shape[1]])
# # store every flat_r array in flat_r_array
# # falt need to divide exposure time
# for num in range(0,len(flat_r)):
#  	flat_temp = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/' + flat_r[num])[0].data
#  	exptime = flat_r[num][flat_r[num].index('_',20)+1:flat_r[num].index('s_')]
#  	flat_r_array[num] = flat_temp/int(exptime)

# # median values in every pixel
# flat_r_median = np.median(flat_r_array,axis=0)
# print('get the median of flat_r')

# # write flat fits file
# flat_r_fits.data = flat_r_median
# flat_r_fits.writeto('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/median_flat_r.fits')
# print('write flat median in flat_r_fits')

# # store every flat_i array in flat_i_array
# flat_i_fits = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/'+flat_i[0])[0]
# flat_i_fits_array = flat_i_fits.data
# flat_i_array = np.zeros([len(flat_i),flat_i_fits_array.shape[0],flat_i_fits_array.shape[1]])
# # store every flat_i array in flat_i_array
# # falt need to divide exposure time
# for num in range(0,len(flat_i)):
#  	flat_temp = fits.open('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/' + flat_i[num])[0].data
#  	exptime = flat_i[num][flat_i[num].index('_',20)+1:flat_i[num].index('s_')]
#  	flat_i_array[num] = flat_temp/int(exptime)

# # median values in every pixel
# flat_i_median = np.median(flat_i_array,axis=0)
# print('get the median of flat_i')

# # write flat fits file
# flat_i_fits.data = flat_i_median
# flat_i_fits.writeto('/Users/lpr/Data/fits/pridata/xin_jiang/FLAT/median_flat_i.fits')
# print('write flat median in flat_i_fits')

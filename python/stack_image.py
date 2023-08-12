from astropy.io import fits
import numpy as np
import os
def rescale(arr):
    arr_min = arr.min()
    arr_max = arr.max()
    return (arr - arr_min) / (arr_max - arr_min)
path = '/Users/lpr/Desktop/image_stack/'
image_list = []
for file in os.listdir(path):
	if file.endswith('.fits'):
		image_list.append(file)
array = np.full([101,101],0,dtype='float64')
for num in range(0,len(image_list)):
	image = fits.open(path+image_list[num])[0].data
	image = rescale(image)
	array += image
print('------------------- multi-dimension array load done -------------------')
data = fits.PrimaryHDU(array)
hdu = fits.HDUList([data])
hdu.writeto(path + 'stack_f160w.fits')
print('done')
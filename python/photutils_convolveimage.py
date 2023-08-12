# # ================
# # photutils and pypher version of image_convolution
# # ================

# # f606 resolution: 0.124"(2.1pixels), f160 resolution: 0.184"(3.1pixels), goodsn pixel scale: 0.059"/pix, 0.135"= 2.288 pixels
import os
from astropy.io import fits
from photutils import CosineBellWindow, create_matching_kernel
from astropy.convolution import convolve

fileList606 = os.listdir('/Users/lpr/Data/fits/pridata/goodsn_f606')
fileList160 = os.listdir('/Users/lpr/Data/fits/pridata/goodsn_f160')
fileList_606 = []
fileList_160 = []
for items in fileList606:
    if items.endswith('.fits'):
        fileList_606.append(items)
for items in fileList160:
    if items.endswith('.fits'):
        fileList_160.append(items)

# calculate kernel
for i in range(0,len(fileList_606)):
    hdu1 = fits.open('/Users/lpr/Data/fits/pridata/goodsn_f606/'+fileList_606[i])[0]
    hdu2 = fits.open('/Users/lpr/Data/fits/pridata/goodsn_f160/goodsn_'+fileList_606[i][7:fileList_606[i].index('f606w')-1]+'_f160w.fits')[0]
    window = CosineBellWindow(alpha=1) # alpha imply the percentage of array values that are tapered(more shark)
    kernel = create_matching_kernel(hdu1.data, hdu2.data,window=window)
    astropy_conv = convolve(hdu1.data,kernel)
    hdu1.data = astropy_conv
    filename = 'photuconv_'+fileList_606[i]
    hdu1.writeto('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/photutils_convolve/'+filename)
    print(filename+' is done')

# following code should be run in terminal
#pypher /Users/lpr/Data/fits/pridata/goodsn_f606/goodsn_1603_f606w.fits /Users/lpr/Data/fits/pridata/goodsn_f160/goodsn_1603_f160w.fits /Users/lpr/Data/fits/temdata/pypher_kernel_1603.fits
# pypher /Users/lpr/Data/fits/pridata/goodsn_f606/goodsn_1092_f606w.fits /Users/lpr/Data/fits/pridata/goodsn_f160/goodsn_1092_f160w.fits /Users/lpr/Data/fits/temdata/pypher_kernel.fits
# from astropy.io import fits
# kernel = fits.open('/Users/lpr/Data/fits/temdata/pypher_kernel_1603.fits')[0]
# hdu = fits.open('/Users/lpr/Data/fits/pridata/goodsn_f606/goodsn_1603_f606w.fits')[0]
# from astropy.convolution import convolve
# astropy_conv = convolve(hdu.data, kernel.data)
# hdu.data = astropy_conv
# hdu.writeto('/Users/lpr/Data/fits/temdata/pypher_606_1603.fits')
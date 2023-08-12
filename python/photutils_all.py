#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 16:58:23 2020

@author: lpr
"""

from astropy.io import fits
from photutils import CosineBellWindow, create_matching_kernel
from astropy.convolution import convolve
import numpy as np

hdu1 = fits.open('/Users/lpr/Data/fits/pridata/goodsn/psf/PSFSTD_WFC3UV_F606W.fits',ignore_missing_end=True)
hdu2 = fits.open('/Users/lpr/Data/fits/pridata/goodsn/psf/PSFSTD_WFC3IR_F160W.fits',ignore_missing_end=True)
window = CosineBellWindow(alpha=1) # alpha imply the percentage of array values that are tapered(more shark)
h1_data = hdu1[0].data
h2_data = hdu2[0].data

# kernel = create_matching_kernel(h1_data[0], h2_data[0],window=window)
# convolution = np.array()
# for num in range(0,10):
#     astropy_conv = convolve_fft(hdu1.data[num:(num+1)*len(kernel)/10],kernel[num:(num+1)*len(kernel)/10],allow_huge=True)
#     np.insert(convolution,num*len(kernel)/10,values=astropy_conv.transpose(),axis=1)
# hdu1.data = convolution
# hdu1.writeto('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/photutils_convolve/all.fits')

# ker = fits.PrimaryHDU(kernel)
# ker_hdu = fits.HDUList([ker])
# ker_hdu.writeto('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/photutils_convolve/psf_kernel_2.fits')

# print('kernel generation is done')
#========convolve image
# ima1 = fits.open('/Users/lpr/Data/fits/pridata/goodsn_f606/goodsn_1092_f606w.fits')
# astropy_conv = convolve(h1_data[0],kernel)

# data = fits.PrimaryHDU(astropy_conv)
# hdu = fits.HDUList([data])
# hdu.writeto('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/photutils_convolve/606_psf_conv.fits')
data1 = fits.PrimaryHDU(h1_data[0])
hdu1 = fits.HDUList([data1])
hdu1.writeto('/Users/lpr/Data/fits/pridata/goodsn/psf/606_psf_conv.fits')
data2 = fits.PrimaryHDU(h2_data[0])
hdu2 = fits.HDUList([data2])
hdu2.writeto('/Users/lpr/Data/fits/pridata/goodsn/psf/160_psf_conv.fits')
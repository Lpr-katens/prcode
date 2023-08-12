#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 11:2size:53 2020

@author: lpr
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# psf = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/psf/psf_160_stackimage.fits')[0]
psf = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/multiband_cutimage/goodsn_f606w/goodsn_f606w_1092.fits')[0]
# 1size pixels in psf stack image
#===============
# annulus flux plot
size = int((psf.data.shape[0]-1)/2)
flux_list = []
for radius in range(0,size):
    flux = 0
    for x in range(0,psf.data.shape[1]):
        for y in range(0,psf.data.shape[0]):
            if radius == 0:
                if np.sqrt((x-size)**2+(y-size)**2) <= radius:
                    flux += psf.data[x][y]
            else:
                if np.sqrt((x-size)**2+(y-size)**2) >= radius-1 and np.sqrt((x-size)**2+(y-size)**2) <= radius:
                    flux += psf.data[x][y]              
    flux_list.append(flux)
    
x_axis = np.arange(0,size,1)
plt.figure()
plt.scatter(x_axis,flux_list,color='r')
plt.xlabel('radius')
plt.ylabel('flux inside annulus')
plt.title('annulus flux')

#===============
# growth curve
flux_list = []
for radius in range(0,size):
    flux = 0
    for x in range(0,psf.data.shape[1]):
        for y in range(0,psf.data.shape[0]):
            if np.sqrt((x-size)**2+(y-size)**2) <= radius:
                flux += psf.data[x][y]
    flux_list.append(flux)
    
x_axis = np.arange(0,size,1)
plt.figure()
plt.scatter(x_axis,flux_list,color='b')
plt.xlabel('radius')
plt.ylabel('flux inside radius')
plt.title('growth curve')
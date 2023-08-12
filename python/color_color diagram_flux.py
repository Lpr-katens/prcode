#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 09:35:55 2020

@author: lpr
"""

# ================
# try to plot color-color diagram
# ================

import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# load the topcat match file
fileList = os.listdir('/Users/lpr/Data/fits/pridata/HSC_wide/new_hsc')
file = []
for names in fileList:
    if names.endswith('.fits') and names[names.index('.fits')-2:names.index('.fits')] != 'ap':
        file.append(names)
print('file load done')

# define a division function, in case divisor=0
def divi(p,r):
    if r!= 0:
        return p/r
    elif r == 0:
        return float('inf')

# define a function to calculate signal to noise ratio
def sigtonoi(s,n):
    snr = []
    for num in range(0,len(s)):
        snr.append(divi(s[num],n[num]))
    return snr

g = []
g_error = []
r = []
r_error = []
i = []
i_error = []
z = []
z_error = []
y = []
y_error = []
for num in range(0,len(file)):
    # fits.open()[1] is fits binary table, fits.open()[0] is fits header
    hdu = fits.open('/Users/lpr/Data/fits/pridata/HSC_wide/new_hsc/'+file[num])[1]
    # because empty list cannot use +=, i dont want to use .append(add items to list),
    # so i make list begin with a if sentence.
    if num == 0: # this if sentence judge the list is empty(num=0) or not empty(num!=0)
        g = hdu.data['g_cmodel_flux'].tolist()
        g_error = hdu.data['g_cmodel_fluxsigma'].tolist()
        r = hdu.data['r_cmodel_flux'].tolist()
        r_error = hdu.data['r_cmodel_fluxsigma'].tolist()
        i = hdu.data['i_cmodel_flux'].tolist()
        i_error = hdu.data['i_cmodel_fluxsigma'].tolist()
        z = hdu.data['z_cmodel_flux'].tolist()
        z_error = hdu.data['z_cmodel_fluxsigma'].tolist()
        y = hdu.data['y_cmodel_flux'].tolist()
        y_error = hdu.data['y_cmodel_fluxsigma'].tolist()
        
    else:
        # if list is not empty, then use +=
        g += hdu.data['g_cmodel_flux'].tolist()
        g_error += hdu.data['g_cmodel_fluxsigma'].tolist()
        r += hdu.data['r_cmodel_flux'].tolist()
        r_error += hdu.data['r_cmodel_fluxsigma'].tolist()
        i += hdu.data['i_cmodel_flux'].tolist()
        i_error += hdu.data['i_cmodel_fluxsigma'].tolist()
        z += hdu.data['z_cmodel_flux'].tolist()
        z_error += hdu.data['z_cmodel_fluxsigma'].tolist()
        y += hdu.data['y_cmodel_flux'].tolist()
        y_error += hdu.data['y_cmodel_fluxsigma'].tolist()
print('data load done')

# calculate s/n in everyband
g_snr = sigtonoi(g,g_error)
r_snr = sigtonoi(r,r_error)
i_snr = sigtonoi(i,i_error)
z_snr = sigtonoi(z,z_error)
y_snr = sigtonoi(y,y_error)
print('signal to noise calculation done')

# g_f on behalf of g_band final list(to store all g_band data which are satisfied 
# S/N>3 ), and so on, r_f, i_f, z_f
g_f = []
r_f = []
i_f = []
z_f = []
y_f = []
# remove data with s/n<3(or in another way to say, save satisfied data in g_f,r_f,i_f,z_f)
for num in range(0,len(g)):
    if g_snr[num]>=3 and r_snr[num]>=3 and i_snr[num]>=3 and z_snr[num]>=3 and y_snr[num]>=3:
        g_f.append(g[num])
        r_f.append(r[num])
        i_f.append(i[num])
        z_f.append(z[num])
        y_f.append(y[num])
print('signal to noise selection done')

# calculate color
gr = [-2.5*np.log10(divi(a,b)) for a,b in zip(g_f,r_f)]
ri = [-2.5*np.log10(divi(a,b)) for a,b in zip(r_f,i_f)]
iz = [-2.5*np.log10(divi(a,b)) for a,b in zip(i_f,z_f)]
iy = [-2.5*np.log10(divi(a,b)) for a,b in zip(i_f,y_f)] 
print('color calculation done')

#=================================
# plot color-color diagram(there are 4 band, so i can plot two color-color diagram)
for num in range(0,3):
    if num==0: # first diagram
        plt.figure()
        plt.scatter(ri,gr, c='b', marker='o', s=0.5)
        plt.xlabel('r-i')
        plt.ylabel('g-r')
        plt.title('color-color diagram')
        plt.savefig('/Users/lpr/Data/fits/expdata/HSC/gr_ri_coco.eps')
        plt.savefig('/Users/lpr/Data/fits/expdata/HSC/gr_ri_coco.jpg')
        print('gr_ri color-color diagram plot done')
    elif num==1: # second diagram
        plt.figure()
        plt.scatter(iz,ri, c='r', marker='o', s=0.5)
        plt.xlabel('i-z')
        plt.ylabel('r-i')
        plt.title('color-color diagram')
        plt.savefig('/Users/lpr/Data/fits/expdata/HSC/ri_iz_coco.eps')
        plt.savefig('/Users/lpr/Data/fits/expdata/HSC/ri_iz_coco.jpg')
        print('ri_iz color-color diagram plot done')
    else:
        plt.figure()
        plt.scatter(iy,ri, c='g', marker='o', s=0.5)
        plt.xlabel('i-y')
        plt.ylabel('r-i')
        plt.title('color-color diagram')
        plt.savefig('/Users/lpr/Data/fits/expdata/HSC/ri_iy_coco.eps')
        plt.savefig('/Users/lpr/Data/fits/expdata/HSC/ri_iy_coco.jpg')
        print('ri_iy color-color diagram plot done')
    #====following code is used to select data with color>1.3
    #   for list in range(0,len(y)):
    #       if y[list] >= 1.3:
    #           yy.append(y[list])
    #           xx.append(specz[list])
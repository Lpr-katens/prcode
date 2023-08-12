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
    if names.endswith('.fits'):
        file.append(names)
print('file import done')

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


# specz = []
# specz_error = []
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
        # specz = hdu.data['specz_redshift'].tolist()
        # specz_error = hdu.data['specz_redshift_err'].tolist()
        g = hdu.data['g_cmodel_mag'].tolist()
        g_error = hdu.data['g_cmodel_magsigma'].tolist()
        r = hdu.data['r_cmodel_mag'].tolist()
        r_error = hdu.data['r_cmodel_magsigma'].tolist()
        i = hdu.data['i_cmodel_mag'].tolist()
        i_error = hdu.data['i_cmodel_magsigma'].tolist()
        z = hdu.data['z_cmodel_mag'].tolist()
        z_error = hdu.data['z_cmodel_magsigma'].tolist()
        y = hdu.data['y_cmodel_mag'].tolist()
        y_error = hdu.data['y_cmodel_magsigma'].tolist()
    else:
        # if list is not empty, then use +=
        # specz += hdu.data['specz_redshift'].tolist()
        # specz_error += hdu.data['specz_redshift_err'].tolist()
        g += hdu.data['g_cmodel_mag'].tolist()
        g_error += hdu.data['g_cmodel_magsigma'].tolist()
        r += hdu.data['r_cmodel_mag'].tolist()
        r_error += hdu.data['r_cmodel_magsigma'].tolist()
        i += hdu.data['i_cmodel_mag'].tolist()
        i_error += hdu.data['i_cmodel_magsigma'].tolist()
        z += hdu.data['z_cmodel_mag'].tolist()
        z_error += hdu.data['z_cmodel_magsigma'].tolist()
        y += hdu.data['y_cmodel_mag'].tolist()
        y_error += hdu.data['y_cmodel_magsigma'].tolist()
print('mag data import done')

# calculate s/n in everyband
g_snr = sigtonoi(g,g_error)
r_snr = sigtonoi(r,r_error)
i_snr = sigtonoi(i,i_error)
z_snr = sigtonoi(z,z_error)
y_snr = sigtonoi(y,y_error)
print('signal to noise done')

# g_f on behalf of g_band final list(to store all g_band data which are satisfied 
# S/N>3 ), and so on, r_f, i_f, z_f, y_f
g_f = []
r_f = []
i_f = []
z_f = []
y_f = []
# remove data with s/n<3(or in another way to say, save satisfied data in g_f,r_f,i_f,z_f)
for num in range(0,len(g)):
    if g_snr[num]>=10 and r_snr[num]>=10 and i_snr[num]>=10 and z_snr[num]>=10 and y_snr[num]>=10:
        g_f.append(g[num])
        r_f.append(r[num])
        i_f.append(i[num])
        z_f.append(z[num])
        y_f.append(y[num])

# calculate color
gr = [a-b for a,b in zip(g_f,r_f)]
ri = [a-b for a,b in zip(r_f,i_f)]
iz = [a-b for a,b in zip(i_f,z_f)]
iy = [a-b for a,b in zip(i_f,y_f)]
print('color calculation done')

#=================================
# plot color-color diagram(there are 4 band, so i can plot two color-color diagram)
# for num in range(0,3):
#     if num==0: # first diagram
#         plt.figure()
#         plt.scatter(ri,gr, c='b', marker='o', s=0.5)
#         plt.xlabel('r-i')
#         plt.ylabel('g-r')
#         plt.title('color-color diagram')
#         plt.savefig('/Users/lpr/Data/fits/expdata/HSC/gr_ri_twocolor.eps')
#     elif num==1: # second diagram
#         plt.figure()
#         plt.scatter(iz,ri, c='r', marker='o', s=0.5)
#         plt.xlabel('i-z')
#         plt.ylabel('r-i')
#         plt.title('color-color diagram')
#         plt.savefig('/Users/lpr/Data/fits/expdata/HSC/ri_iz_twocolor.eps')
#     else:
#         plt.figure()
#         plt.scatter(iy,ri, c='g', marker='o', s=0.5)
#         plt.xlabel('i-y')
#         plt.ylabel('r-i')
#         plt.title('color-color diagram')
#         plt.savefig('/Users/lpr/Data/fits/expdata/HSC/ri_iy_twocolor.eps')
for num in range(0,3):
    if num==0: # first diagram
        plt.figure()
        plt.scatter(ri,gr, c='b', marker='o', s=0.5)
        plt.xlabel('r-i')
        plt.ylabel('g-r')
        plt.title('color-color diagram')
        plt.savefig('/Users/lpr/Data/fits/expdata/HSC/gr_ri_twocolor.jpg')
    elif num==1: # second diagram
        plt.figure()
        plt.scatter(iz,ri, c='r', marker='o', s=0.5)
        plt.xlabel('i-z')
        plt.ylabel('r-i')
        plt.title('color-color diagram')
        plt.savefig('/Users/lpr/Data/fits/expdata/HSC/ri_iz_twocolor.jpg')
    else:
        plt.figure()
        plt.scatter(iy,ri, c='g', marker='o', s=0.5)
        plt.xlabel('i-y')
        plt.ylabel('r-i')
        plt.title('color-color diagram')
        plt.savefig('/Users/lpr/Data/fits/expdata/HSC/ri_iy_twocolor.jpg')
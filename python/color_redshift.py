#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 09:35:55 2020

@author: lpr
"""

# ================
# try to plot specz-color diagram
# ================

import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# load the topcat match file
fileList = os.listdir('/Users/lpr/Data/fits/expdata/HSC')
file = []
for names in fileList:
    if names.endswith('.fits'):
        file.append(names)

def divi(p,r):
    if r!= 0:
        return p/r
    elif r == 0:
        return 1e3
# define a function to calculate signal to noise ratio
def sigtonoi(s,n):
    snr = []
    for num in range(0,len(s)):
        snr.append(divi(s[num],n[num]))
    return snr

specz = []
specz_error = []
g = []
g_error = []
r = []
r_error = []
i = []
i_error = []
z = []
z_error = []
for num in range(0,len(file)):
    # fits.open()[1] is fits binary table, fits.open()[0] is fits header
    hdu = fits.open('/Users/lpr/Data/fits/expdata/HSC/'+file[num])[1]
    if num == 0:
        specz = hdu.data['specz_redshift'].tolist()
        specz_error = hdu.data['specz_redshift_err'].tolist()
        g = hdu.data['g_apertureflux_15_flux'].tolist()
        g_error = hdu.data['g_apertureflux_15_fluxsigma'].tolist()
        r = hdu.data['r_apertureflux_15_flux'].tolist()
        r_error = hdu.data['r_apertureflux_15_fluxsigma'].tolist()
        i = hdu.data['i_apertureflux_15_flux'].tolist()
        i_error = hdu.data['i_apertureflux_15_fluxsigma'].tolist()
        z = hdu.data['z_apertureflux_15_flux'].tolist()
        z_error = hdu.data['z_apertureflux_15_fluxsigma'].tolist()
    else:
        specz += hdu.data['specz_redshift'].tolist()
        specz_error += hdu.data['specz_redshift_err'].tolist()
        g += hdu.data['g_apertureflux_15_flux'].tolist()
        g_error += hdu.data['g_apertureflux_15_fluxsigma'].tolist()
        r += hdu.data['r_apertureflux_15_flux'].tolist()
        r_error += hdu.data['r_apertureflux_15_fluxsigma'].tolist()
        i += hdu.data['i_apertureflux_15_flux'].tolist()
        i_error += hdu.data['i_apertureflux_15_fluxsigma'].tolist()
        z += hdu.data['z_apertureflux_15_flux'].tolist()
        z_error += hdu.data['z_apertureflux_15_fluxsigma'].tolist()

# calculate s/n in everyband
g_snr = sigtonoi(g,g_error)
r_snr = sigtonoi(r,r_error)
i_snr = sigtonoi(i,i_error)
z_snr = sigtonoi(z,z_error)
specz_snr = sigtonoi(specz,specz_error)

# g_f on behalf of g_band final list(to store all g_band data which are satisfied 
# S/N>3 ), and so on, r_f, i_f, z_f
g_f = []
r_f = []
i_f = []
z_f = []
specz_f = []
# remove data with s/n<3(or in another way to say, save satisfied data in g_f,r_f,i_f,z_f)
for num in range(0,len(g)):
    if g_snr[num]>=3 and r_snr[num]>=3 and i_snr[num]>=3 and z_snr[num]>=3 and specz_snr[num]>=3:
        g_f.append(g[num])
        r_f.append(r[num])
        i_f.append(i[num])
        z_f.append(z[num])
        specz_f.append(specz[num])
#============= plot z-color diagram
y = [2.5*np.log10(divi(a,b)) for a,b in zip(z_f,i_f)]
plt.scatter(specz_f, y , c='b', marker='o', s=0.5)
plt.xlim(0.1,1)
plt.xlabel('redshift')
plt.ylabel('i-z')
plt.title('color-redshift diagram')
plt.savefig('/Users/lpr/Data/fits/expdata/HSC/i-z_redshift_snr3.eps')

# #============== select sources with i-z>1.3
# # yy = []
# # xx = []
# # for list in range(0,len(y)):
# #     if y[list] >= 1.3:
# #         yy.append(y[list])
# #         xx.append(specz[list])
# plt.scatter(specz, y , c='b', marker='o', s=0.5)
# plt.xlim(0.1,1)
# # plt.ylim(1.3,)
# plt.xlabel('redshift')
# plt.ylabel('i-z(>1.3)')
# plt.title('color-redshift diagram')
# plt.savefig('/Users/lpr/Data/fits/expdata/HSC/i-z_redshift_1.eps')



















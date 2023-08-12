#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 15:32:05 2020

@author: lpr
"""

# pixel sacel of all image in HST goodsn is 0.059"
# ================
# plot FWHM-magnitude diagram
# ================

import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess as sp

# get all HST goodsn image
image_list = [] # all wavelength image of whole field
image_folder_file = os.listdir('/Users/lpr/Data/fits/pridata/goodsn')
for items in image_folder_file:
    if items.endswith('_drz.fits'):
        image_list.append(items)
print('goodsn file load done')

os.chdir('/Users/lpr/Data/SExtractor/test')
print('current word directory is: '+os.getcwd())

for num in range(0,len(image_list)):
    # run SExtractor command in terminal
    band = image_list[num][image_list[num].index('_f')+1:image_list[num].index('_060mas')]
    comline = '/usr/local/astromatic/SEx/bin/sex /Users/lpr/Data/fits/pridata/goodsn/'+image_list[num]+' -c default.sex -CATALOG_NAME '+band+'.cat -PARAMETERS_NAME sigma.param -CHECKIMAGE_NAME '+band+'_seg.fits'
    sp.run(comline,shell=True,check=True)
    print(band+' SExtractor command execution done')
    # open sextractor results file
    fwhm_mag = open(band+'.cat')
    wholefile = fwhm_mag.readlines() # including first 2 lines
    del wholefile[:3] # remove first 2 lines
    fwhm = [] # store all fwhm in variable fwhm
    mag = [] # store all magnitude in variable mag
    # find FWHM index in every line
    for num in range(0,len(wholefile)):
        temp1 = wholefile[num].split()[0]
        temp2 = wholefile[num].split()[1]
        fwhm.append(float(temp1))
        mag.append(float(temp2))         
    print('fwhm-mag list append done')
    
    # plot FWHM-magnitude diagram
    plt.figure()
    fig=plt.subplot(1,1,1)
    fig.scatter(mag,fwhm,color='b',marker='o',s=0.5)
    plt.xlim(-15,10)
    plt.ylim(-4,15)
    x = np.linspace(-15,10,30)
    y = np.linspace(-4,15,30)
    plt.xticks(x,rotation=90)
    plt.yticks(y)
    plt.xlabel('mag')
    plt.ylabel('FWHM')
    plt.title('FWHM-magnitude')
    fig.grid(color='r',linewidth=0.1)
    plot_name = 'fwhm_mag_'+band
    plt.savefig('/Users/lpr/Data/fits/expdata/HST/goodsn_all/SEx_select_source/'+plot_name+'.eps')
    plt.savefig('/Users/lpr/Data/fits/expdata/HST/goodsn_all/SEx_select_source/'+plot_name+'.jpg')
    print('all steps of '+band+' done')

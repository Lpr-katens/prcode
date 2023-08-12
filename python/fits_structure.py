#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 09:22:51 2020

@author: lpr
"""
from astropy.io import fits

hdulist = fits.open('/Users/lpr/Data/fits/pridata/goodsn_f606/goodsn_597_f606w.fits')
hdulist.info()
prihdr = hdulist[0].header
print(prihdr) #print all fits headers, if want a specific header information, use index or string, for example ['history'].

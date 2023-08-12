#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 15:31:42 2020

@author: lpr
"""

# ================
# plot I-Z--Z and R-I--Z diagram
# ================

import os
import matplotlib.pyplot as plt

# use 'for' to get all fits file in list filelist
filelist=[]
for names in os.listdir('/Users/lpr/Data/fits/pridata/HSC-wide/S18A_starmask'):
  if names.endswith(".txt"):
    filelist.append(names)

r = []
g = []
z = []
i = []
for k in range(0,len(filelist)):
    print(filelist[k])
    hdu = open('/Users/lpr/Data/fits/pridata/HSC-wide//S18A_starmask/'+filelist[k])
    wholefile = hdu.readlines() # including first 6 lines
    neededfile = []
    if 'g' in filelist[k]:
        g.append([])
        g.append([])
    elif 'r' in filelist[k]:
        r.append([])
        r.append([])
    elif 'i' in filelist[k]:
        i.append([])
        i.append([])
    elif 'z' in filelist[k]:
        z.append([])
        z.append([])
    for j in range(0,len(wholefile)):
        if j > 5 and j%3 == 0:
            neededfile.append(wholefile[j]) # write needed data in neededfile
    for j in range(0,len(neededfile)):
        index1 = neededfile[j].index('DR2')
        index2 = neededfile[j].index(', mag')
        if 'r' in filelist[k]:
            # r[1] is objects' index; r[2] is objects' magnitude
            r[0].append(int(neededfile[j][index1+4:index2]))
            r[1].append(float(neededfile[j][index2+6:]))
        elif 'g' in filelist[k]:
            g[0].append(int(neededfile[j][index1+4:index2]))
            g[1].append(float(neededfile[j][index2+6:]))
        elif 'z' in filelist[k]:
            z[0].append(int(neededfile[j][index1+4:index2]))
            z[1].append(float(neededfile[j][index2+6:]))
        elif 'i' in filelist[k]:
            i[0].append(int(neededfile[j][index1+4:index2]))
            i[1].append(float(neededfile[j][index2+6:]))

x = [] # horizontal axis data in x
y = [] # vertical axis data in y

# match r, g with objects' index
for k in range(0,len(r[0])):
    if r[0][k] in g[0]:
        x.append(r[1][k])
        y.append(g[1][g[0].index(r[0][k])] - r[1][k])
        
# plot g-r - r diagram
plt.plot(x,y,color='red')
plt.xlabel('r_band magnitude')
plt.ylabel('g-r magnitude')
plt.title('g-r - r diagram')

# match i, z with objects' index
for k in range(0,len(z[0])):
    if z[0][k] in i[0]:
        x.append(z[1][k])
        y.append(i[1][i[0].index(z[0][k])] - z[1][k])
        
# plot g-r - r diagram
plt.plot(x,y,color='red')
plt.xlabel('z_band magnitude')
plt.ylabel('i-z magnitude')
plt.title('i-z - z diagram')

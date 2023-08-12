import numpy as np
from astropy.modeling.models import Sersic2D
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.io import fits
# from photutils.datasets import make_noise_image # if need noise
from scipy.integrate import quad
from scipy.special import gammaincinv
import os

print('         *              *\n'
	  '        ***            ***\n'
	  '       *****          *****\n'
	  '      **********************\n'
	  '    **************************\n'
	  '   **********  code  **********\n'
	  '   *********   begin  *********\n'
	  '    ********   now    ********\n'
	  '      **********************                         ***\n'
	  '        ******************                       **********\n'
	  '          **************                    ******************\n'
	  '            **********              ******************************\n'
	  '            *********************************************************\n'
	  '          ************************************************************\n'
	  '        ***************************************************************\n'
	  '      ******************************************************************\n'
	  '      *******      *******                        *******     *******\n'
	  '       *****        *****                          *****       *****\n'
	  '        ***          ***                            ***         ***\n'
	  '         *            *                              *           *\n')

def Intensity_r(r,re,Ire,n):
	# sersic profile function
	# input:
	# 	Ire: surface brightness in re
	# 	re: effective radius
	# 	n: sersic index
	# 	r: radius
	bn = gammaincinv(2.*n,0.5)
	Intensity_r = Ire*np.exp(-bn*((r/re)**(1/n)-1))
	return Intensity_r

def integral_Intensity_r(r,re,Ire,n):
	# i need to do some adjustment to calculate luminosity of sersic profile
	# input:
	# 	Ire: surface brightness in re
	# 	re: effective radius
	# 	n: sersic index
	# 	r: radius
	bn = gammaincinv(2.*n,0.5)
	L_r = 2*np.pi*r*Ire*np.exp(-bn*((r/re)**(1/n)-1))
	return L_r

def sersic_profile(Ire,re,n,x0,y0,x,y):
	# calculate surface brightness at radius r
	# input:
	# 	Ire: surface brightness in re
	# 	re: effective radius
	# 	n: sersic index
	# 	x0,y0: galaxy center(in pixel)
	# 	x,y: position you want to calculate(in pixel)
	r = np.sqrt((x-x0)**2+(y-y0)**2)
	Ir = Intensity_r(r,re,Ire,n)
	return Ir

def sersic_luminosity(n,re,Ire):
	# calculate total luminosity at radius r
	# input:
	# 	Ire: surface brightness in re
	# 	re: effective radius
	# 	n: sersic 
	L_total = quad(integral_Intensity_r,0,float('inf'),args=(re,Ire,n))[0]
	return L_total

def modify_center_value(shape,Ire,re,n):
	# modify the value of Sersic2D center pixel
	# Because one CCD reciever can't tell 
	# the real sersic distribution of galaxy.
	# So, we need to modify center pixel value to
	# justify CCD condition.
	# input:
	# 	shape: cut center pixel into subpixel of shape, elements in shape should be odd
	# 	Ire: surface brightness in re
	# 	re: effective radius
	# 	n: sersic index
	temp = np.zeros(shape)
	for x in range(0,shape[0]):
		for y in range(0,shape[1]):
			temp[x,y] = sersic_profile(Ire,re,n,0.5,0.5,x/(shape[0]-1),y/(shape[1]-1))
	return np.mean(temp)

def effective_radius_total(hdu,x0,y0):
	# Calculate the effective radius of input image array.(in pixels)
	# input:
	# 	hdu: cut center pixel into subpixel of shape, elements in shape should be odd
	# 	x0,y0: galaxy center(in pixels)
	flux = 0
	for re in np.arange(1,np.sqrt(((hdu.shape[0]-1)/2)**2+((hdu.shape[1]-1)/2)**2),1):
		if flux >= 0.5:
			print(flux)
			break
		for x in range(int((hdu.shape[0]-1)/2),hdu.shape[0]):
			for y in range(int((hdu.shape[1]-1)/2),hdu.shape[1]):
				if np.sqrt((x-x0)**2+(y-y0)**2)>re-1 and np.sqrt((x-x0)**2+(y-y0)**2)<=re:
					if x == int((hdu.shape[0]-1)/2) and y == int((hdu.shape[1]-1)/2):
						flux += hdu[x,y]
					elif x == int((hdu.shape[0]-1)/2) or y == int((hdu.shape[1]-1)/2):
						flux += hdu[x,y]*2
					else:
						flux += hdu[x,y]*4
	return re-1

export_path = '/Users/lpr/Data/fits/expdata/BTratio/'
# export name in order of B/T,ellip,devau-re,exp-re,re_total
shape = (101,101) # image size 2001x2001 pixels
# noise_mean=5.0 # with no noise
x,y = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
# ----------------------------------------------------------------------
# --------------------------- galaxy profile ---------------------------
# ----------------------------------------------------------------------
RE_surbri = 10 # surface brightness in Re
brightness_show = 100 # set brightness to more visible
ellip_step = 0.1
pa_step = 2*np.pi/40
re_devau_size = [41,1]
re_exp_size = [41,1]
# for re_devau in np.arange(1,re_devau_size[0],re_devau_size[1]): # --------------------------- de vaucouleurs profile ---------------------------
# 	model_devau = Sersic2D(amplitude = RE_surbri,r_eff = re_devau,n=4,x_0=int((shape[0]-1)/2),y_0=int((shape[1]-1)/2))
# 	image_devau = model_devau(x,y)
# 	image_devau[int((shape[0]-1)/2),int((shape[1]-1)/2)] = modify_center_value([51,51],RE_surbri,re_devau,4)
# 	image_devau = image_devau/np.sum(image_devau) # normalization to total = 1
# 	hdu_bulge = fits.PrimaryHDU(image_devau)
# 	name_bulge = export_path+'bulge/'+str(np.around(re_devau,1))
# 	hdu_bulge.writeto(name_bulge+'.fits', overwrite=True)
# 	print(str(np.around(re_devau,1))+'bulge.fits export done')
for ellip in np.arange(0,1.1,ellip_step):
	for pa in np.arange(0,2*np.pi,pa_step):
		for re_exp in np.arange(1,re_exp_size[0],re_exp_size[1]): # --------------------------- de vaucouleurs profile ---------------------------
			model_exp = Sersic2D(amplitude = RE_surbri,r_eff = re_exp,n=1,x_0=int((shape[0]-1)/2),y_0=int((shape[1]-1)/2),ellip=ellip,theta=pa)
			image_exp = model_exp(x,y)
			image_exp[int((shape[0]-1)/2),int((shape[1]-1)/2)] = modify_center_value([51,51],RE_surbri,re_exp,1)
			image_exp = image_exp/np.sum(image_exp) # normalization to total = 1
			hdu_disk = fits.PrimaryHDU(image_exp)
			name_disk = export_path+'disk/'+str(np.around(ellip,1))+'-'+str(np.around(pa,3))+'-'+str(np.around(re_exp,1))
			hdu_disk.writeto(name_disk+'.fits', overwrite=True)
			print(str(ellip)+'-'+str(np.around(re_exp,1))+'exp.fits export done')
# catalog = np.zeros([int((re_devau_size[0]-1)*(re_exp_size[0]-1)*(1-0)/(re_devau_size[1]*re_exp_size[1]*0.1)),5]) # in order of ellip,re_devau,re_exp,B/T,re_total
# catalog_count = 0
# bulge_file_list = os.listdir(export_path+'bulge')
# bulge_file = []
# for file in bulge_file_list:
# 	if file.endswith('.fits'):
# 		bulge_file.append(file)
# disk_file_list = os.listdir(export_path+'disk')
# disk_file = []
# for file in disk_file_list:
# 	if file.endswith('.fits'):
# 		disk_file.append(file)
# for BT_ratio in np.arange(0.001,1.01,0.1):
# 	for ellip in np.arange(0.4,1,ellip_step):
# 		temp_devau = []
# 		temp_exp = []
# 		for file_devau in bulge_file:
# 			ellip_devau = file_devau[0:file_devau.index('-')]
# 			if ellip_devau == str(np.around(ellip,1)):
# 				temp_devau.append(file_devau)
# 		for file_exp in disk_file:
# 			ellip_exp = file_exp[0:file_exp.index('-')]
# 			if ellip_exp == ellip_devau:
# 				temp_exp.append(file_exp)
# 		for file_devau in temp_devau:
# 			re_devau = int(file_devau[file_devau.index('-')+1:file_devau.index('.fits')])
# 			image_devau = fits.open(export_path+'bulge/'+file_devau)[0].data
# 			image_devau = image_devau*BT_ratio
# 			for file_exp in temp_exp:
# 				re_exp = int(file_exp[file_exp.index('-')+1:file_exp.index('.fits')])
# 				image_exp = fits.open(export_path+'disk/'+file_exp)[0].data
# 				image_exp = image_exp*(1-BT_ratio)
# 				image_galaxy = image_devau+image_exp
# 				hdu = fits.PrimaryHDU(image_galaxy)
# 				name = export_path+str(np.around(BT_ratio,2))+'-'+str(np.around(ellip,1))+'-'+str(re_devau)+'-'+str(re_exp)
# 				hdu.writeto(name+'.fits', overwrite=True)
# 				print(str(np.around(BT_ratio,2))+'-'+str(np.around(ellip,1))+'-'+str(re_devau)+'-'+str(re_exp)+' image is done')
# 				re_total = effective_radius_total(image_galaxy,int((image_galaxy.shape[0]-1)/2),int((image_galaxy.shape[1]-1)/2))
# 				print('re_total = '+str(np.around(re_total,2)))
# 				catalog[catalog_count] = [BT_ratio,ellip,re_devau,re_exp,re_total]
# 				catalog_count += 1
# col1 = fits.Column(name='BT_ratio', array=catalog[:,0], format='D')
# col2 = fits.Column(name='ellipticity',array=catalog[:,1], format='D')
# col3 = fits.Column(name='re_devau',array=catalog[:,2], format='D')
# col4 = fits.Column(name='re_exp',array=catalog[:,3], format='D')
# col5 = fits.Column(name='re_total',array=catalog[:,4], format='D')
# table = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5])
# table.writeto(export_path+'BT_ratio.fits')
# print('CATALOG done')
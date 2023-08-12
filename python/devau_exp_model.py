import numpy as np
from astropy.modeling.models import Sersic2D
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.io import fits
# from photutils.datasets import make_noise_image # if need noise
from scipy.integrate import quad
from scipy.special import gammaincinv

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

def effective_radius_total(re_devau,re_exp,Ire,brightness,bt_ratio,hdu,x0,y0):
	# Calculate the effective radius of input image array.(in pixels)
	# input:
	# 	hdu: cut center pixel into subpixel of shape, elements in shape should be odd
	# 	x0,y0: galaxy center(in pixels)
	flux_total = brightness
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
# export name in order of devau-re,devau-ellip,exp-re,exp-ellip
shape = (2001,2001) # image size 2001x2001 pixels
# noise_mean=5.0 # with no noise
x,y = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
# ----------------------------------------------------------------------
# --------------------------- galaxy profile ---------------------------
# ----------------------------------------------------------------------
RE_surbri = 10 # surface brightness in Re
brightness_show = 100 # set brightness to more visible
ellip_step = 0.2
BT_ratio_step = 0.05
re_devau_size = [41,1]
re_exp_size = [51,1]
catalog = np.zeros([int((re_devau_size[0]-1)*(re_exp_size[0]-1)*(1-0)/(re_devau_size[1]*re_exp_size[1]*0.1)),5]) # in order of ellip,re_devau,re_exp,B/T,re_total
# catalog = np.zeros(5)
catalog_count = 0
for BT_ratio in np.arange(0.001,1.01,0.05):
	for ellip in np.arange(0,1,ellip_step):
		for re_devau in np.arange(5,re_devau_size[0],re_devau_size[1]): # --------------------------- de vaucouleurs profile ---------------------------
			model_devau = Sersic2D(amplitude = RE_surbri,r_eff = re_devau,n=4,x_0=int((shape[0]-1)/2),y_0=int((shape[1]-1)/2),ellip=ellip,theta=0)
			image_devau = model_devau(x,y)
			image_devau[int((shape[0]-1)/2),int((shape[1]-1)/2)] = modify_center_value([101,101],RE_surbri,re_devau,4)
			image_devau = image_devau/np.sum(image_devau) # normalization to total = 1
			hdu_bulge = fits.PrimaryHDU(image_devau)
			name_bulge = export_path+'bulge/'+str(ellip)+'-'+str(np.around(re_devau,1))
			hdu_bulge.writeto(name_bulge+'.fits', overwrite=True)
			for re_exp in np.arange(re_devau,re_exp_size[0],re_exp_size[1]): # --------------------------- exponential profile ---------------------------
				model_exp = Sersic2D(amplitude = RE_surbri,r_eff = re_exp,n=1,x_0=int((shape[0]-1)/2),y_0=int((shape[1]-1)/2),ellip=ellip,theta=0)
				image_exp = model_exp(x,y)
				image_exp[int((shape[0]-1)/2),int((shape[1]-1)/2)] = modify_center_value([101,101],RE_surbri,re_exp,1)
				image_exp = image_exp/np.sum(image_exp) # normalization to total = 1
				hdu_disk = fits.PrimaryHDU(image_exp)
				name_disk = export_path+'disk/'+str(ellip)+'-'+str(np.around(re_exp,1))
				hdu_disk.writeto(name_disk+'.fits', overwrite=True)
				image_galaxy = image_devau*brightness_show*BT_ratio + image_exp*brightness_show*(1-BT_ratio)
				re_total = effective_radius_total(re_devau,re_exp,RE_surbri,brightness_show,BT_ratio,image_galaxy,int((shape[0]-1)/2),int((shape[1]-1)/2))
				print('B/T: '+str(BT_ratio)+'\nellip: '+str(ellip)+'\nre_devau: '+str(re_devau)+'\nre_exp:'+str(re_exp)+'\nre_total:'+str(re_total))
				catalog[catalog_count] = [BT_ratio,np.around(re_devau,2),np.around(re_exp,1),np.around(re_total,1),ellip]
				catalog_count += 1
				hdu = fits.PrimaryHDU(image_galaxy)
				name = export_path+str(np.around(BT_ratio,2))+'-'+str(np.around(re_total,2))+'-'+str(np.around(re_devau,2))+'-'+str(np.around(re_exp,2))
				hdu.writeto(name+'.fits', overwrite=True)
col1 = fits.Column(name='BT_ratio', array=catalog[:,0], format='D')
col2 = fits.Column(name='re_devau',array=catalog[:,1], format='D')
col3 = fits.Column(name='re_exp',array=catalog[:,2], format='D')
col4 = fits.Column(name='re_total',array=catalog[:,3], format='D')
col5 = fits.Column(name='ellipticity',array=catalog[:,4], format='D')
table = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5])
table.writeto(export_path+'BT_ratio.fits')
print('CATALOG done')
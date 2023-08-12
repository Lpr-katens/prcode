# ===================================================
# Generate a catalog, which contains B/T, re_total, re_devau, re_exp
# ===================================================
import numpy as np
from astropy.modeling.models import Sersic2D
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.io import fits
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
	  '      **********************                 ***\n'
	  '        ******************               **********\n'
	  '          **************            ******************\n'
	  '            **********      ******************************\n'
	  '            *************************************************\n'
	  '          ****************************************************\n'
	  '        *******************************************************\n'
	  '      **********************************************************\n'
	  '      *******      *******                *******     *******\n'
	  '       *****        *****                  *****       *****\n'
	  '        ***          ***                    ***         ***\n'
	  '         *            *                      *           *\n')

def Intensity_r(r,re,Ire,n):
	bn = gammaincinv(2.*n,0.5)
	Intensity_r = Ire*np.exp(-bn*((r/re)**(1/n)-1))
	return Intensity_r

def integral_Intensity_r(r,re,Ire,n):
	bn = gammaincinv(2.*n,0.5)
	L_r = 2*np.pi*r*Ire*np.exp(-bn*((r/re)**(1/n)-1))
	return L_r

def re_from_equation(re_devau,re_exp,bt_ratio):
	Ire = 10
	luminosity_exp = sersic_luminosity(1,re_exp,Ire)
	luminosity_devau = sersic_luminosity(4,re_devau,Ire)
	for r in np.arange(1,max(re_exp,re_devau)+1,0.05):
		lr = quad(integral_Intensity_r,0,r,args=(re_exp,Ire,1))[0]*(1-bt_ratio)/luminosity_exp+quad(integral_Intensity_r,0,r,args=(re_devau,Ire,4))[0]*bt_ratio/luminosity_devau
		if lr>0.5:
			break
	return r

def sersic_profile(Ire,re,n,x0,y0,x,y):
	r = np.sqrt((x-x0)**2+(y-y0)**2)
	Ir = Intensity_r(r,re,Ire,n)
	return Ir

def sersic_luminosity(n,re,Ire):
	L_total = quad(integral_Intensity_r,0,float('inf'),args=(re,Ire,n))[0]
	return L_total

def modify_center_value(shape,Ire,re,n):
	temp = np.zeros(shape)
	for x in range(0,shape[0]):
		for y in range(0,shape[1]):
			temp[x,y] = sersic_profile(Ire,re,n,0.5,0.5,x/(shape[0]-1),y/(shape[1]-1))
	return np.mean(temp)

# def effective_radius_total(re_devau,re_exp,Ire,brightness,bt_ratio,hdu,x0,y0):
# 	# Calculate the effective radius of input image array.(in pixels)
# 	# input:
# 	# 	hdu: cut center pixel into subpixel of shape, elements in shape should be odd
# 	# 	x0,y0: galaxy center(in pixels)
# 	flux_total = brightness
# 	flux = 0
# 	for re in np.arange(1,np.sqrt(((hdu.shape[0]-1)/2)**2+((hdu.shape[1]-1)/2)**2),0.5):
# 		if flux >= flux_total/2:
# 			break
# 		for x in range(0,hdu.shape[0]):
# 			for y in range(0,hdu.shape[1]):
# 				if np.sqrt((x-x0)**2+(y-y0)**2)>re-0.5 and np.sqrt((x-x0)**2+(y-y0)**2)<=re:
# 					flux += hdu[x,y]
# 	return re-0.5
export_path = '/Users/lpr/Data/fits/expdata/BTratio/'
# export name in order of devau-re,devau-ellip,exp-re,exp-ellip
shape = (2001,2001) # image size 2001x2001 pixels
# # noise_mean=5.0 # with no noise
# x,y = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
# # ----------------------------------------------------------------------
# # --------------------------- galaxy profile ---------------------------
# # ----------------------------------------------------------------------
RE_surbri = 10 # surface brightness in Re
# brightness_show = 1000 # set brightness to more visible
re_devau_size = [51,2]
re_exp_size = [61,2]
catalog = np.zeros([int((re_devau_size[0]-1)*(re_exp_size[0]-1)*(1-0)/(re_devau_size[1]*re_exp_size[1]*0.05)),4]) # in order of ellip,re_devau,re_exp,B/T,re_total
# catalog = np.zeros(5)
catalog_count = 0
brightness = 100
for BT_ratio in np.arange(0.05,1.01,0.05):
	# for ellip in np.arange(0,1,ellip_step):
	# ellip=0 # not take ellipticity into account
	for re_devau in np.arange(1,5,re_devau_size[1]):
		for re_exp in np.arange(re_devau,re_exp_size[0],re_exp_size[1]):
	# for re_devau in np.arange(5,re_devau_size[0],re_devau_size[1]):
	# 	for re_exp in np.arange(re_devau,re_exp_size[0],re_exp_size[1]):
			re_total = re_from_equation(re_devau,re_exp,BT_ratio)
			catalog[catalog_count] = [BT_ratio,np.around(re_devau,4),np.around(re_exp,4),np.around(re_total,4)]
			catalog_count += 1
			print('B/T = '+str(np.around(BT_ratio,2))+'\nre_devau = '+str(np.around(re_devau,2))+'\nre_exp = '+str(np.around(re_exp,2)))
			# image = np.zeros(shape)
			# x_0,y_0 = 1000,1000
			# I_re = 10
			# luminosity_devau = sersic_luminosity(4,re_devau,I_re)
			# luminosity_exp = sersic_luminosity(1,re_exp,I_re)
			# for x in range(0,shape[0]):
			# 	for y in range(0,shape[1]):
			# 		image[x,y] = brightness*(sersic_profile(I_re,re_devau,4,x_0,y_0,x,y)*BT_ratio/luminosity_devau+sersic_profile(I_re,re_exp,1,x_0,y_0,x,y)*(1-BT_ratio)/luminosity_exp)
			# print('image center original value is: '+str(np.around(image[1000,1000],2)))
			# image[1000,1000] = brightness*(modify_center_value([100,100],I_re,re_devau,4)*BT_ratio/luminosity_devau+modify_center_value([100,100],I_re,re_exp,1)*(1-BT_ratio)/luminosity_exp)
			# print('image center modified value is: '+str(np.around(image[1000,1000],2)))
			# hdu = fits.PrimaryHDU(image)
			# name = export_path+'equation/'+str(np.around(BT_ratio,2))+'-'+str(np.around(re_devau,2))+'-'+str(np.around(re_exp,2))
			# hdu.writeto(name+'.fits', overwrite=True)
			# print('======== image is written ========')
col1 = fits.Column(name='BT_ratio', array=catalog[:,0], format='D')
col2 = fits.Column(name='re_devau',array=catalog[:,1], format='D')
col3 = fits.Column(name='re_exp',array=catalog[:,2], format='D')
col4 = fits.Column(name='re_total',array=catalog[:,3], format='D')
table = fits.BinTableHDU.from_columns([col1,col2,col3,col4])
table.writeto(export_path+'equation/B2T_equation_05_app.fits',overwrite=True)
print('CATALOG done')
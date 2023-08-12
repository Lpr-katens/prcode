
import numpy as np
from astropy.modeling.models import Sersic2D
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
	  '         *            *                              *           *)\n')

def Intensity_r(r,re,Ire,n):
	# sersic profile function
	# input:
	# 	Ire: surface brightness in re
	# 	re: effective radius
	# 	n: sersic index
	# 	r: radius
	bn = gammaincinv(2.*n,0.5)
	Intensity_r = Ire*np.exp(-bn*(r/re)**(1/n)-1)
	return Intensity_r

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
	# 	n: sersic index
	L_total = quad(Intensity_r,0,float('inf'),args=(re,Ire,n))[0]
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

def effective_radius_total(re_devau,re_exp,Ire,hdu,x0,y0):
	# Calculate the effective radius of input image array.(in pixels)
	# input:
	# 	hdu: cut center pixel into subpixel of shape, elements in shape should be odd
	# 	x0,y0: galaxy center(in pixels)
	flux_total = sersic_luminosity(1,re_exp,Ire)+sersic_luminosity(4,re_devau,Ire)
	print('re_exp = '+str(re_exp)+'\nre_devau = '+str(re_devau)+'\nIre = '+str(Ire))
	print(hdu[100,100])
	flux_re = 0
	flux = 0
	print(flux_total)
	for re in np.arange(1,np.sqrt(((hdu.shape[0]-1)/2)**2+((hdu.shape[1]-1)/2)**2),2):
		if flux >= flux_total/2:
			break
		for r in range(1,int(re)):
			for x in range(0,hdu.shape[0]):
				for y in range(0,hdu.shape[1]):
					if np.sqrt((x-x0)**2+(y-y0)**2) <= r:
						print(x,y,r)
						flux += hdu[x,y]
						print(flux)
	return re-2
export_path = '/Users/lpr/Data/fits/expdata/BTratio/'
# export name in order of devau-re,devau-ellip,exp-re,exp-ellip
shape = (201,201) # image size 200x200 pixels
# noise_mean=5.0 # with no noise
x,y = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
# ----------------------------------------------------------------------
# --------------------------- galaxy profile ---------------------------
# ----------------------------------------------------------------------
RE_surbri = 10 # surface brightness in Re
re_devau_size = [51,2]
ellip_step = 0.2
re_exp_size = [61,2]
catalog = np.zeros([int((re_devau_size[0]-1)*(re_exp_size[0]-1)*(1-0)/(re_devau_size[1]*re_exp_size[1]*0.1)),5]) # in order of ellip,re_devau,re_exp,B/T,re_total
# catalog = np.zeros(5)
catalog_count = 0
for ellip in np.arange(0,1,ellip_step):
	for re_devau in np.arange(5,re_devau_size[0],re_devau_size[1]): # --------------------------- de vaucouleurs profile ---------------------------
		model_devau = Sersic2D(amplitude = RE_surbri,r_eff = re_devau,n=4,x_0=100,y_0=100,ellip=ellip,theta=0)
		image_devau = model_devau(x,y)
		image_devau[int((shape[0]-1)/2),int((shape[1]-1)/2)] = modify_center_value([101,101],RE_surbri,re_devau,4)
		mag_b = 25.96-2.5*np.log10(sersic_luminosity(4,re_devau,RE_surbri))
		for re_exp in np.arange(re_devau,re_exp_size[0],re_exp_size[1]): # --------------------------- exponential profile ---------------------------
			model_exp = Sersic2D(amplitude = RE_surbri,r_eff = re_exp,n=1,x_0=100,y_0=100,ellip=ellip,theta=0)
			image_exp = model_exp(x,y)
			mag_d = 25.96-2.5*np.log10(sersic_luminosity(1,re_exp,RE_surbri))
			DB_ratio = 10**((mag_b-mag_d)/2.5)
			BT_ratio = 1/(DB_ratio+1)
			image_galaxy = image_devau + image_exp
			re_total = effective_radius_total(re_devau,re_exp,RE_surbri,image_galaxy,100,100)
			catalog[catalog_count] = [ellip,re_devau,re_exp,BT_ratio,re_total]
			catalog_count += 1
			# hdu = fits.PrimaryHDU(image_galaxy)
			# name = export_path+str(ellip)+'-'+str(re_devau)+'-'+str(re_exp)
			# hdu.writeto(name+'.fits', overwrite=True)
			# log_image = np.log10(image_galaxy)
			# print('------------------ imshow ------------------')
			# plt.figure(figsize=(20,20))
			# plt.imshow(log_image,cmap='gray',interpolation=None,origin='lower')
			# plt.colorbar()
			# plt.savefig(name+'.eps')
			# plt.close()
			print('B/T: '+str(BT_ratio)+'\nellip: '+str(ellip)+'\nre_devau: '+str(re_devau)+'\nre_exp:'+str(re_exp)+'\nre_total:'+str(re_total))
# ######================================================================================
# ######|||		The first version is color gradient with multiple annuli.
# ######|||		The data cut off at annulus with s/n<10.
# ######================================================================================
# import matplotlib.pyplot as plt
# from lpr.image import photometry as pm
# from astropy.io import fits
# import numpy as np
# import os
# from scipy.integrate import quad
# import matplotlib.gridspec as gridspec
# from scipy.optimize import curve_fit
# from astropy.cosmology import FlatLambdaCDM as flcdm
# import matplotlib.colors as colors
# from matplotlib.patches import Ellipse, Circle
# def kpc_per_arcsec(z):
#     angular_distance = flcdm(H0=70,Om0=0.3).angular_diameter_distance(z).to_value() # to_value transform the value with unit quantity to the value with unit float
#     arc_scale = angular_distance*np.pi*1000/(180*3600)
#     return(arc_scale)
# def rescale(arr):
#     arr_min = arr.min()
#     arr_max = arr.max()
#     return (arr - arr_min) / (arr_max - arr_min)
# def gaussian(x,mu,sigma):
#     return np.exp(-(x-mu)**2/(2*sigma**2))/(sigma*np.sqrt(2*np.pi))

# short_wave_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/convolveGALAXYimage_gaussian/goodsn_f606w/'
# long_wave_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/convolveGALAXYimage_gaussian/goodsn_f125w/'
# expath = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/color_gradient/'
# cata_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/goodsn_Huangall_radecmatch_modifyz_kpc.fits'
# mask_path = '/Users/lpr/Data/SExtractor/goodsn_f125w/'
# img_list = []
# for file in os.listdir(long_wave_path):
# 	if file.endswith('.fits'):
# 		img_list.append(file)
# cata = fits.open(cata_path)[1].data
# cata = cata[np.where((cata['z_used']>0.8)&(cata['z_used']<1.3))]
# sed_type = ['AGN','Composite','Star-Forming','Quiescent','Blue-compact']
# for num1 in range(0,len(cata)):
# 	idx = cata[num1]['ID']
# 	print(idx)
# 	if 'goodsn_'+str(idx)+'_f125w.fits' in img_list:
# 		z = cata[num1]['z_used']
# 		lx = np.around(cata[num1]['LX'],2)
# 		l45 = np.around(np.log10(cata[num1]['L4P5EX']),2)
# 		sed = int(cata[num1]['TMP_CLASS'][0])-1
# 		print(sed)
# 		re = np.around(cata[num1]['kpc'],2)
# 		n = np.around(cata[num1]['n_f160w'],1)
# 		mass = np.around(cata[num1]['MASS'],2)
# 		J_band_img = fits.open(long_wave_path+'goodsn_'+str(idx)+'_f125w.fits')[0]
# 		J_band_imgdata = J_band_img.data
# 		histo_ydata_J,t_J = np.histogram(J_band_imgdata.flatten(),bins=100)
# 		histo_xdata_J = t_J[0:-1]
# 		try:
# 			popt_J,pcov_J = curve_fit(gaussian,histo_xdata_J,histo_ydata_J)
# 		except RuntimeError:
# 			popt_J = [np.median(np.sort(J_band_imgdata.flatten())[0:100]),np.std(J_band_imgdata)]
# 		J_band_seg_img = fits.open(mask_path+str(idx)+'_seg.fits')[0].data
# 		J_band_seg_cat = fits.open(mask_path+str(idx)+'_cat.fits')[1].data
# 		J_distance = np.sqrt((J_band_seg_cat['X_IMAGE']-50)**2 + (J_band_seg_cat['Y_IMAGE']-50)**2)
# 		ind = np.where(J_distance==np.min(J_distance))[0]
# 		J_background = np.random.normal(loc=popt_J[0],scale=popt_J[1],size=J_band_imgdata.shape)
# 		J_band_imgdata[np.where((J_band_seg_img!=ind+1)&(J_band_seg_img!=0))] = J_background[np.where((J_band_seg_img!=ind+1)&(J_band_seg_img!=0))]
# 		J_x,J_y,J_yerr,J_aper = pm.aperture_flux(J_band_imgdata,ccdgain=J_band_img.header['CCDGAIN'],cx=50,cy=50,bins=16,poisson=False)
# 		J_ZP = -2.5*np.log10(J_band_img.header['PHOTFLAM'])-21.1-5*np.log10(J_band_img.header['PHOTPLAM'])+18.692
# 		J_yerr = -2.5*np.array(J_yerr)/(np.array(J_y)*np.log(10))
# 		J_y = J_ZP-2.5*np.log10(J_y)
# 		V_band_img = fits.open(short_wave_path+'goodsn_'+str(idx)+'_f606w.fits')[0]
# 		V_band_imgdata = V_band_img.data
# 		histo_ydata_V,t_V = np.histogram(V_band_imgdata.flatten(),bins=100)
# 		histo_xdata_V = t_V[0:-1]
# 		try:
# 			popt_V,pcov_V = curve_fit(gaussian,histo_xdata_V,histo_ydata_V)
# 		except RuntimeError:
# 			popt_V = [np.median(np.sort(V_band_imgdata.flatten())[0:100]),np.std(V_band_imgdata)]
# 		V_background = np.random.normal(loc=popt_V[0],scale=popt_V[1],size=V_band_imgdata.shape)
# 		V_band_imgdata[np.where((J_band_seg_img!=ind+1)&(J_band_seg_img!=0))] = V_background[np.where((J_band_seg_img!=ind+1)&(J_band_seg_img!=0))]
# 		V_x,V_y,V_yerr,V_aper = pm.aperture_flux(V_band_imgdata,ccdgain=V_band_img.header['CCDGAIN'],cx=50,cy=50,bins=16,poisson=False)
# 		V_ZP = -2.5*np.log10(V_band_img.header['PHOTFLAM'])-21.1-5*np.log10(V_band_img.header['PHOTPLAM'])+18.692
# 		V_yerr = -2.5*np.array(V_yerr)/(np.array(V_y)*np.log(10))
# 		V_y = V_ZP-2.5*np.log10(V_y)
# 		J_x = (np.array(J_x)*J_band_img.header['CD2_2']*3600)*kpc_per_arcsec(z)
# 		color = V_y-J_y
# 		colorerr = np.sqrt(V_yerr**2+J_yerr**2)
# 		sn = color/colorerr
# 		sn_idx = np.where(abs(sn)>10)
# 		J_x = J_x[sn_idx]
# 		color = color[sn_idx]
# 		colorerr = colorerr[sn_idx]
# 		textstr = '\n'.join((
#     		r'$L_{X}=$'+str(lx),
#     		r'$L_{4.5}^{EXC}=$'+str(l45),
#     		r'$r_e=$'+str(re),
#     		r'MASS='+str(mass),
#     		r'n='+str(n)))
# 		if len(color) != 0:
# 			fig = plt.figure(figsize=[10,12])
# 			gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[2, 3])
# 			ax = fig.add_subplot(gs[1,:])
# 			ax.set_title(str(idx)+' color gradient')
# 			ax.errorbar(J_x,color,yerr=colorerr,color='blue',ecolor='r')
# 			ax.set_ylim(np.nanmin(color)-1,np.nanmax(color)+1)
# 			ax.set_ylabel('V-J')
# 			ax.set_xlabel('radius (kpc)')
# 			ax = fig.add_subplot(gs[0,0])
# 			ax.imshow(rescale(J_band_imgdata),cmap='Greys',norm=colors.PowerNorm(gamma=0.3,vmin=np.min(rescale(J_band_imgdata)),vmax=np.max(rescale(J_band_imgdata))))
# 			ax.set_title('J band image')
# 			for num2 in range(0,len(J_aper)):
# 				J_aper[num2].plot(color='red',lw=0.1,alpha=0.1)
# 			circle_theta = np.arange(0,2*np.pi,0.01)
# 			circle_x = 50+(10/kpc_per_arcsec(z))/0.06*np.cos(circle_theta)
# 			circle_y = 50+(10/kpc_per_arcsec(z))/0.06*np.sin(circle_theta)
# 			ax.plot(circle_x,circle_y,color='blue')
# 			ax = fig.add_subplot(gs[0,1])
# 			ax.imshow(rescale(V_band_imgdata),cmap='Greys',norm=colors.PowerNorm(gamma=0.4,vmin=np.min(rescale(V_band_imgdata)),vmax=np.max(rescale(V_band_imgdata))))
# 			ax.set_title('V band image')
# 			for num2 in range(0,len(J_aper)):
# 				J_aper[num2].plot(color='red',lw=0.1,alpha=0.1)
# 			ax.plot(circle_x,circle_y,color='blue')
# 			fig.text(0.02,0.89,textstr,fontsize=15,color='red')
# 			fig.text(0.02,0.02,sed_type[sed],fontsize=15,color='red')
# 			fig.text(0.02,0.04,'z='+str(np.around(z,2)),fontsize=15,color='red')
# 			plt.savefig(expath+str(idx)+'_V-J.pdf')
# 			plt.close()
# 		print('=========================== '+str(idx)+' is done ===========================')
# print('=========================== GOODSN DONE ===========================')
# #####================================================================================
# #####|||		The second version is color gradient with two annuli.
# #####|||		The bulge size I use is come from Peletier et. al.(1996), which is 0.5 R_eff.
# #####|||		The disk size I use is annulus with s/n<10.
# #####================================================================================
import matplotlib.pyplot as plt
from lpr.image import photometry as pm
from astropy.io import fits
import numpy as np
import os
from scipy.integrate import quad
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from astropy.cosmology import FlatLambdaCDM as flcdm
import matplotlib.colors as colors
from matplotlib.patches import Ellipse, Circle
import numpy as np
from scipy.optimize import curve_fit
from photutils.utils import calc_total_error as photoerr
from photutils.aperture import CircularAnnulus as can
from photutils.aperture import CircularAperture as cap
from photutils.aperture import aperture_photometry as ap
def kpc_per_arcsec(z):
    angular_distance = flcdm(H0=70,Om0=0.3).angular_diameter_distance(z).to_value() # to_value transform the value with unit quantity to the value with unit float
    arc_scale = angular_distance*np.pi*1000/(180*3600)
    return(arc_scale)
def rescale(arr):
    arr_min = arr.min()
    arr_max = arr.max()
    return (arr - arr_min) / (arr_max - arr_min)
def gaussian(x,mu,sigma):
    return np.exp(-(x-mu)**2/(2*sigma**2))/(sigma*np.sqrt(2*np.pi))
def r_determine(img,cx,cy,totalerr):
	for r in range(1,int(img.shape[0]/2)-1):
		aperture = can([cx,cy],r,r+1)
		aper_flux = ap(img,aperture,totalerr,method='exact')
		if aper_flux['aperture_sum'][0]/aper_flux['aperture_sum_err'][0]<10:
			break;
	return r

short_wave_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/convolveGALAXYimage_gaussian/goodsn_f606w/'
long_wave_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/convolveGALAXYimage_gaussian/goodsn_f125w/'
expath = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/color_difference/'
cata_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/goodsn_Huangall_radecmatch_modifyz_kpc.fits'
mask_path = '/Users/lpr/Data/SExtractor/goodsn_f125w/'
wht_606_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/singleGALAXYimage_wht/goodsn_f606w/'
wht_125_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/singleGALAXYimage_wht/goodsn_f125w/'
img_list = []
for file in os.listdir(long_wave_path):
	if file.endswith('.fits'):
		img_list.append(file)
cata = fits.open(cata_path)[1].data
cata = cata[np.where((cata['z_used']>0.8)&(cata['z_used']<1.3)&(cata['re_f160w']>0))]
sed_type = ['AGN','Composite','Star-Forming','Quiescent','Blue-compact']
for num1 in range(0,len(cata)):
	idx = cata[num1]['ID']
	print(idx)
	if 'goodsn_'+str(idx)+'_f125w.fits' in img_list:
		J_flux = []
		J_error = []
		J_aperture_list = []
		V_flux = []
		V_error = []
		V_aperture_list = []
		z = cata[num1]['z_used']
		lx = np.around(cata[num1]['LX'],2)
		l45 = np.around(np.log10(cata[num1]['L4P5EX']),2)
		sed = int(cata[num1]['TMP_CLASS'][0])-1
		re = np.around(cata[num1]['kpc'],2)
		n = np.around(cata[num1]['n_f160w'],1)
		mass = np.around(cata[num1]['MASS'],2)
		# --------------------------------------------------
		J_band_img = fits.open(long_wave_path+'goodsn_'+str(idx)+'_f125w.fits')[0]
		J_band_imgdata = J_band_img.data
		V_band_img = fits.open(short_wave_path+'goodsn_'+str(idx)+'_f606w.fits')[0]
		V_band_imgdata = V_band_img.data
		J_band_whtimg = fits.open(wht_125_path+'goodsn_'+str(idx)+'_f125w_wht.fits')[0].data
		V_band_whtimg = fits.open(wht_606_path+'goodsn_'+str(idx)+'_f606w_wht.fits')[0].data
		J_totalerr = photoerr(J_band_imgdata,np.sqrt(1/np.median(J_band_whtimg.flatten())),J_band_img.header['EXPTIME'])
		V_totalerr = photoerr(V_band_imgdata,np.sqrt(1/np.median(V_band_whtimg.flatten())),V_band_img.header['EXPTIME'])
		r_bulge = cata[num1]['re_f160w']*0.5/(J_band_img.header['CD2_2']*3600)
		# r_disk_J,r_disk_V = r_determine(J_band_img.data,cx=50,cy=50,totalerr=J_totalerr),r_determine(V_band_img.data,cx=50,cy=50,totalerr=V_totalerr)
		# r_disk = min(r_disk_J,r_disk_V)
		r_disk = r_determine(J_band_img.data,cx=50,cy=50,totalerr=J_totalerr)
		if r_bulge < r_disk:
			histo_ydata_J,t_J = np.histogram(J_band_imgdata.flatten(),bins=100)
			histo_xdata_J = t_J[0:-1]
			try:
				popt_J,pcov_J = curve_fit(gaussian,histo_xdata_J,histo_ydata_J)
			except RuntimeError:
				popt_J = [np.median(np.sort(J_band_imgdata.flatten())[0:100]),np.std(J_band_imgdata)]
			J_band_seg_img = fits.open(mask_path+str(idx)+'_seg.fits')[0].data
			J_band_seg_cat = fits.open(mask_path+str(idx)+'_cat.fits')[1].data
			J_distance = np.sqrt((J_band_seg_cat['X_IMAGE']-50)**2 + (J_band_seg_cat['Y_IMAGE']-50)**2)
			ind = np.where(J_distance==np.min(J_distance))[0]
			J_background = np.random.normal(loc=popt_J[0],scale=popt_J[1],size=J_band_imgdata.shape)
			J_band_imgdata[np.where((J_band_seg_img!=ind+1)&(J_band_seg_img!=0))] = J_background[np.where((J_band_seg_img!=ind+1)&(J_band_seg_img!=0))]
			aperture = cap([50,50],r_bulge)
			aper_flux = ap(J_band_imgdata,aperture,J_totalerr,method='exact')
			J_flux.append(aper_flux['aperture_sum'][0])
			J_error.append(aper_flux['aperture_sum_err'][0])
			J_aperture_list.append(aperture)
			aperture = can([50,50],r_bulge,r_disk)
			aper_flux = ap(J_band_imgdata,aperture,J_totalerr,method='exact')
			J_flux.append(aper_flux['aperture_sum'][0])
			J_error.append(aper_flux['aperture_sum_err'][0])
			J_aperture_list.append(aperture)
			J_ZP = -2.5*np.log10(J_band_img.header['PHOTFLAM'])-21.1-5*np.log10(J_band_img.header['PHOTPLAM'])+18.692
			J_yerr = abs(-2.5*np.array(J_error)/(np.array(J_flux)*np.log(10)))
			J_flux = J_ZP-2.5*np.log10(J_flux)
			histo_ydata_V,t_V = np.histogram(V_band_imgdata.flatten(),bins=100)
			histo_xdata_V = t_V[0:-1]
			try:
				popt_V,pcov_V = curve_fit(gaussian,histo_xdata_V,histo_ydata_V)
			except RuntimeError:
				popt_V = [np.median(np.sort(V_band_imgdata.flatten())[0:100]),np.std(V_band_imgdata)]
			V_background = np.random.normal(loc=popt_V[0],scale=popt_V[1],size=V_band_imgdata.shape)
			V_band_imgdata[np.where((J_band_seg_img!=ind+1)&(J_band_seg_img!=0))] = V_background[np.where((J_band_seg_img!=ind+1)&(J_band_seg_img!=0))]
			aperture = cap([50,50],r_bulge)
			aper_flux = ap(V_band_imgdata,aperture,V_totalerr,method='exact')
			V_flux.append(aper_flux['aperture_sum'][0])
			V_error.append(aper_flux['aperture_sum_err'][0])
			V_aperture_list.append(aperture)
			aperture = can([50,50],r_bulge,r_disk)
			aper_flux = ap(V_band_imgdata,aperture,V_totalerr,method='exact')
			V_flux.append(aper_flux['aperture_sum'][0])
			V_error.append(aper_flux['aperture_sum_err'][0])
			V_aperture_list.append(aperture)
			V_ZP = -2.5*np.log10(V_band_img.header['PHOTFLAM'])-21.1-5*np.log10(V_band_img.header['PHOTPLAM'])+18.692
			V_yerr = abs(-2.5*np.array(V_error)/(np.array(V_flux)*np.log(10)))
			V_flux = V_ZP-2.5*np.log10(V_flux)
			x_axis = [r_bulge,r_disk]
			x_axis = (np.array(x_axis)*J_band_img.header['CD2_2']*3600)*kpc_per_arcsec(z)
			if len(J_flux) > 0 and len(V_flux) > 0:
				color = V_flux-J_flux
				colorerr = np.sqrt(V_yerr**2+J_yerr**2)
			else:
				color = []
			textstr = '\n'.join((
    			r'$L_{X}=$'+str(lx),
    			r'$L_{4.5}^{EXC}=$'+str(l45),
    			r'$r_e=$'+str(re),
    			r'MASS='+str(mass),
    			r'n='+str(n)))
			if len(color) != 0:
				fig = plt.figure(figsize=[10,12])
				gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[2, 3])
				ax = fig.add_subplot(gs[1,:])
				ax.set_title(str(idx)+' color difference')
				ax.errorbar(x_axis,color,yerr=colorerr,color='blue',ecolor='r')
				ax.set_ylim(np.nanmin(color)-0.1,np.nanmax(color)+0.1)
				ax.set_ylabel('F606W-F125W')
				ax.set_xlabel('radius (kpc)')
				ax = fig.add_subplot(gs[0,0])
				ax.imshow(rescale(J_band_imgdata),cmap='Greys',norm=colors.PowerNorm(gamma=0.3,vmin=np.min(rescale(J_band_imgdata)),vmax=np.max(rescale(J_band_imgdata))))
				ax.set_title('J band image')
				for num2 in range(0,len(J_aperture_list)):
					J_aperture_list[num2].plot(color='red',lw=0.3,alpha=1)
				circle_theta = np.arange(0,2*np.pi,0.01)
				circle_x = 50+(10/kpc_per_arcsec(z))/0.06*np.cos(circle_theta)
				circle_y = 50+(10/kpc_per_arcsec(z))/0.06*np.sin(circle_theta)
				ax.plot(circle_x,circle_y,color='blue')
				ax = fig.add_subplot(gs[0,1])
				ax.imshow(rescale(V_band_imgdata),cmap='Greys',norm=colors.PowerNorm(gamma=0.4,vmin=np.min(rescale(V_band_imgdata)),vmax=np.max(rescale(V_band_imgdata))))
				ax.set_title('V band image')
				for num2 in range(0,len(J_aperture_list)):
					J_aperture_list[num2].plot(color='red',lw=0.3,alpha=1)
				ax.plot(circle_x,circle_y,color='blue')
				fig.text(0.02,0.89,textstr,fontsize=15,color='red')
				fig.text(0.02,0.02,sed_type[sed],fontsize=15,color='red')
				fig.text(0.02,0.04,'z='+str(np.around(z,2)),fontsize=15,color='red')
				plt.savefig(expath+str(idx)+'_V-J.pdf')
				plt.close()
			print('=========================== '+str(idx)+' is done ===========================')
print('=========================== GOODSN DONE ===========================')
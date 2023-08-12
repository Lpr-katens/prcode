import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import os
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit

##===##===##===function beginning##===##===##===
# center determination from catalog
##===##===##===##===##===##===##===##===##===##===
# define a function to determine the center of every galaxy from a exist catalog
def catalog_center(image_direc_path,cata_path):
    catalog = fits.open(cata_path)[1].data # catalog contains coordinates and redshift
    col = catalog.columns
    z_id,ra_id,dec_id,flux_id,re160_id,re125_id,n160_id,n125_id = 0,0,0,0,0,0,0,0
    for temp in range(0,len(col)):
            if col[temp].name == 'zbest':
                z_id = temp
            elif col[temp].name == 'RA_1':
                ra_id = temp
            elif col[temp].name == 'DEC_1':
                dec_id = temp
            elif col[temp].name == 'FLUX_APER_1_F160W':
                flux_id = temp
            elif col[temp].name == 're_f160w':
                re160_id = temp
            elif col[temp].name == 're_f125w':
                re125_id = temp
            elif col[temp].name == 'n_f160w':
                n160_id = temp
            elif col[temp].name == 'n_f125w':
                n125_id = temp
    filelist = os.listdir(image_direc_path)
    file = []
    for name in filelist:
        if name.endswith('.fits'):
            file.append(name)
    coordlist = np.zeros([len(file),12])
    for num in range(0,len(file)):
        # get the WCS  FK5 of image.maxpixel
        # index2: ra, index3:dec
        gala_idx = file[num][file[num].index('_',11)+1:file[num].index('.fits')]
        [ra,dec] = galaxy_center(image_direc_path + '/' + file[num])
        ra_cata,dec_cata,z_cata,flux_cata,re160,re125,n160,n125 = 0,0,0,0,0,0,0,0
        for num2 in range(0,len(catalog)):
            if num2 == 0:
                distance = np.sqrt((catalog[num2][ra_id]-ra)**2 + (catalog[num2][dec_id]-dec)**2)
                ra_cata = catalog[num2][ra_id]
                dec_cata = catalog[num2][dec_id]
                z_cata = catalog[num2][z_id]
                flux_cata = catalog[num2][flux_id]
                re160 = catalog[num2][re160_id]
                re125 = catalog[num2][re125_id]
                n160 = catalog[num2][n160_id]
                n125 = catalog[num2][n125_id]
            elif np.sqrt((catalog[num2][ra_id]-ra)**2 + (catalog[num2][dec_id]-dec)**2) <= distance:
                distance = np.sqrt((catalog[num2][ra_id]-ra)**2 + (catalog[num2][dec_id]-dec)**2)
                ra_cata = catalog[num2][ra_id]
                dec_cata = catalog[num2][dec_id]
                z_cata = catalog[num2][z_id]
                flux_cata = catalog[num2][flux_id]
                re160 = catalog[num2][re160_id]
                re125 = catalog[num2][re125_id]
                n160 = catalog[num2][n160_id]
                n125 = catalog[num2][n125_id]
        coordlist[num] = int(gala_idx),ra,dec,distance,ra_cata,dec_cata,z_cata,flux_cata,re160,re125,n160,n125
    return coordlist
##===##===##===function ending##===##===##===

##===##===##===function beginning##===##===##===
# center determination of image
##===##===##===##===##===##===##===##===##===##===
# define a function to determine the center of every galaxy
def galaxy_center(image):
    hdu_temp = fits.open(image)[0]
    hdu_temp.header['CTYPE1'] = 'RA---TAN-SIP'
    hdu_temp.header['CTYPE2'] = 'DEC--TAN-SIP'
    his = hdu_temp.header['HISTORY'][1]
    cut_range = np.array([int(his[his.index('[')+2:his.index(':',16)]),int(his[his.index(',')+3:his.index(':',his.index(','))])])
    x,y = 50+cut_range[0],50+cut_range[1]
    [ra,dec] = WCS(hdu_temp.header).wcs_pix2world(x,y,0)
    return [ra,dec]
# def galaxy_center(image):
#     hdu_temp = fits.open(image)[0]
#     hdu = np.transpose(hdu_temp.data) # astropy transpose x and y axis of fits
#     rho_x = np.zeros(hdu.shape[0]) # marginal sum of x axis, namely, sum along j
#     rho_y = np.zeros(hdu.shape[1]) # marginal sum of y axis, namely, sum along i
#     for x in range(0,hdu.shape[0]):
#         sum = 0
#         for y in range(0,hdu.shape[1]):
#             sum += hdu[x,y]
#         rho_x[x] = sum
#     for y in range(0,hdu.shape[1]):
#         sum = 0
#         for x in range(0,hdu.shape[0]):
#             sum += hdu[x,y]
#         rho_y[y] = sum
#     meaninten_x = np.sum(rho_x)/hdu.shape[0]
#     meaninten_y = np.sum(rho_y)/hdu.shape[1]
#     #============
#     # find new center coor_x
#     for x in range(0,hdu.shape[0]):
#         if rho_x[x] < meaninten_x:
#             rho_x[x] = 0
#         else:
#             rho_x[x] = rho_x[x] - meaninten_x
#     numerator = 0
#     for x in range(0,len(rho_x)):
#         numerator += rho_x[x]*(x+1) # make image begin with index 1
#     coor_x = numerator/np.sum(rho_x)
#     #============
#     # find new center coor_y
#     for y in range(0,hdu.shape[1]):
#         if rho_y[y] < meaninten_y:
#             rho_y[y] = 0
#         else:
#             rho_y[y] = rho_y[y] - meaninten_y
#     numerator = 0
#     for y in range(0,len(rho_y)):
#         numerator += rho_y[y]*(y+1) # make image begin with indey 1
#     coor_y = numerator/np.sum(rho_y)
#     #============
#     # convert new x_image and y_image into ra,dec
#     hdu_temp.header['CTYPE1'] = 'RA---TAN-SIP'
#     hdu_temp.header['CTYPE2'] = 'DEC--TAN-SIP'
#     [ra,dec] = WCS(hdu_temp.header).wcs_pix2world(coor_x,coor_y,0)
#     return [ra,dec]
##================================================

# ##===##===##===function beginning##===##===##===
# # fwhm determination of image
# ##===##===##===##===##===##===##===##===##===##===
# # define a function to determine the fwhm of every galaxy
# def fwhm(image):
#     hdu = np.transpose(image)
#     y = hdu[int(hdu.shape[0]/2)]
#     y_max = np.max(y)
#     x_index = []
#     for num in range(1,len(y)):
#         if y[num-1] < y_max/2 and y_max/2 <= y[num]:
#             x_index.append((2*num+1)/2)
#         elif y[num-1] >= y_max/2 and y_max/2 > y[num]:
#             x_index.append((2*num+1)/2)
#     if len(x_index) > 1:
#         fwhm = x_index[1]-x_index[0]
#     elif len(x_index) <= 1:
#         fwhm = 1
#     return fwhm
# ##================================================

# ##===##===##===function beginning##===##===##===
# #===============================
# # gaussian function
# #===============================
# def gaussian(x,a,mu,sigma):
#     return a*np.exp(-((x-mu)/sigma)**2)
# ###############################===============================#######################

##===##===##===function beginning##===##===##===
#===============================
# subtract background noise
#===============================
# define a function to subtract background noise from image
# def subt_backnoise(image):
#     hdu = np.transpose(image.data)
#     hdu1 = hdu[np.where((hdu>-0.01)&(hdu<0.01))]
#     n,bins,patches = plt.hist(hdu1,bins=1000,range=[-0.01,0.01])
#     try:
#         popt,pcov = curve_fit(gaussian,bins[1:],n)
#     except RuntimeError:
#         hdu2 = hdu
#         return(hdu2,0)
#     else:
#         plt.plot(bins,gaussian(bins,popt[0],popt[1],popt[2]))
#         hdu2 = hdu-popt[1]
#         return(hdu2,popt[1])
def subt_backnoise(image):
    hdu = np.transpose(image.data)
    hdu1 = hdu[np.where((hdu>-0.01)&(hdu<0.01))]
    noise = np.median(hdu1)
    return(hdu-noise,noise)
###############################===============================#######################

image_folder_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/multiband_cutimage'
filelist = os.listdir(image_folder_path)
dirlist = []
for directory in filelist:
    if directory.startswith('goodsn'):
        dirlist.append(directory)
        
# we use stdima to get image index, such as 597, 1092......
filelist = os.listdir(image_folder_path+'/goodsn_f098m')
stdima_order = []
for file in filelist:
    if file.endswith('.fits'):
        stdima_order.append(file)

image_path = image_folder_path+'/goodsn_f160w'
goodsn_catalog_path = '/Users/lpr/Data/fits/pridata/gdn_all.fits'
# match image location with candels catalog
center_of_catalog = catalog_center(image_path,goodsn_catalog_path)

# save the match table to fits
tablepath = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/'
col1 = fits.Column(name='galaxy_id', array=center_of_catalog[:,0], format='K')
col2 = fits.Column(name='image_ra',array=center_of_catalog[:,1], format='D')
col3 = fits.Column(name='image_dec',array=center_of_catalog[:,2], format='D')
col4 = fits.Column(name='distance_offset',array=center_of_catalog[:,3], format='D')
col5 = fits.Column(name='catalog_ra',array=center_of_catalog[:,4], format='D')
col6 = fits.Column(name='catalog_dec',array=center_of_catalog[:,5], format='D')
col7 = fits.Column(name='catalog_redshift',array=center_of_catalog[:,6], format='D')
col8 = fits.Column(name='catalog_flux_aper_1_f160w',array=center_of_catalog[:,7],format='D')
col9 = fits.Column(name='re_f160w',array=center_of_catalog[:,8],format='D')
col10 = fits.Column(name='re_f125w',array=center_of_catalog[:,9],format='D')
col11 = fits.Column(name='n_f160w',array=center_of_catalog[:,10],format='D')
col12 = fits.Column(name='n_f125w',array=center_of_catalog[:,11],format='D')
table = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12])
table.writeto(tablepath+'galaxy_match.fits')
print('done')
# # # draw SED diagram
# # # aperture set to 1.5"
# # flux_list = {}
# # for stdfile in stdima_order:
# #     image_index = stdfile[stdfile.index('_',7)+1:stdfile.index('.fits')]
# #     galaxy_name = 'goodsn_'+image_index+'_f160w.fits'
# #     image_160 = fits.open(image_path+'/'+galaxy_name)[0]
# #     elec_to_flux_160 = float(image_160.header['PHOTFNU'])*10**6
# #     image_160.header['CTYPE1'] = 'RA---TAN-SIP'
# #     image_160.header['CTYPE2'] = 'DEC--TAN-SIP'
# #     center_indices = WCS(image_160.header).wcs_world2pix(center_of_catalog[np.where(center_of_catalog == int(image_index))[0],1][0],center_of_catalog[np.where(center_of_catalog == int(image_index))[0],2][0],0)
# #     flux_160 = 0
# #     aftersubt_image_160,sky_160 = subt_backnoise(image_160.data, fwhm(image_160.data),center_indices[0],center_indices[1])
# #     # from last step, image array has been transposed
# #     pixel_number = 0
# #     for num_x in range(0,aftersubt_image_160.shape[0]):
# #         for num_y in range(0,aftersubt_image_160.shape[1]):
# #             if np.sqrt((num_x-float(center_indices[0]))**2+(num_y-float(center_indices[1]))**2)<=25:
# #                 flux_160 += aftersubt_image_160[num_x][num_y]
# #                 pixel_number += 1
# #     err_160 = np.sqrt(pixel_number)*sky_160
# #     flux_list[1600] = flux_160*elec_to_flux_160,err_160*elec_to_flux_160
# #     #============ append band in x axis of SED 
# #     x_axis = np.zeros(0) # band list
# #     for band in dirlist:
# #         if '1' not in band[band.index('f')+1:band.index('f')+2] and '0'not in band[band.index('f')+1:band.index('f')+2]:
# #             x_axis = np.append(x_axis,int(band[band.index('f')+1:-1]))
# #         else:
# #             x_axis = np.append(x_axis,int(band[band.index('f')+1:-1])*10)
# #         #============ next few steps is to do photometry in circle apertures
# #         if 'f160' not in band:
# #             image = fits.open(image_folder_path+'/'+band+'/'+band+'_'+image_index+'.fits')[0]
# #             elec_to_flux = float(image.header['PHOTFNU'])*10**6
# #             aftersubt_image,sky = subt_backnoise(image.data, fwhm(image.data),center_indices[0],center_indices[1])
# #             flux = 0
# #             for num_x in range(0,len(aftersubt_image)):
# #                 for num_y in range(0,len(aftersubt_image)):
# #                     if np.sqrt((num_x-float(center_indices[0]))**2+(num_y-float(center_indices[1]))**2)<=25:
# #                         flux += aftersubt_image[num_x][num_y]
# #             err = np.sqrt(pixel_number)*sky
# #             if '1' not in band[band.index('f')+1:band.index('f')+2] and '0'not in band[band.index('f')+1:band.index('f')+2]:
# #                 flux_list[int(band[band.index('f')+1:-1])] = flux*elec_to_flux,err*elec_to_flux
# #             else:
# #                 flux_list[int(band[band.index('f')+1:-1])*10] = flux*elec_to_flux,err*elec_to_flux
# #     #============ plot SED
# #     y_axis = np.zeros(0) # flux of every band
# #     errorbar = np.zeros(0)
# #     x_axis = np.sort(x_axis)
# #     for wavelength in x_axis:
# #         if wavelength in flux_list:
# #             y_axis = np.append(y_axis,flux_list[wavelength][0])
# #             errorbar = np.append(errorbar,flux_list[wavelength][1])
# #     # remove band without data
# #     x_axis = np.delete(x_axis,np.where(np.isnan(y_axis)))
# #     errorbar = np.delete(errorbar,np.where(np.isnan(y_axis)))
# #     y_axis = np.delete(y_axis,np.where(np.isnan(y_axis)))
# #     x_axis = np.delete(x_axis,np.where(y_axis<np.e**-3))
# #     errorbar = np.delete(errorbar,np.where(y_axis<np.e**-3))
# #     y_axis = np.delete(y_axis,np.where(y_axis<np.e**-3))
# #     # y_axis = -2.5*np.log(y_axis)+23.93
# #     plt.figure(figsize=(10,10))
# #     plt.errorbar(x_axis,y_axis,yerr=errorbar)
# #     plt.xticks(x_axis,rotation=90)
# #     plt.yticks()
# #     plt.xlabel('$\lambda (nm)$')
# #     plt.ylabel('$f_{\lambda}(\mu Jy)$')
# #     plt.title(image_index+' SED')
# #     plt.savefig('/Users/lpr/Data/fits/expdata/HST/goodsn_all/SED/'+image_index+'.eps')
# #     plt.close()
# #     print(stdfile+' SED is done')

# #=============================================================
# radec_table = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/radec_table.fits')[1]
# #=============================================================
# from scipy.integrate import quad
# from astropy.cosmology import FlatLambdaCDM as flcdm
# def dc(x,lambda0,m0):
#     hz=np.sqrt(lambda0+m0*(1+x)**3)
#     return 1/hz
# H0=70 # km/s/Mpc, nowadays commonly used
# Lambda0=0.7
# M0=0.3
# def kpc_per_arcsec(z):
#     angular_distance = flcdm(H0=70,Om0=0.3).angular_diameter_distance #((3*10**5)/H0)*quad(dc,0,z,args=(Lambda0,M0))[0]/(1+z)
#     arc_scale = angular_distance*np.pi*1000/(180*3600)
#     return(arc_scale)

# size = 25 # 1.5"
# for stdfile in stdima_order:
#     image_index = stdfile[stdfile.index('_',7)+1:stdfile.index('.fits')]
#     ima_160 = fits.open(image_folder_path+'/goodsn_f160w/goodsn_f160w_'+image_index+'.fits')[0]
#     index = np.where(radec_table.data['galaxy_id']==int(image_index))[0]
#     z = radec_table.data['catalog_redshift'][index]
#     ima_160.header['CTYPE1'] = 'RA---TAN-SIP'
#     ima_160.header['CTYPE2'] = 'DEC--TAN-SIP'
#     his = ima_160.header['HISTORY'][1]
#     cut_range = np.array([int(his[his.index('[')+2:his.index(':',16)]),int(his[his.index(',')+3:his.index(':',his.index(','))])])
#     ra,dec = float(radec_table.data['catalog_ra'][index]),float(radec_table.data['catalog_dec'][index])
#     center = np.array(WCS(ima_160.header).wcs_world2pix(ra,dec,0))-cut_range
#     scale = kpc_per_arcsec(z)
#     # pixelscale = float(ima_160.header['CD2_2'])
#     # ima_160.header['CTYPE1'] = 'RA---TAN-SIP'
#     # ima_160.header['CTYPE2'] = 'DEC--TAN-SIP'
#     # exptime_ = float(ima_160.header['EXPTIME'])
    
#     # try:
#     #     elec_to_flux = float(ima_160.header['PHOTFNU'])*10**6
#     # else KeyError:
#     #     elec_to_flux = float(ima_160.header['PHOTFLAM'])*10**29/
#     # aftersubt_image,sky = subt_backnoise(ima_160)
#     y_list = {}
#     error_list = {}
#     x_list = {}
#     # pixel_area = scale**2
#     # y_axis = np.zeros(0,dtype=float)
#     # error = np.zeros(0,dtype=float)
#     # np.seterr(divide='ignore')
#     # for radius in range(1,size+1):
#     #     flux = 0
#     #     pixel_number = 0
#     #     for num_x in range(0,aftersubt_image.shape[0]):
#     #         for num_y in range(0,aftersubt_image.shape[1]):
#     #             if np.sqrt((num_x-center[0])**2+(num_y-center[1])**2) >= radius-1 and np.sqrt((num_x-center[0])**2+(num_y-center[1])**2) <= radius:
#     #                 flux += aftersubt_image[num_x][num_y]
#     #                 pixel_number += 1
#     #     fsb = flux*elec_to_flux/(pixel_number*pixel_area)
#     #     sigmasb = sky/(np.sqrt(pixel_number)*pixel_area)
#     #     if fsb/sigmasb >= 5.:
#     #         y_axis = np.append(y_axis,fsb)
#     #         error = np.append(error,sigmasb/(fsb*np.log(10)))
#     # y_list['f160w'] = np.log10(y_axis)
#     # error_list['f160w'] = error
#     # x_list['f160w'] = np.arange(1,y_axis.shape[0]+1,1)
#     # band_list = ['f160w']
#     # for band in dirlist:
#     #     if 'f160' not in band:
#     #         band_list.append(band[7:])
#     #     if 'f160' not in band:
#     #         ima = fits.open(image_folder_path+'/'+band+'/'+band+'_'+image_index+'.fits')[0]
#     #         try:
#     #             elec_to_flux = float(ima_160.header['PHOTFNU'])*10**6
#     #         else KeyError:
#     #             elec_to_flux = float(ima_160.header['PHOTFLAM'])*10**29/
#     #         aftersubt_image,sky = subt_backnoise(ima)
#     #         y_axis = np.zeros(0,dtype=float)
#     #         error = np.zeros(0,dtype=float)
#     #         for radius in range(1,size+1):
#     #             flux = 0
#     #             pixel_number = 0
#     #             for num_x in range(0,aftersubt_image.shape[0]):
#     #                 for num_y in range(0,aftersubt_image.shape[1]):
#     #                     if np.sqrt((num_x-center[0])**2+(num_y-center[1])**2) >= radius-1 and np.sqrt((num_x-center[0])**2+(num_y-center[1])**2) <= radius:
#     #                         flux += aftersubt_image[num_x][num_y]
#     #                         pixel_number += 1
#     #             fsb = flux*elec_to_flux/(pixel_number*pixel_area)
#     #             sigmasb = sky/(np.sqrt(pixel_number)*pixel_area)
#     #             if fsb/sigmasb >= 5.:
#     #                 y_axis = np.append(y_axis,fsb)
#     #                 error = np.append(error,sigmasb/(fsb*np.log(10)))
#     #         if np.isnan(y_axis).any() or np.max(y_axis,initial=0)-np.min(y_axis,initial=0)<10**(-5) or y_axis.size == 0:
#     #             band_list.remove(band[7:])
#     #         else:
#     #             y_list[band[7:]] = np.log10(y_axis)
#     #             x_list[band[7:]] = np.arange(1,y_axis.shape[0]+1,1)
#     #             error_list[band[7:]] = error
#     band_list = []
#     for band in dirlist:
#         band_list.append(band[7:])
#         ima = fits.open(image_folder_path+'/'+band+'/'+band+'_'+image_index+'.fits')[0]
#         pixelscale = float(ima.header['CD2_2'])*3600
#         pixel_area = pixelscale*scale**2
#         elec_to_flux = ima.header['PHOTFLAM']*ima.header['PHOTPLAM']**2*1e6/(3e-5) # turn flux density to Jy
#         aftersubt_image,sky = subt_backnoise(ima)
#         y_axis = np.zeros(0,dtype=float)
#         error = np.zeros(0,dtype=float)
#         for radius in range(1,size+1):
#             flux = 0
#             pixel_number = 0
#             for num_x in range(0,aftersubt_image.shape[0]):
#                 for num_y in range(0,aftersubt_image.shape[1]):
#                     if np.sqrt((num_x-center[0])**2+(num_y-center[1])**2) >= radius-1 and np.sqrt((num_x-center[0])**2+(num_y-center[1])**2) <= radius:
#                         flux += aftersubt_image[num_x][num_y]
#                         pixel_number += 1
#             fsb = flux*elec_to_flux/(pixel_number*pixel_area)
#             sigmasb = sky/(np.sqrt(pixel_number)*pixel_area)
#             if fsb/sigmasb >= 5.:
#                 y_axis = np.append(y_axis,fsb)
#                 error = np.append(error,sigmasb/(fsb*np.log(10)))
#         if np.isnan(y_axis).any() or np.max(y_axis,initial=0)-np.min(y_axis,initial=0)<10**(-5) or y_axis.size == 0:
#             band_list.remove(band[7:])
#         else:
#             y_list[band[7:]] = np.log10(y_axis)
#             x_list[band[7:]] = np.arange(1,y_axis.shape[0]+1,1)*pixelscale*scale
#             error_list[band[7:]] = error
#     #===============
#     # annulus flux plot
#     c = ['k','b','r','g','c','m','y','sienna','gray','cyan','olive','orange','purple','pink','skyblue']
#     c_b = {}
#     c_n = 0
#     for band in dirlist:
#         c_b[band[7:]] = c[c_n]
#         c_n += 1
#     plt.figure(figsize=[10,6])
#     o = 0
#     for band in band_list:
#         plt.errorbar(x_list[band],y_list[band],yerr=error_list[band],color=c_b[band],label=band,lw=1,elinewidth=1,capsize=1,capthick=1)
#         o += 1
#     plt.legend(loc='upper left',bbox_to_anchor=(1,1))
#     plt.xlabel('r(kpc)')
#     plt.ylabel('annulus brightness($log10(\mu Jy/kpc^2$)')
#     plt.title(image_index+'annulus flux')
#     plt.savefig('/Users/lpr/Data/fits/expdata/HST/goodsn_all/SED/annulus'+image_index+'.eps')
#     plt.close()
#     print(image_index+' annulun is done')
#     # #===============
#     # # growth curve
#     # growth_list = {}
#     # for band in band_list:
#     #     growth = np.zeros(0,dtype=float)
#     #     for num in range(0,size):
#     #         if num == 0:
#     #             growth = np.append(growth,y_list[band][num])
#     #         else:
#     #             growth = np.append(growth,growth[num-1]+y_list[band][num])
#     #     growth_list[band] = growth
#     # plt.figure(figsize=[8,6])
#     # o = 0
#     # for band in band_list:
#     #     plt.plot(x_axis*pixelscale*3600*scale,growth_list[band],color=c[o],label=band)
#     #     o += 1
#     # plt.legend()
#     # plt.xlabel('r(kpc)')
#     # plt.ylabel('accumulated brightness($log10(\mu Jy/kpc^2$)')
#     # plt.xlim(0,6)
#     # plt.title('growth curve')
#     # plt.savefig('/Users/lpr/Data/fits/expdata/HST/goodsn_all/SED/growth'+image_index+'.eps')
#     # plt.close()
#     # print(image_index+' growth is done')
#     ########==============following steps write photometry of mine to fits
# # image_folder_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/multiband_cutimage'
# # image_path = image_folder_path+'/goodsn_f160w'
# # filelist = os.listdir(image_path)
# # image_file = []
# # for file in filelist:
# #     if file.endswith('.fits'):
# #         image_file.append(file)

# # radec_table = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/radec_table.fits')[1]
# # from scipy.integrate import quad
# # def dc(x,lambda0,m0):
# #     hz=np.sqrt(lambda0+m0*(1+x)**3)
# #     return 1/hz
# # H0=70 # km/s/Mpc, nowadays commonly used
# # Lambda0=0.7
# # M0=0.3
# # def kpc_per_arcsec(z):
# #     angular_distance=((3*10**5)/H0)*quad(dc,0,z,args=(Lambda0,M0))[0]/(1+z)
# #     arc_scale = angular_distance*np.pi*1000/(180*3600)
# #     return(arc_scale)
# # size = 25
# # flux_list = np.zeros([len(image_file),3])
# # num = 0
# # for file in image_file:
# #     image_index = file[file.index('_',7)+1:file.index('.fits')]
# #     galaxy_name = 'goodsn_f160w_'+image_index+'.fits'
# #     cata_index = np.where(radec_table.data['galaxy_id']==int(image_index))[0]
# #     image = fits.open(image_path+'/'+galaxy_name)[0]
# #     elec_to_flux = float(image.header['PHOTFNU'])*10**6
# #     image.header['CTYPE1'] = 'RA---TAN-SIP'
# #     image.header['CTYPE2'] = 'DEC--TAN-SIP'
# #     his = image.header['HISTORY'][1]
# #     cut_range = np.array([int(his[his.index('[')+2:his.index(':',16)]),int(his[his.index(',')+3:his.index(':',his.index(','))])])
# #     ra,dec = float(radec_table.data['catalog_ra'][cata_index]),float(radec_table.data['catalog_dec'][cata_index])
# #     center = np.array(WCS(image.header).wcs_world2pix(ra,dec,0))-cut_range
# #     flux = 0
# #     aftersubt_image,sky = subt_backnoise(image)
# #     # from last step, image array has been transposed
# #     pixel_number = 0
# #     for num_x in range(0,aftersubt_image.shape[0]):
# #         for num_y in range(0,aftersubt_image.shape[1]):
# #             if np.sqrt((num_x-float(center[0]))**2+(num_y-float(center[1]))**2)<=1.47:
# #                 flux += aftersubt_image[num_x][num_y]
# #                 pixel_number += 1
# #     err = np.sqrt(pixel_number)*sky*elec_to_flux
# #     flux = flux*elec_to_flux
# #     flux_list[num] = int(image_index),flux,err
# #     num += 1
# # col1 = fits.Column(name='galaxy_id', array=flux_list[:,0], format='K')
# # col2 = fits.Column(name='flux_160', array=flux_list[:,1], format='D')
# # col3 = fits.Column(name='fluxerr_160', array=flux_list[:,2], format='D')
# # flux_table = fits.BinTableHDU.from_columns([col1,col2,col3])
# # flux_table.writeto('/Users/lpr/Data/fits/expdata/HST/goodsn_all/flux_table.fits')
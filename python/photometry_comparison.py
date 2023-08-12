import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import os
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit

def subt_backnoise(image):
    hdu = np.transpose(image.data)
    hdu1 = hdu[np.where((hdu>-0.01)&(hdu<0.01))]
    noise = np.median(hdu1)
    return(hdu-noise,noise)
image_folder_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/multiband_cutimage'
image_path = image_folder_path+'/goodsn_f160w'
sed_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/SED/'
filelist = os.listdir(image_path)
image_file = []
for file in filelist:
    if file.endswith('.fits'):
        image_file.append(file)

radec_table = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/radec_table.fits')[1]
from scipy.integrate import quad
def dc(x,lambda0,m0):
    hz=np.sqrt(lambda0+m0*(1+x)**3)
    return 1/hz
H0=70 # km/s/Mpc, nowadays commonly used
Lambda0=0.7
M0=0.3
def kpc_per_arcsec(z):
    angular_distance=((3*10**5)/H0)*quad(dc,0,z,args=(Lambda0,M0))[0]/(1+z)
    arc_scale = angular_distance*np.pi*1000/(180*3600)
    return(arc_scale)
size = 25
flux_list = np.zeros([len(file),3])
num = 0
for file in image_file:
    image_index = file[file.index('_',7)+1:file.index('.fits')]
    galaxy_name = 'goodsn_f160w_'+image_index+'.fits'
    cata_index = np.where(radec_table.data['galaxy_id']==int(image_index))[0]
    image = fits.open(image_path+'/'+galaxy_name)[0]
    elec_to_flux = float(image.header['PHOTFNU'])*10**6
    image.header['CTYPE1'] = 'RA---TAN-SIP'
    image.header['CTYPE2'] = 'DEC--TAN-SIP'
    ra,dec = float(radec_table.data['catalog_ra'][cata_index]),float(radec_table.data['catalog_dec'][cata_index])
    center = WCS(image.header).wcs_world2pix(ra,dec,0)
    flux = 0
    aftersubt_image,sky = subt_backnoise(image)
    # from last step, image array has been transposed
    pixel_number = 0
    for num_x in range(0,aftersubt_image.shape[0]):
        for num_y in range(0,aftersubt_image.shape[1]):
            if np.sqrt((num_x-float(center[0]))**2+(num_y-float(center[1]))**2)<=1.47:
                flux += aftersubt_image[num_x][num_y]
                pixel_number += 1
    err = np.sqrt(pixel_number)*sky*elec_to_flux
    flux = flux*elec_to_flux
    flux_list[num] = int(image_index),flux,err
col1 = fits.Column(name='galaxy_id', array=flux_list[:,0], format='K')
col2 = fits.Column(name='flux_160', array=flux_list[:,1], format='D')
col3 = fits.Column(name='fluxerr_160', array=flux_list[:,2], format='D')
flux_table = fits.BinTableHDU.from_columns([col1,col2,col3])
flux_table.writeto('/Users/lpr/Data/fits/expdata/HST/goodsn_all/flux_table.fits')
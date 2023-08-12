import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import os
import matplotlib.pyplot as plt

def catalog_center(image_direc_path,cata_path):
    catalog = fits.open(cata_path)[1].data # catalog contains coordinates and redshift
    col = catalog.columns
    z_id = 0
    ra_id = 0
    dec_id = 0
    flux_id = 0
    for temp in range(0,len(col)):
            if col[temp].name == 'zbest':
                z_id = temp
            elif col[temp].name == 'RA_1':
                ra_id = temp
            elif col[temp].name == 'DEC_1':
                dec_id = temp
            elif col[temp].name == 'FLUX_APER_1_F160W':
                flux_id = temp
    filelist = os.listdir(image_direc_path)
    file = []
    for name in filelist:
        if name.endswith('.fits'):
            file.append(name)
    coordlist = np.zeros([len(file),8])
    for num in range(0,len(file)):
        # get the WCS  FK5 of image.maxpixel
        # index2: ra, index3:dec
        gala_idx = file[num][file[num].index('_')+1:file[num].index('_',7)]
        [ra,dec] = galaxy_center(image_direc_path + '/' + file[num])
        coord_catalog = np.array([0.,0.,0.,0.])
        for num2 in range(0,len(catalog)):
            if num2 == 0:
                distance = np.sqrt((catalog[num2][ra_id]-ra)**2 + (catalog[num2][dec_id]-dec)**2)
                coord_catalog[0] = catalog[num2][ra_id]
                coord_catalog[1] = catalog[num2][dec_id]
                coord_catalog[2] = catalog[num2][z_id]
                coord_catalog[3] = catalog[num2][flux_id]
            elif np.sqrt((catalog[num2][ra_id]-ra)**2 + (catalog[num2][dec_id]-dec)**2) <= distance:
                distance = np.sqrt((catalog[num2][ra_id]-ra)**2 + (catalog[num2][dec_id]-dec)**2)
                coord_catalog[0] = catalog[num2][ra_id]
                coord_catalog[1] = catalog[num2][dec_id]
                coord_catalog[2] = catalog[num2][z_id]
                coord_catalog[3] = catalog[num2][flux_id]
        coordlist[num] = int(gala_idx),coord_catalog[0],coord_catalog[1],ra,dec,coord_catalog[2],distance,coord_catalog[3]
    return coordlist

def galaxy_center(image):
    hdu_temp = fits.open(image)[0]
    hdu = np.transpose(hdu_temp.data) # astropy transpose x and y axis of fits
    rho_x = np.zeros(hdu.shape[0]) # marginal sum of x axis, namely, sum along j
    rho_y = np.zeros(hdu.shape[1]) # marginal sum of y axis, namely, sum along i
    for x in range(0,hdu.shape[0]):
        sum = 0
        for y in range(0,hdu.shape[1]):
            sum += hdu[x,y]
        rho_x[x] = sum
    for y in range(0,hdu.shape[1]):
        sum = 0
        for x in range(0,hdu.shape[0]):
            sum += hdu[x,y]
        rho_y[y] = sum
    meaninten_x = np.sum(rho_x)/hdu.shape[0]
    meaninten_y = np.sum(rho_y)/hdu.shape[1]
    #============
    # find new center coor_x
    for x in range(0,hdu.shape[0]):
        if rho_x[x] < meaninten_x:
            rho_x[x] = 0
        else:
            rho_x[x] = rho_x[x] - meaninten_x
    numerator = 0
    for x in range(0,len(rho_x)):
        numerator += rho_x[x]*(x+1) # make image begin with index 1
    coor_x = numerator/np.sum(rho_x)
    #============
    # find new center coor_y
    for y in range(0,hdu.shape[1]):
        if rho_y[y] < meaninten_y:
            rho_y[y] = 0
        else:
            rho_y[y] = rho_y[y] - meaninten_y
    numerator = 0
    for y in range(0,len(rho_y)):
        numerator += rho_y[y]*(y+1) # make image begin with indey 1
    coor_y = numerator/np.sum(rho_y)
    #============
    # convert new x_image and y_image into ra,dec
    hdu_temp.header['CTYPE1'] = 'RA---TAN-SIP'
    hdu_temp.header['CTYPE2'] = 'DEC--TAN-SIP'
    [ra,dec] = WCS(hdu_temp.header).wcs_pix2world(coor_x,coor_y,0)
    return [ra,dec]

def subt_backnoise(image):
    hdu = np.transpose(image.data)
    hdu1 = hdu[np.where((hdu>-0.01)&(hdu<0.01))]
    noise = np.median(hdu1)
    return(hdu-noise,noise)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
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

size = 25 # 1.5"
for stdfile in stdima_order[:1]:
    image_index = stdfile[stdfile.index('_',7)+1:stdfile.index('.fits')]
    ima_160 = fits.open(image_folder_path+'/goodsn_f160w/goodsn_f160w_'+image_index+'.fits')[0]
    index = np.where(radec_table.data['galaxy_id']==int(image_index))[0]
    z = radec_table.data['catalog_redshift'][index]
    ima_160.header['CTYPE1'] = 'RA---TAN-SIP'
    ima_160.header['CTYPE2'] = 'DEC--TAN-SIP'
    his = ima_160.header['HISTORY'][1]
    cut_range = np.array([int(his[his.index('[')+2:his.index(':',16)]),int(his[his.index(',')+3:his.index(':',his.index(','))])])
    ra,dec = float(radec_table.data['catalog_ra'][index]),float(radec_table.data['catalog_dec'][index])
    center = np.array(WCS(ima_160.header).wcs_world2pix(ra,dec,0))-cut_range
    scale = kpc_per_arcsec(z)
    y_list = {}
    error_list = {}
    x_list = {}
    band_list = []
    for band in dirlist:
        band_list.append(band[7:])
        ima = fits.open(image_folder_path+'/'+band+'/'+band+'_'+image_index+'.fits')[0]
        pixelscale = float(ima.header['CD2_2'])*3600
        pixel_area = (pixelscale*scale)**2
        exptime = float(ima_160.header['EXPTIME'])
        try:
            elec_to_flux = float(ima_160.header['PHOTFNU'])*10**6
        except KeyError:
            elec_to_flux = float(ima_160.header['PHOTFLAM'])*10**29/exptime
        aftersubt_image,sky = subt_backnoise(ima)
        y_axis = np.zeros(0,dtype=float)
        error = np.zeros(0,dtype=float)
        for radius in range(1,size+1):
            flux = 0
            pixel_number = 0
            for num_x in range(0,aftersubt_image.shape[0]):
                for num_y in range(0,aftersubt_image.shape[1]):
                    if np.sqrt((num_x-center[0])**2+(num_y-center[1])**2) >= radius-1 and np.sqrt((num_x-center[0])**2+(num_y-center[1])**2) <= radius:
                        flux += aftersubt_image[num_x][num_y]
                        pixel_number += 1
            fsb = flux*elec_to_flux/(pixel_number*pixel_area)
            sigmasb = sky/(np.sqrt(pixel_number)*pixel_area)
            if fsb/sigmasb >= 5.:
                y_axis = np.append(y_axis,fsb)
                error = np.append(error,sigmasb/(fsb*np.log(10)))
        if np.isnan(y_axis).any() or np.max(y_axis,initial=0)-np.min(y_axis,initial=0)<10**(-5) or y_axis.size == 0:
            band_list.remove(band[7:])
        else:
            y_list[band[7:]] = np.log10(y_axis)
            x_list[band[7:]] = np.arange(1,y_axis.shape[0]+1,1)*pixelscale*scale
            error_list[band[7:]] = error
    #===============
    # annulus flux plot
    c = ['k','b','r','g','c','m','y','sienna','gray','cyan','olive','orange','purple','pink','skyblue']
    c_b = {}
    c_n = 0
    for band in dirlist:
        c_b[band[7:]] = c[c_n]
        c_n += 1
    plt.figure(figsize=[10,6])
    o = 0
    for band in band_list:
        plt.errorbar(x_list[band],y_list[band],yerr=error_list[band],color=c_b[band],label=band,lw=1,elinewidth=1,capsize=1,capthick=1)
        o += 1
    plt.legend(loc='upper left',bbox_to_anchor=(1,1))
    plt.xlabel('r(kpc)')
    plt.ylabel('annulus brightness($log10(\mu Jy/kpc^2$)')
    plt.title(image_index+'annulus flux')
    plt.savefig('/Users/lpr/Data/fits/expdata/HST/goodsn_all/SED/annulus'+image_index+'.eps')
    plt.close()
    print(image_index+' annulun is done')


    mag_z_i = -2.5*alog10(convolve_i/convolve_denominator_i)
        mag_z_Y = -2.5*alog10(convolve_Y/convolve_denominator_Y)
        color = mag_z_i - mag_z_Y
        y_axis[num3] = color
        x_axis[num3] = z
# ==========================================================
# cut whole field image into single galaxy size with number
# drz.fits file is dizzled fits
# wht.fits file is sigma image(inverse variance)
# ==========================================================

from astropy.io import fits
import os
import datetime
import numpy as np                                                            
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs.utils import skycoord_to_pixel
import datetime
#from astropy.io.fits.Header import copy

# pripath = '/Users/lpr/Data/fits/pridata/CANDELS_image'
# expath = '/Users/lpr/Data/fits/expdata/HST/goodsn_all'
pripath = '/Users/lpr/Data/fits/pridata/CANDELS_image'
expath = '/Users/lpr/Data/fits/expdata/HST/goodss_all'
# firstly, i cut all drz.fits and sigma image into small images which have one galaxy in one image
drz_image = [] # all drizzled image of whole field
wht_image = [] # all sigma image of whole field
gdn_image = os.listdir(pripath+'/goodsn')
for items in gdn_image:
    if items.endswith('drz.fits'):
        drz_image.append(items)
    if items.endswith('wht.fits')
        wht_image.append(items)
print('===============drizzled images and sigma images load done===============')

# if dont want to use Cutout2D, then i need to use the header information of images in Huang(2021)
# smallfilelist = [] # smallfilelist is LIRGs images of Huang(2021)
# smalloriginal_file = os.listdir(pripath+'/goodsn_f125')
# for items in smalloriginal_file:
#     if items.endswith('.fits'):
#         smallfilelist.append(items)
# print('===============small image load done===============')

# if use Cutout2D, then i only need galaxy center in ra,dec unit
catalog = fits.open(expath+'/gdnFINALtable.fits')[1] # this binary table contains all LIRGs' ra,dec
CENarr = np.concatenate((np.reshape(catalog.data['ID_1'],(len(catalog.data),1)),np.reshape(catalog.data['RA_1'],(len(catalog.data),1)),np.reshape(catalog.data['DEC_1'],(len(catalog.data),1))),axis=1)
print('===============galaxy center catalog load done===============')

for big in drz_image[2:3]:
    band = big[big.index('f',18):big.index('_',19)] # HST band
    big_image = fits.open(pripath+'/goodsn/'+big)[0] # open every drizzled HST image in loop
    os.system('mkdir '+expath+'/singleGALAXYimage/gdn_'+band) # make new directory of every band
    print('==============='+band+' directory generation done===============')
    big_image.header['CTYPE1'] = 'RA---TAN-SIP'
    big_image.header['CTYPE2'] = 'DEC--TAN-SIP'
    fits.Header.remove(big_image.header,'CPDIS1',ignore_missing=True)
    fits.Header.remove(big_image.header,'CPDIS2',ignore_missing=True)
    big_image_copy = fits.ImageHDU.copy(big_image)
    for num in range(0,len(CENarr)): #
        gal_id = CENarr[num][0]
        ra = CENarr[num][1]
        dec = CENarr[num][2]
        big_image_copy.header = fits.Header.copy(big_image.header)
        position = SkyCoord(frame='fk5',unit='deg',ra=ra,dec=dec)
        x,y = skycoord_to_pixel(position,WCS(big_image.header),mode='wcs')
        cutout = Cutout2D(big_image.data,position,size=(101,101),wcs=WCS(big_image.header),mode='trim')
        big_image_copy.data = cutout.data
        big_image_copy.header.update(cutout.wcs.to_header())
        big_image_copy.header['HISTORY'] = 'EXTRACT: MON March 22nd 2021, ' + str(datetime.datetime.now().time())
        big_image_copy.header['HISTORY'] = 'EXTRACT range: '+str(int(x)-50)+':'+str(int(x)+50)+','+str(int(y)-50)+':'+str(int(y)+50)
        big_image_copy.writeto(expath+'/singleGALAXYimage/gdn_'+band+'/gdn_'+band+'_'+str(gal_id)+'.fits')
        print('==============='+band+'_'+str(gal_id)+' done===============')
    # for small in smallfilelist:
    #     small_image = fits.open('/Users/lpr/Data/fits/pridata/goodsn_f125/'+small)
    #     small_image_header = fits.getheader('/Users/lpr/Data/fits/pridata/goodsn_f125/'+small,0)
    #     small_image[0].header = big_image[0].header
    #     small_image[0].header.append('HISTORY',end=True)
    #     small_image[0].header['HISTORY'] = 'HEXTRACT: Mon Oct 9 2020'+ str(datetime.datetime.now().time())
    #     x_range = small_image_header['HISTORY'][2][small_image_header['HISTORY'][2].index('[')+1:small_image_header['HISTORY'][2].index(',')]
    #     y_range = small_image_header['HISTORY'][2][small_image_header['HISTORY'][2].index(',')+1:small_image_header['HISTORY'][2].index(']')]
    #     small_image[0].header['HISTORY'] = 'Extracted Image: '+ str([x_range,y_range])
    #     cut_image = big_image[0].data[int(y_range[0:y_range.index(':')]):int(y_range[y_range.index(':')+1:])+1,int(x_range[0:x_range.index(':')]):int(x_range[x_range.index(':')+1:])+1]
    #     small_image[0].data = cut_image
    #     small_image.writeto('/Users/lpr/Data/fits/expdata/HST/goodsn_all/multiband_cutimage/goodsn_'+band+'/goodsn_'+band+'_'+small[small.index('_')+1:small.index('_f125w')]+'.fits')
    #     del small_image[0].header['HISTORY']
    #     print(band+small[small.index('_')+1:small.index('_f125w')]+'.fits cut done')

# # f160w
# big_image = fits.open('/Users/lpr/Data/fits/pridata/goodsn/goodsn_all_wfc3_ir_f160w_060mas_v1.0_drz.fits')
# os.system('mkdir /Users/lpr/Data/fits/expdata/HST/goodsn_all/multiband_cutimage/goodsn_f160w')
# for small in smallfilelist:
#     small_image = fits.open('/Users/lpr/Data/fits/pridata/goodsn_f125/'+small)
#     small_image_header = fits.getheader('/Users/lpr/Data/fits/pridata/goodsn_f125/'+small,0)
#     small_image[0].header = big_image[0].header
#     small_image[0].header['HISTORY'] = 'HEXTRACT: Mon Oct 9 2020'+ str(datetime.datetime.now().time())
#     x_range = small_image_header['HISTORY'][2][small_image_header['HISTORY'][2].index('[')+1:small_image_header['HISTORY'][2].index(',')]
#     y_range = small_image_header['HISTORY'][2][small_image_header['HISTORY'][2].index(',')+1:small_image_header['HISTORY'][2].index(']')]
#     small_image[0].header['HISTORY'] = 'Extracted Image: '+ str([x_range,y_range])
#     cut_image = big_image[0].data[int(y_range[0:y_range.index(':')]):int(y_range[y_range.index(':')+1:])+1,int(x_range[0:x_range.index(':')]):int(x_range[x_range.index(':')+1:])+1]
#     small_image[0].data = cut_image
#     small_image.writeto('/Users/lpr/Data/fits/expdata/HST/goodsn_all/multiband_cutimage/goodsn_f160w/goodsn_f160w_'+small[small.index('_')+1:small.index('_f125w')]+'.fits')
#     del small_image[0].header['HISTORY']
#     print(small+' done')
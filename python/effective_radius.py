from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
def subt_backnoise(image):
    hdu = np.transpose(image.data)
    hdu1 = hdu[np.where((hdu>-0.01)&(hdu<0.01))]
    noise = np.median(hdu1)
    return(hdu-noise,noise)
# calculate total flux
hdu = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/multiband_cutimage/goodsn_f160w/goodsn_f160w_1092.fits')[0]
radec_table = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/tem_re_sersic2.fits')[1]
cata_index = np.where(radec_table.data['galaxy_id']==1092)[0]
ra,dec = float(radec_table.data['catalog_ra'][cata_index]),float(radec_table.data['catalog_dec'][cata_index])
his = hdu.header['HISTORY'][1]
cut_range = np.array([int(his[his.index('[')+2:his.index(':',16)]),int(his[his.index(',')+3:his.index(':',his.index(','))])])
center = np.array(WCS(hdu.header).wcs_world2pix(ra,dec,0))-cut_range
array,sky = subt_backnoise(hdu)
flux = 0
for num_x in range(0,len(array)):
    for num_y in range(0,len(array)):
        if np.sqrt((num_x-float(center[0]))**2+(num_y-float(center[1]))**2)<=28:
            flux += array[num_x,num_y]
# output effective radius
flux_re = 0
print('total flux: '+str(flux))
for re in range(0,len(array)):
    for x in range(0,len(array)):
        for y in range(0,len(array)):
            if np.sqrt((x-float(center[0]))**2+(y-float(center[1]))**2)<=re:
                flux_re += array[x,y]
    print(flux_re)
    if flux_re >= flux/2:
        break
print(re,flux_re)
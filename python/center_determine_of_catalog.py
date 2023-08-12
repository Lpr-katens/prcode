#======================
# find the matched van der Wel structure parameters from Huang coordinates
#======================

from astropy.io import fits
import os
import numpy as np
path_van = '/Users/lpr/Data/fits/pridata/CANDELS_catalog/3dhst/'
path_huang = '/Users/lpr/Data/fits/expdata/HST/'
field_list = ['goodsn','goodss','aegis']
for field in filelist:
    vanderwel = fits.open(path_van + field_list[field]+'_3dhst.v4.1_f160wf125w.fits')[1].data
    huang = fits.open(path_huang + field_list[field]+'all/' + field_list[field] + '_Huangall_van.fits')[1].data
    
    for num in range(0,len(file)):
        image = fits.open(cutimage_path + '/' + file[num])[0]
        # get the WCS  FK5 of image.maxpixel
        # index2: ra, index3:dec
        x,y = np.where(image.data == np.max(image.data))[0][0],np.where(image.data == np.max(image.data))[1][0]
        [ra,dec] = WCS(image.header).wcs_pix2world(y,x,0)
        coord_catalog = np.array([0.,0.])
        for num2 in range(0,len(catalog)):
            if num2 == 0:
                distance = np.sqrt((catalog[num2][2]-ra)**2 + (catalog[num2][3]-dec)**2)
                coord_catalog[0] = catalog[num2][2]
                coord_catalog[1] = catalog[num2][3]
            elif np.sqrt((catalog[num2][2]-ra)**2 + (catalog[num2][3]-dec)**2) <= distance:
                distance = np.sqrt((catalog[num2][2]-ra)**2 + (catalog[num2][3]-dec)**2)
                coord_catalog[0] = catalog[num2][2]
                coord_catalog[1] = catalog[num2][3]
        coordlist[file[num]] = np.array([[coord_catalog[0],coord_catalog[1]],[ra,dec],distance])
        print(file[num]+' ra,dec is append to file')
    
    # save center of catalog
    file = open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/center_of_cat.cat','w')
    file.write(str(coordlist))
    file.close()
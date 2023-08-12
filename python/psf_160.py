from astropy.io import fits
import math
import matplotlib.pyplot as plt
import numpy as np

# ori_ima = fits.open('/Users/lpr/Data/lirg_project/intake/CANDELS/goodsn/goodsn_all_acs_wfc_f814w_060mas_v2.0_drz.fits')[0]
catalog = fits.open('/Users/lpr/Data/SExtractor/test/goodsn_f814w_psf.fits')[1].data
# # ind = np.where(catalog['CLASS_STAR']>0.9)
# # catalog = catalog[ind]
# # coor_line = coor_file.readlines()
# # del coor_line[:4] # now coor_line store all FWHM, MAG, COORDINATE
# coor_x = [] # x coordinate of point source
# coor_y = [] # y coordinate of point source
# fwhm = [] # full width at half-maximum of point sources
# mag = [] # magnitude of point source
# idx = []
# for num in range(0,len(catalog)):
#     x = catalog[num]['X_IMAGE'] # store the third line of .cat
#     y = catalog[num]['Y_IMAGE'] # store the forth line of .cat
#     fwhm_1 = catalog[num]['FWHM_IMAGE']
#     mag_1 = catalog[num]['MAG_AUTO']
#     idx_1 = catalog[num]['NUMBER']
#     if 1.5<fwhm_1<2.4 and -6.3<mag_1<-3: 
#             coor_x.append(float(x))
#             coor_y.append(float(y))
#             fwhm.append(float(fwhm_1))
#             mag.append(float(mag_1))
#             idx.append(idx_1)
# print('point source coordinates are done')
# # FWHM = 2*sqrt(2*ln2), cut to 10*sigma
# # FWHM = 3.655
# # sigma = FWHM/(2*np.sqrt(2*np.log(2)))
# # ini_size = 10*int(sigma)*2+1
ini_size = 69
# # image = np.zeros([len(idx),ini_size,ini_size]) # image date((10*sigma*2+1)x(10*sigma*2+1) pixels) from every point source
# #====== get all image data of point source
# for num in range(0,len(idx)):
#     # +1 is because coordinate's index is different with python, begin with 1
#     x_range = [int(coor_y[num])-int((ini_size-1)/2),int(coor_y[num])+int((ini_size-1)/2)+1]
#     y_range = [int(coor_x[num])-int((ini_size-1)/2),int(coor_x[num])+int((ini_size-1)/2)+1]
#     image = ori_ima.data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
#     ini_image = fits.PrimaryHDU(image)
#     ini_hdu = fits.HDUList([ini_image])
#     ini_hdu_name = str(idx[num])+'_ps.fits'
#     ini_hdu.writeto('/Users/lpr/Data/SExtractor/test/goodsn_psf/'+ini_hdu_name,overwrite=True)
#     plt.figure()
#     ax = plt.axes(projection ='3d')
#     x = np.arange(0,image.shape[0])
#     y = np.arange(0,image.shape[1])
#     [x,y] = np.meshgrid(x,y)
#     threeD_surf = str(idx[num])+'_surface.eps'
#     ax.plot_surface(x,y,image)
#     plt.savefig('/Users/lpr/Data/SExtractor/test/goodsn_psf/'+threeD_surf)
#     plt.close()
#     # image[num] = ori_ima.data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
#     # fits.open.data transform the x and y axis of image
# print('image cut is done')

# #====== export all point sources selected
# for num in range(0,len(idx)):
#     ini_image = fits.PrimaryHDU(image[num])
#     ini_hdu = fits.HDUList([ini_image])
#     ini_hdu_name = str(num)+'_ps.fits'
#     ini_hdu.writeto('/Users/lpr/Data/SExtractor/test/goodsn_psf/'+ini_hdu_name,overwrite=True)
# print('selected images export done')

# for num in range(0,len(image)):
#     # plot brigtness 3D surface of every selected point source
#     plt.figure()
#     ax = plt.axes(projection ='3d')
#     x = np.arange(0,image[num].shape[0])
#     y = np.arange(0,image[num].shape[0])
#     [x,y] = np.meshgrid(x,y)
#     threeD_surf = str(num)+'_surface.eps'
#     ax.plot_surface(x,y,image[num])
#     plt.savefig('/Users/lpr/Data/SExtractor/test/goodsn_psf/'+threeD_surf)
#     plt.close()
# print('surface image export done')

#====== selected 20 sources
indices = [784,1454,7056,8277,8345,9444,9919,10250,10599,12117,12965,13479,13937,14341,14939,18439,19417,20991]
subsize=69
#====== difine a list of fixed number and size
def gensubpixel(number,size):
    new_list = np.zeros([size[0],size[1]])
    for x in range(0,size[0]):
          for y in range(0,size[1]):
              new_list[x][y] = number/(size[0]*size[1])
    return new_list

#====== get new image data(size is ((10*sigma*2+1)*101)^2) of 
#====== corrected point source center
new_image = np.zeros([len(indices),(ini_size)*subsize-(subsize-1),(ini_size)*subsize-(subsize-1)])
selected_num = 0
for num in indices:
  pixel_sub = np.zeros([(ini_size)*subsize,(ini_size)*subsize]) # store every sub pixel data
  image = fits.open('/Users/lpr/Data/SExtractor/test/goodsn_psf/'+str(num)+'_ps.fits')[0].data
  coor_x = catalog[np.where(catalog['NUMBER']==num)]['X_IMAGE']
  coor_y = catalog[np.where(catalog['NUMBER']==num)]['Y_IMAGE']
  for x in range(0,image.shape[0]):
          for y in range(0,image.shape[1]):
              pixel_sub[x*subsize:(x+1)*subsize,y*subsize:(y+1)*subsize]= gensubpixel(image[x][y],[subsize,subsize])
  # math.modf()=[fractional part,integer part]
  center_x = ((ini_size-1)/2)*subsize+int(math.modf(coor_y)[0]*subsize)
  center_y = ((ini_size-1)/2)*subsize+int(math.modf(coor_x)[0]*subsize)
  new_image[selected_num] = pixel_sub[int(center_x)-int((ini_size-1)/2)*subsize:int(center_x)+int((ini_size-1)/2)*subsize+1,int(center_y)-int((ini_size-1)/2)*subsize:int(center_y)+int((ini_size-1)/2)*subsize+1]
  selected_num += 1
print('new image with right center is done')

# plot fwhm distribution diagram
# fwhm_final = []
# for num in indices:
#     fwhm_final.append(fwhm[num])
# plt.figure()
# plt.hist(fwhm_final)
# plt.xlabel('fwhm')
# plt.ylabel('number')
# plt.savefig('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/goodsn_all/psf/psf_160_histgram.eps')
# plt.close()
# print('histgram image export is done')

#====== stack image
stack_image = np.zeros([new_image.shape[1],new_image.shape[2]])
for num in range(0,len(new_image)):
    stack_image += new_image[num]/np.sum(new_image[num])
print('image stack is done')

#====== change image size back to 101*101
change_orimage = np.zeros([subsize,subsize])
# ini_index = 1010-50-9*101
for x_pixel in range(0,subsize):
    for y_pixel in range(0,subsize):
        change_orimage[x_pixel,y_pixel] = np.sum(stack_image[subsize*x_pixel:subsize*(x_pixel+1),subsize*y_pixel:subsize*(y_pixel+1)])
#====== save psf image as fits
data = fits.PrimaryHDU(change_orimage)
hdu = fits.HDUList([data])
hdu.writeto('/Users/lpr/Data/SExtractor/test/goodsn_psf/goodsn_f814w_psf.fits')
print('psf with .fits is done')
#======= show psf image
# plt.figure()
# plt.imshow(change_orimage,interpolation=None,origin='lower')
# plt.xlabel('x axis')
# plt.ylabel('y axis')
# plt.colorbar()
# plt.savefig('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/goodsn_all/psf/psf_160_stackimage.eps')
# plt.close()
# print('psf stack image is done')

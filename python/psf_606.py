from astropy.io import fits
import math
import matplotlib.pyplot as plt
import numpy as np

ori_ima = fits.open('/Users/lpr/Data/fits/pridata/goodsn/goodsn_all_acs_wfc_f606w_060mas_v2.0_drz.fits')[0]
coor_file = open('/Users/lpr/Data/SExtractor/test/606_psf.cat')
coor_line = coor_file.readlines()
del coor_line[:4] # now coor_line store all FWHM, MAG, COORDINATE

coor_x = [] # x coordinate of point source
coor_y = [] # y coordinate of point source
fwhm = [] # full width at half-maximum of point sources
mag = [] # magnitude of point source
for num in range(0,len(coor_line)):
    x = coor_line[num].split()[2] # store the third line of .cat
    y = coor_line[num].split()[3] # store the forth line of .cat
    fwhm_1 = coor_line[num].split()[0]
    mag_1 = coor_line[num].split()[1]
    if 1.793<float(coor_line[num].split()[0])<2.414 and -7.14<float(coor_line[num].split()[1])<-3.03 : 
            coor_x.append(float(x))
            coor_y.append(float(y))
            fwhm.append(float(fwhm_1))
            mag.append(float(mag_1))
print('point source coordinates are done')

# FWHM = 2*sqrt(2*ln2), cut to 10*sigma
FWHM = 2.414
sigma = FWHM/(2*np.sqrt(2*np.log(2)))
ini_size = 10*int(sigma)*2+1
image = np.zeros([len(coor_x),ini_size,ini_size]) # image date((10*sigma*2+1)x(10*sigma*2+1) pixels) from every point source
#====== get all image data of point source
for num in range(0,len(coor_x)):
    # +1 is because coordinate's index is different with python, begin with 1
    x_range = [int(coor_y[num])-int((ini_size-1)/2),int(coor_y[num])+int((ini_size-1)/2)+1]
    y_range = [int(coor_x[num])-int((ini_size-1)/2),int(coor_x[num])+int((ini_size-1)/2)+1]
    image[num] = ori_ima.data[x_range[0]:x_range[1],y_range[0]:y_range[1]]
    # fits.open.data transform the x and y axis of image
print('image cut is done')

#====== export all point sources selected
for num in range(0,len(image)):
    ini_image = fits.PrimaryHDU(image[num])
    ini_hdu = fits.HDUList([ini_image])
    ini_hdu_name = str(num)+'_point_source_606.fits'
    ini_hdu.writeto('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/goodsn_all/psf/ini_point_source/606/'+ini_hdu_name)
print('selected images export done')

for num in range(0,len(image)):
    # plot brigtness 3D surface of every selected point source
    plt.figure()
    ax = plt.axes(projection ='3d')
    x = np.arange(0,image[num].shape[0])
    y = np.arange(0,image[num].shape[0])
    [x,y] = np.meshgrid(x,y)
    threeD_surf = str(num)+'_surface.eps'
    ax.plot_surface(x,y,image[num])
    plt.savefig('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/goodsn_all/psf/ini_point_surface/606/'+threeD_surf)
    plt.close()
print('surface image export done')


# #====== selected 20 sources
indices = [167,180,179,177,176,175,173,171,169,160,153,151,141,138,137,135,134,129,124,121,112,106,98,95,89,77,68,11]

#====== difine a list of fixed number and size
def gensubpixel(number,size):
    new_list = np.zeros([size[0],size[1]])
    for x in range(0,size[0]):
          for y in range(0,size[1]):
              new_list[x][y] = number/(size[0]*size[1])
    return new_list

#====== get new image data(size is ((10*sigma*2+1)*101)^2) of 
#====== corrected point source center
new_image = np.zeros([len(indices),(ini_size)*101-100,(ini_size)*101-100])
selected_num = 0
for num in indices:
  pixel_sub = np.zeros([(ini_size)*101,(ini_size)*101]) # store every sub pixel data
  for x in range(0,len(image[num])):
          for y in range(0,len(image[num])):
              pixel_sub[x*101:(x+1)*101,y*101:(y+1)*101]= gensubpixel(image[num][x][y],[101,101])
  # math.modf()=[fractional part,integer part]
  center_x = ((ini_size-1)/2)*101+int(math.modf(coor_y[num])[0]*101)
  center_y = ((ini_size-1)/2)*101+int(math.modf(coor_x[num])[0]*101)
  new_image[selected_num] = pixel_sub[int(center_x)-int((ini_size-1)/2)*101:int(center_x)+int((ini_size-1)/2)*101+1,int(center_y)-int((ini_size-1)/2)*101:int(center_y)+int((ini_size-1)/2)*101+1]
  selected_num += 1
print('new image with right center is done')

# plot fwhm distribution diagram
fwhm_final = []
for num in indices:
    fwhm_final.append(fwhm[num])
plt.figure()
plt.hist(fwhm_final)
plt.xlabel('fwhm')
plt.ylabel('number')
plt.savefig('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/goodsn_all/psf/psf_606_histgram.eps')
plt.close()
print('histgram image export is done')

#====== stack image
stack_image = np.zeros([new_image.shape[1],new_image.shape[2]])
for num in range(0,len(new_image)):
    stack_image += new_image[num]/np.sum(new_image[num])
print('image stack is done')

#====== change image size back to 101*101
change_orimage = np.zeros([19,19])
ini_index = 1010-50-9*101
for x_pixel in range(0,19):
    for y_pixel in range(0,19):
        change_orimage[x_pixel,y_pixel] = np.sum(stack_image[ini_index+101*x_pixel:ini_index+101*(x_pixel+1),ini_index+101*y_pixel:ini_index+101*(y_pixel+1)])
#====== save psf image as fits
data = fits.PrimaryHDU(change_orimage)
hdu = fits.HDUList([data])
hdu.writeto('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/goodsn_all/psf/psf_606_stackimage.fits')
print('psf with .fits is done')
#======= show psf image
plt.figure()
plt.imshow(change_orimage,interpolation=None,origin='lower')
plt.xlabel('x axis')
plt.ylabel('y axis')
plt.colorbar()
plt.savefig('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/goodsn_all/psf/psf_606_stackimage.eps')
plt.close()
print('psf stack image is done')

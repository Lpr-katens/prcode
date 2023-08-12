from astropy.io import fits
import numpy as np

hdu = fits.open('/Users/lpr/Data/fits/pridata/goodsn_f160/goodsn_1092_f160w.fits')[0].data
hdu = np.transpose(hdu)
y = hdu[50]
y_max = np.max(y)
x_max = np.where(y==y_max)
x_index = []
print(y_max/2,x_max)
for num in range(1,len(y)):
    if y[num-1] < y_max/2 and y_max/2 <= y[num]:
        x_index.append((2*num+1)/2)
        print(y[num-1],y[num])
        print((2*num+1)/2)
    elif y[num-1] >= y_max/2 and y_max/2 > y[num]:
        x_index.append((2*num+1)/2)
        print(y[num-1],y[num])
        print((2*num+1)/2)
print('FWHM='+str(x_index[1]-x_index[0]))
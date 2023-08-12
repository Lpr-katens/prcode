import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

image = '/Users/lpr/Data/fits/pridata/goodsn_f160/goodsn_1092_f160w.fits'
fwhm = 8
cen_x,cen_y = 50,50
hdu1 = fits.open(image)[0]
hdu = np.transpose(hdu1.data)
rad_low = int(fwhm/2)*6
rad_high = int(fwhm/2)*7
noi = []
for x in range(0,hdu.shape[0]):
    for y in range(0,hdu.shape[1]):
        if np.sqrt((x-cen_x)**2+(y-cen_y)**2) >= rad_low and np.sqrt((x-cen_x)**2+(y-cen_y)**2) <= rad_high:
            noi.append(hdu[x,y])
ydata,t = np.histogram(noi,bins=100)
plt.hist(noi,bins=100)
xdata = t[1:] # t include begining and end of range, we use end instead of beginning
def gaussian(x,mu,sigma,param):
    return param*np.exp(-((x-mu)/sigma)**2)/(sigma*np.sqrt(2*np.pi))
popt,pcov = curve_fit(gaussian,xdata,ydata) # popt is the optimal values for the parameters
# take mu as background noise
print(popt)
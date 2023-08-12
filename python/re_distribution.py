from scipy.integrate import quad
import matplotlib.pyplot as plt
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
from astropy.io import fits
import numpy as np
hdu = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/gdn_match.fits')[1].data
hdu = hdu[np.where(hdu['re_f160w']>0)]
for num in range(0,len(hdu)):
    hdu[num]['re_f160w'] = hdu[num]['re_f160w']*kpc_per_arcsec(hdu[num]['zbest'])
agn = hdu[np.where(hdu['TMP_CLASS'][:,0] == 1)]
composite = hdu[np.where(hdu['TMP_CLASS'][:,0] == 2)]
sf = hdu[np.where(hdu['TMP_CLASS'][:,0] == 3)]
quiescent = hdu[np.where(hdu['TMP_CLASS'][:,0] == 4)]
class5 = hdu[np.where(hdu['TMP_CLASS'][:,0] == 5)]
plt.hist(hdu['re_f160w'],bins=25,histtype='step',color='r',label='all')
plt.hist(agn['re_f160w'],bins=25,histtype='step',color='b',label='agn')
plt.hist(composite['re_f160w'],bins=25,histtype='step',color='g',label='composite')
plt.hist(sf['re_f160w'],bins=25,histtype='step',color='y',label='star forming')
plt.hist(quiescent['re_f160w'],bins=25,histtype='step',color='cyan',label='quiescent')
plt.hist(class5['re_f160w'],bins=25,histtype='step',color='k',label='class5')
plt.legend()
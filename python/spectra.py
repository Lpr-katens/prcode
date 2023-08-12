import os
from astropy.io import fits
import numpy as np
from lpr.spectrum import transframe
import matplotlib.pyplot as plt
# # ==============================================================
# # modify spectra from observed frame to rest frame
# # ==============================================================
path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/spectral/'
expath = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/spectral/restframe/'
catalog = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/goodsn_Huang_van_3dhst.fits')[1].data
spectra_list = []
for file in os.listdir(expath):
    if file.endswith('.fits') and file.startswith('res_mod_'):
        spectra_list.append(file)
for file in spectra_list:
    spectrum = fits.open(path+file)
    wave_obs = spectrum[1].data['wave']
    flux = spectrum[1].data['flux_mean3']
    idx = int(file[file.index('_',13)+1:file.rindex('_')])
    z = catalog[np.where(catalog['ID_van']==idx)[0][0]]['REDSHIFT']
    wave_res = transframe.obs2res(wave_obs,z)
    col = fits.Column(name='wave_res',array=wave_res,format='D',unit='ANGSTROMS')
    table = fits.BinTableHDU.from_columns(spectrum[1].data.columns+col)
    hdu = fits.HDUList([fits.PrimaryHDU(),table])
    hdu.writeto(path+'/restframe/res_'+file,overwrite=True)
    spectrum.close()
    hdu.close()
    plt.figure(figsize=[10,7])
    plt.plot(wave_res,flux,linewidth=1)
    plt.text(1,1,'z='+str(np.around(z,3)),color='red',transform=plt.gca().transAxes)
    plt.xlabel('$Å$',fontsize=10)
    plt.ylabel('electrons/s',fontsize=10)
    plt.title('spectrum')
    plt.savefig(path+'/restframe/res_'+file[:-4]+'pdf')
    plt.close()
    print('----------------- '+file+' is done--------------')
# ==============================================================
# add H_alpha indicator to plot
# ==============================================================
# path = '/Users/lpr/Data/fits/expdata/HST/goodss_all/spectral/restframe/'
# catalog = fits.open('/Users/lpr/Data/fits/expdata/HST/goodss_all/goodss_Huang_van_3dhstFITemissionline.fits')[1].data
# temp = os.listdir(path)
# spectra_list = []
# for file in temp:
#     if file.endswith('.fits'):
#         spectra_list.append(file)
H_alpha = 6562.5 # anstrong
v = 0.01 # velocity is 1000 km/s relative to central wavelength
lambda1 = H_alpha*(3/(3+v))
lambda2 = H_alpha*(3/(3-v))
for file in spectra_list:
    idx = int(file[file.index('_',17)+1:file.rindex('_')])
    print(idx)
    z = catalog[np.where(catalog['ID_3dhst']==idx)[0][0]]['z_used']
    spectrum = fits.open(expath+file)[1].data
    wave = spectrum['wave_res']
    flux = spectrum['flux_mean3']
    if lambda1>wave[0] and lambda1<wave[-1] and lambda2>wave[0] and lambda2<wave[-1]:
        plt.figure(figsize=[10,7])
        plt.plot(wave,flux,linewidth=1,color='blue')
        plt.text(1,1,'z='+str(np.around(z,3)),color='red',transform=plt.gca().transAxes)
        plt.xlabel('$Å$',fontsize=10)
        plt.ylabel('electrons/s',fontsize=10)
        plt.title('spectrum')
        plt.axvline(x=lambda1,color='red',linewidth=0.5)
        plt.axvline(x=lambda2,color='red',linewidth=0.5)
        print(file[:-5]+'_ha.pdf')
        plt.savefig(expath+file[:-5]+'_ha.pdf')
        plt.close()
        print('----------------- '+file+' is done -----------------')

for file in spectra_list:
    spectrum = fits.open(path+file)[1].data
    wave = spectrum['wave']
    flux = spectrum['flux']
    array = np.full_like(flux,-999)
    for num1 in range(0,len(flux)):
        if num1 == 0:
            array[num1] = flux[num1]
        elif num1 == len(flux) - 1:
            array[num1] = flux[num1]
        else:
            array[num1] = np.mean(flux[num1]+flux[num1-1]+flux[num1+1])
    col1 = fits.Column(name='wave',array=wave,format='K')
    col2 = fits.Column(name='flux_mean3',array=array,format='D')
    hdu = fits.BinTableHDU.from_columns([col1,col2])
    hdu.writeto(path+'/mod_'+file,overwrite=True)
    print(file+' is done')
H_alpha = 6562.5 # anstrong
v = 0.01 # velocity is 1000 km/s relative to central wavelength
lambda1 = H_alpha*(3/(3+v))
lambda2 = H_alpha*(3/(3-v))
for file in spectra_list:
    idx = int(file[file.index('_',13)+1:file.rindex('_')])
    z = catalog[np.where(catalog['ID_3dhst']==idx)[0][0]]['z_used']
    spectrum = fits.open(expath+file)[1].data
    wave = spectrum['wave']
    flux = spectrum['flux']
    if lambda1>wave[0] and lambda1<wave[-1] and lambda2>wave[0] and lambda2<wave[-1]:
        plt.figure(figsize=[10,7])
        plt.plot(wave,flux,linewidth=1,color='blue')
        plt.text(1,1,'z='+str(np.around(z,3)),color='red',transform=plt.gca().transAxes)
        plt.xlabel('$Å$',fontsize=10)
        plt.ylabel('electrons/s',fontsize=10)
        plt.title('spectrum')
        plt.axvline(x=lambda1,color='red',linewidth=0.5)
        plt.axvline(x=lambda2,color='red',linewidth=0.5)
        print(file[:-5]+'_ha.eps')
        plt.savefig(path+file[:-5]+'_ha.eps')
        plt.close()
        print('----------------- '+file+' is done -----------------')
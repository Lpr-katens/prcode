#==========
# do pypher, generate kernel
import os
import subprocess as sp

filelist_606 = os.listdir('/Users/lpr/Data/fits/pridata/goodsn_f606')
filelist_160 = os.listdir('/Users/lpr/Data/fits/pridata/goodsn_f160')

file_606 = []
file_160 = []
for file in filelist_160:
 	if file.endswith('.fits'):
         file_160.append(file)
for file in filelist_606:
 	if file.endswith('.fits'):
         file_606.append(file)

for num in range(0,len(file_160)):
 	tar = '/Users/lpr/Data/fits/pridata/goodsn_f160/' + file_160[num]
 	sou = '/Users/lpr/Data/fits/pridata/goodsn_f606/' + 'goodsn_' + file_160[num][7:file_160[num].index('f160w')-1] + '_f606w.fits'
 	ker = '/Users/lpr/Data/fits/expdata/CONVOLIMAGE/pypher/kernel_' + file_160[num][7:file_160[num].index('f160w')-1] + '.fits'
 	pyp = 'pypher ' + sou +  ' ' + tar + ' ' + ker
 	sp.run(pyp,shell=True,check=True)

#==========
# do convolution
from astropy.io import fits
from astropy.convolution import convolve

filelist_ker = os.listdir('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/pypher')
kernel_list = []
for file in filelist_ker:
    if file.endswith('.fits'):
        kernel_list.append(file)
for num in range(0,len(file_606)):
    sou = fits.open('/Users/lpr/Data/fits/pridata/goodsn_f606/'+file_606[num])[0]
    objid = 'kernel_'+file_606[num][7:file_606[num].index('f606w')-1]+'.fits'
    index = kernel_list.index(objid)
    kernel= fits.open('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/pypher/'+kernel_list[index])[0]
    conv = convolve(sou.data,kernel.data)
    data = fits.PrimaryHDU(conv)
    hdu = fits.HDUList([data])
    name = 'conv_'+file_606[num][7:file_606[num].index('f606w')-1]+'.fits'
    hdu.writeto('/Users/lpr/Data/fits/expdata/CONVOLIMAGE/pypher/'+name)
    print(name+' is done')
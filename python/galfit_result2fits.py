# # ==================================================
# # this code store galfit results in a fits file
# # ==================================================
# import os
# import numpy as np
# from astropy.io import fits

# path = '/Users/lpr/Data/galfit/LIRGproject'
# temp = os.listdir(path)
# gala = []
# for file in temp:
# 	if file.endswith('.DS_Store'):
# 		print(file)
# 	else:
# 		gala.append(file)
# print('-------------------------- ready to do for loop --------------------------')
# array = np.ones([len(gala),7]) # array contains galaxy id, bulge mag, re, disk mag, rs, B/T
# for num in range(0,len(gala)):
# 	temp = os.listdir(path+'/'+gala[num])
# 	for file in temp:
# 		if file == 'galfit.02' or file == 'galfit.01':
# 			f = open(path+'/'+gala[num]+'/'+file)
# 			data = f.readlines()
# 			for num2 in range(0,len(data)): 
# 				if 'Component number: 1' in data[num2]: 
# 				    mag_b = data[num2+3] 
# 				    re = data[num2+4] 
# 				elif 'Component number: 2' in data[num2]: 
# 				    mag_d = data[num2+3] 
# 				    rs = data[num2+4]
# 			MAG_B = float(mag_b[mag_b.index(')')+2:mag_b.index('    1')-1])
# 			RE = float(re[re.index(')')+2:re.index('    1')-1])
# 			MAG_D = float(mag_d[mag_d.index(')')+2:mag_d.index('    1')-1])
# 			RS = float(rs[rs.index(')')+2:rs.index('    1')-1])
# 			DB = 10**((MAG_B-MAG_D)/2.5)
# 			BT = 1/(DB+1)
# 			array[num]=int(gala[num]),MAG_B,RE,MAG_D,RS,DB,BT
# print('-------------------------- ready to save FITS file --------------------------')
# col1 = fits.Column(name='galaxy_id',array=array[:,0],format='K')
# col2 = fits.Column(name='MAG_bulge',array=array[:,1],format='D')
# col3 = fits.Column(name='re',array=array[:,2],format='D')
# col4 = fits.Column(name='MAG_disk',array=array[:,3],format='D')
# col5 = fits.Column(name='rs',array=array[:,4],format='D')
# col6 = fits.Column(name='disk2bulge',array=array[:,5],format='D')
# col7 = fits.Column(name='bulge2total',array=array[:,6],format='D')
# hdu = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7])
# hdu.writeto(path+'/galfit_result.fits')
# print('-------------------------- galfit result has been saved --------------------------')
# for num2 in range(0,len(data)): 
# 				if 'devauc' in data[num2]: 
# 				    MAG_B = float(data[num2][data[num2].index(')')+4:data[num2].index(')')+9].replace('*',' '))
# 				    MAG_B_err = float(data[num2+1][data[num2+1].index(')')+4:data[num2+1].index(')')+9].replace('*',' '))
# 				    re = float(data[num2][data[num2].index(')')+10:data[num2].index('.',42)+5].replace('*',' ')) 
# 				    re_err = float(data[num2+1][data[num2+1].index(')')+10:data[num2+1].index('.',42)+5].replace('*',' '))
# 				elif 'disk' in data[num2]: 
# 				    MAG_D = float(data[num2][data[num2].index(')')+4:data[num2].index(')')+9].replace('*',' '))
# 				    MAG_D_err = float(data[num2+1][data[num2+1].index(')')+4:data[num2+1].index(')')+9].replace('*',' '))
# 				    rs = float(data[num2][data[num2].index(')')+10:data[num2].index('.',42)+5].replace('*',' ')) 
# 				    rs_err = float(data[num2+1][data[num2+1].index(')')+10:data[num2+1].index('.',42)+5].replace('*',' '))
# 				elif 'Chi^2/nu' in data[num2]:
# 					chi_nu = float(data[num2][data[num2].index('=')+2:].replace('*',' '))
# -------------------------------------------- use fit.log instead of galfit.01 --------------------------------------------
import os
import numpy as np
from astropy.io import fits

path = '/Users/lpr/Data/galfit/LIRGproject/goodss'
temp = os.listdir(path)
gala = []
for file in temp:
	if file.endswith('.DS_Store') or file.endswith('.fits'):
		print(file)
	else:
		gala.append(file)
print('-------------------------- ready to do for loop --------------------------')
array = np.ones([len(gala),16]) # array contains galaxy id, bulge mag, re, disk mag, rs, B/T
for num in range(0,len(gala)):
	temp = os.listdir(path+'/'+gala[num])
	for file in temp:
		if file == 'fit.log':
			f = open(path+'/'+gala[num]+'/'+file)
			data = f.readlines()
			for num2 in range(0,len(data)): 
				if 'devauc' in data[num2]:
					x_b = float(data[num2].split()[3].replace('*',' ')[:-1])
					y_b = float(data[num2].split()[4].replace('*',' ')[:-1])
					MAG_B = float(data[num2].split()[5].replace('*',' '))
					re = float(data[num2].split()[6].replace('*',' '))
					if data[num2+1].split()[0] == '(':
						MAG_B_err = float(data[num2+1].split()[3].replace('*',' '))
						re_err = float(data[num2+1].split()[4].replace('*',' '))
					else:
						MAG_B_err = float(data[num2+1].split()[2].replace('*',' '))
						re_err = float(data[num2+1].split()[3].replace('*',' '))
				elif 'disk' in data[num2]: 
					x_d = float(data[num2].split()[3].replace('*',' ')[:-1])
					y_d = float(data[num2].split()[4].replace('*',' ')[:-1])
					MAG_D = float(data[num2].split()[5].replace('*',' '))
					rs = float(data[num2].split()[6].replace('*',' '))
					if data[num2+1].split()[0] == '(':
						MAG_D_err = float(data[num2+1].split()[3].replace('*',' '))
						rs_err = float(data[num2+1].split()[4].replace('*',' '))
					else:
						MAG_D_err = float(data[num2+1].split()[2].replace('*',' '))
						rs_err = float(data[num2+1].split()[3].replace('*',' '))
				elif 'Chi^2/nu' in data[num2]:
					chi_nu = float(data[num2].split()[2].replace('*',' '))
			DB = 10**((MAG_B-MAG_D)/2.5)
			BT = 1/(DB+1)
			array[num]=int(gala[num]),x_b,y_b,MAG_B,MAG_B_err,re,re_err,x_d,y_d,MAG_D,MAG_D_err,rs,rs_err,DB,BT,chi_nu
print('-------------------------- ready to save FITS file --------------------------')
col1 = fits.Column(name='ID_Huang',array=array[:,0],format='K')
col2 = fits.Column(name='x_bulge',array=array[:,1],format='D')
col3 = fits.Column(name='y_bulge',array=array[:,2],format='D')
col4= fits.Column(name='MAG_bulge',array=array[:,3],format='D')
col5= fits.Column(name='MAG_bulge_err',array=array[:,4],format='D')
col6= fits.Column(name='re',array=array[:,5],format='D')
col7= fits.Column(name='re_err',array=array[:,6],format='D')
col8 = fits.Column(name='x_disk',array=array[:,7],format='D')
col9 = fits.Column(name='y_disk',array=array[:,8],format='D')
col10= fits.Column(name='MAG_disk',array=array[:,9],format='D')
col11= fits.Column(name='MAG_disk_err',array=array[:,10],format='D')
col12= fits.Column(name='rs',array=array[:,11],format='D')
col13= fits.Column(name='rs_err',array=array[:,12],format='D')
col14 = fits.Column(name='disk2bulge',array=array[:,13],format='D')
col15 = fits.Column(name='bulge2total',array=array[:,14],format='D')
col16 = fits.Column(name='chi_nu',array=array[:,15],format='D')
hdu = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16])
hdu.writeto(path+'/galfit_result.fits',overwrite=True)
print('-------------------------- galfit result has been saved --------------------------')
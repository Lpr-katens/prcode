from astropy.io import fits
import numpy as np
import os
import matplotlib.pyplot as plt
root_path = '/Users/lpr/Data/galfit/test/goodsn/'
temp = os.listdir(root_path)
galaxy_list = []
for file in temp:
	if file.endswith('.DS_Store') or file.endswith('.fits'):
		print(file)
	else:
		galaxy_list.append(file)
for file in galaxy_list:
	galaxy_id = int(file)
	os.chdir(root_path)
	os.chdir(file)
	f = open('fit.log')
	data = f.readlines()
	devau_iteration = np.empty((0,4))
	exp_iteration = np.empty((0,6))
	for num2 in range(0,len(data)): 
		if 'devauc' in data[num2]:
			x_b = data[num2].split()[3].replace('*',' ')[:-1]
			y_b = data[num2].split()[4].replace('*',' ')[:-1]
			MAG_B = data[num2].split()[5].replace('*',' ')
			re = data[num2].split()[6].replace('*',' ')
			devau_iteration = np.append(devau_iteration,np.array([[x_b,y_b,MAG_B,re]]),axis=0)
		elif 'disk' in data[num2]: 
			x_d = data[num2].split()[3].replace('*',' ')[:-1]
			y_d = data[num2].split()[4].replace('*',' ')[:-1]
			MAG_D = data[num2].split()[5].replace('*',' ')
			rs = data[num2].split()[6].replace('*',' ')
			ba = data[num2].split()[7].replace('*',' ')
			pa = data[num2].split()[8].replace('*',' ')
			exp_iteration = np.append(exp_iteration,np.array([[x_d,y_d,MAG_D,rs,ba,pa]]),axis=0)
	col1 = fits.Column(name='x_b',array=devau_iteration[:,0],format='D')
	col2 = fits.Column(name='y_b',array=devau_iteration[:,1],format='D')
	col3 = fits.Column(name='mag_b',array=devau_iteration[:,2],format='D')
	col4= fits.Column(name='re',array=devau_iteration[:,3],format='D')
	col5= fits.Column(name='x_d',array=exp_iteration[:,0],format='D')
	col6= fits.Column(name='y_d',array=exp_iteration[:,1],format='D')
	col7= fits.Column(name='mag_d',array=exp_iteration[:,2],format='D')
	col8 = fits.Column(name='rs',array=exp_iteration[:,3],format='D')
	col9 = fits.Column(name='ba',array=exp_iteration[:,4],format='D')
	col10= fits.Column(name='pa',array=exp_iteration[:,5],format='D')
	hdu = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10])
	hdu.writeto('iter_result.fits',overwrite=True)
	plt.figure(figsize=(15,10))
	plt.subplot(2,5,1)
	plt.plot(devau_iteration[:,0],label='x_b',color='olive',linewidth=2)
	plt.legend(loc=9,bbox_to_anchor=(0.5,1.1))
	plt.subplot(2,5,2)
	plt.plot(devau_iteration[:,1],label='y_b',color='green',linewidth=2)
	plt.legend(loc=9,bbox_to_anchor=(0.5,1.1))
	plt.subplot(2,5,3)
	plt.plot(devau_iteration[:,2],label='mag_b',color='cyan',linewidth=2)
	plt.legend(loc=9,bbox_to_anchor=(0.5,1.1))
	plt.subplot(2,5,4)
	plt.plot(devau_iteration[:,3],label='re',color='blue',linewidth=2)
	plt.legend(loc=9,bbox_to_anchor=(0.5,1.1))
	plt.subplot(2,5,5)
	plt.plot(exp_iteration[:,0],label='x_d',color='black',linewidth=2)
	plt.legend(loc=9,bbox_to_anchor=(0.5,1.1))
	plt.subplot(2,5,6)
	plt.plot(exp_iteration[:,1],label='y_d',color='brown',linewidth=2)
	plt.legend(loc=9,bbox_to_anchor=(0.5,1.1))
	plt.subplot(2,5,7)
	plt.plot(exp_iteration[:,2],label='mag_d',color='red',linewidth=2)
	plt.legend(loc=9,bbox_to_anchor=(0.5,1.1))
	plt.subplot(2,5,8)
	plt.plot(exp_iteration[:,3],label='rs',color='orange',linewidth=2)
	plt.legend(loc=9,bbox_to_anchor=(0.5,1.1))
	plt.subplot(2,5,9)
	plt.plot(exp_iteration[:,4],label='ba',color='grey',linewidth=2)
	plt.legend(loc=9,bbox_to_anchor=(0.5,1.1))
	plt.subplot(2,5,10)
	plt.plot(exp_iteration[:,5],label='pa',color='pink',linewidth=2)
	plt.legend(loc=9,bbox_to_anchor=(0.5,1.1))
	plt.suptitle(file)
	plt.savefig('iter_result.eps')
	plt.close()
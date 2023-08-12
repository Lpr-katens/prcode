from astropy.io import fits
from photutils.aperture import CircularAperture, aperture_photometry
import photutils
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from lpr.image import display
import subprocess as sp
import os

ctg_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
img_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
fields = ['goodsn']#,'goodss','egs'
output_path = '/Users/lpr/Data/lirg_project/output/van_re_half_light_radius/'
detection_threshold = 3
npixels = 21
pixel_scale = 0.06
# 定义一个圆，用来画re的位置
p = np.linspace(0,1,100)*(2*np.pi)
segm_path = '/Users/lpr/Data/lirg_project/output/van_re_half_light_radius/segmentation/'

for field in fields:
	ctg = fits.getdata(ctg_path+field+'_Huangall_candels_radec_van_modifyz.fits',1)
	idx_used = fits.getdata(ctg_path+field+'_id.fits',1)
	idx_list = np.full(len(idx_used),-99)
	array = np.full([len(idx_used),2],-99.)
	for num in range(0,len(idx_used)):#
		idx = idx_used[num]['id']
		if idx == 32411:
			# ra = ctg[ctg['id']==idx]['ra_candels'][0]
			# dec = ctg[ctg['id']==idx]['dec_candels'][0]
			# re = ctg[ctg['id']==idx]['re_f160w'][0]
			# flag = ctg[ctg['id']==idx]['f_f160w']
			# # # 如果出现图像不够大的情况，需要重新cut图像，做半光半径的话
			# # img = fits.getdata('/Users/lpr/Data/lirg_project/output/van_re_half_light_radius/goodsn_f160w_32411.fits',0)
			# # hdr = fits.getheader('/Users/lpr/Data/lirg_project/output/van_re_half_light_radius/goodsn_f160w_32411.fits',0)
			# # 如果图像大小足够
			# img = fits.getdata(img_path+field+'/'+field+'_f160w/'+field+'_f160w_'+str(idx)+'.fits',0)
			# hdr = fits.getheader(img_path+field+'/'+field+'_f160w/'+field+'_f160w_'+str(idx)+'.fits',0)
			# x,y = WCS(hdr).wcs_world2pix(ra,dec,0)
			# radii = np.linspace(1,50,50)#半径
			# fs = []#半径内的流量
			# # 给图像算背景
			# noise = np.median(img[(img>-0.01)&(img<0.01)].flatten())
			# noise_sigma = np.std(img[(img>-0.01)&(img<0.01)].flatten())
			# noise_map = np.random.normal(noise,noise_sigma,img.shape)#
			# 做sextractor
			os.chdir(segm_path)
			print('yes')
			sp.run('sex '+img_path+field+'/'+field+'_f160w/'+field+'_f160w_'+str(idx)+'.fits'+' -c default.sex -CATALOG_NAME '+str(idx)+'.cat -CHECKIMAGE_NAME '+str(idx)+'.fits -DETECT_MINAREA 21',shell=True,check=True)# -DETECT_THRESH 1 -ANALYSIS_THRESH 1 -DETECT_THRESH 5
			# # 如果出现图像不够大的情况，需要重新cut图像，做半光半径的话
			# sp.run('sex /Users/lpr/Data/lirg_project/output/van_re_half_light_radius/goodsn_f160w_32411.fits'+' -c default.sex -CATALOG_NAME '+str(idx)+'.cat -CHECKIMAGE_NAME '+str(idx)+'.fits',shell=True,check=True)
	# 		segmap = fits.getdata(str(idx)+'.fits',0)
	# 		img[(segmap != segmap[int((img.shape[0]-1)/2),int((img.shape[1]-1)/2)])&(segmap != 0)] = noise_map[(segmap != segmap[int((img.shape[0]-1)/2),int((img.shape[1]-1)/2)])&(segmap != 0)]
	# 		img = img - noise
	# 		apers = [CircularAperture([x,y], r = r) for r in radii]
	# 		phot_table = aperture_photometry(img, apers)
	# 		print('begin to plot')
	# 		fig = plt.figure(figsize=[10,5])
	# 		gs = fig.add_gridspec(1,2)
	# 		fig.suptitle(field+' '+str(idx))
	# 		# 开始画生长曲线
	# 		ax = fig.add_subplot(gs[0,0])
	# 		for r in radii:
	# 			fs.append(phot_table['aperture_sum_'+str(int(r-1))][0])
	# 			ax.scatter(r,phot_table['aperture_sum_'+str(int(r-1))],color='red')
	# 		f = np.median(fs[len(fs)-9:len(fs)])
	# 		inter_func = interpolate.interp1d(fs, radii) #流量和半径的插值函数
	# 		hlr = inter_func(0.5*f) #给一半流量的地方插值出来半径
	# 		print('get half light radius')
	# 		ax.plot([hlr,hlr],[0,0.5*f],color='black')
	# 		ax.plot([0,hlr],[0.5*f,0.5*f],color='black')
	# 		ax.set_xlabel('pixels')
	# 		ax.set_ylabel('arbitrary unit')
	# 		# 开始画星系的图像
	# 		ax = fig.add_subplot(gs[0,1])
	# 		ax.imshow(display.logstretch(img,a=500),cmap='binary')
	# 		ax.plot(x+(re/pixel_scale)*np.cos(p),y+(re/pixel_scale)*np.sin(p),color='green',label='van')
	# 		ax.plot(x+hlr*np.cos(p),y+hlr*np.sin(p),color='red',label='half light radius')
	# 		ax.set_xticks([])
	# 		ax.set_yticks([])
	# 		ax.legend()
	# 		ax.text(1,3,'re='+str(np.around(re,2))+'" '+'hlr='+str(np.around(hlr*pixel_scale,2))+'"',color='red')
	# 		plt.savefig(output_path+field+'_'+str(flag)+'_'+str(idx)+'.png')
	# 		array[num] = [re,hlr*pixel_scale] #hlr是pixel的单位，re是角秒的单位，pixel scale是0.06
	# 		idx_list[num] = idx
	# 		print(re,hlr*pixel_scale)
	# # # col1 = fits.Column(name='id',array=idx_list,format='K')
	# # # col2 = fits.Column(name='re',array=array[:,0],format='D')
	# # # col3 = fits.Column(name='hlr',array=array[:,1],format='D')
	# # # hdu = fits.BinTableHDU.from_columns([col1,col2,col3])
	# # # hdu.writeto(output_path+field+'_re_hlr.fits',overwrite=True)
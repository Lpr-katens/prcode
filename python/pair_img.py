# 把有pair的星系在图上标出来，并且标出来它们的pair的位置，spec的pair用红色，phot的pair用蓝色
from astropy.io import fits
import numpy as np
from lpr.image import display
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from shutil import copy
import os

ctg_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
c_ctg_path = '/Users/lpr/Data/lirg_project/intake/CANDELS/catalog/JFang_CANDELS_Data/'
img_path = '/Users/lpr/Data/lirg_project/output/pair/'
fields = ['goodsn','goodss','egs']#
c_ctg_name = {'goodsn':'gdn','goodss':'gds','egs':'egs'}
flag_color = ['red','blue']

for field in fields:
	p_ctg = fits.getdata(ctg_path+field+'_pair_1000kms.fits',1) #储存pair信息的表
	h_ctg = fits.getdata(ctg_path+field+'_Huangall_radec_candels.fits',1) #储存样本与candels匹配信息的表
	c_ctg = fits.getdata(c_ctg_path+c_ctg_name[field]+'_all.fits',1) #candels表
	for num in range(0,len(p_ctg)):
		id_h = p_ctg[num]['id'] #目标星系
		f_h = p_ctg[num]['obj_z_flag']
		if p_ctg[num]['count_spec']+p_ctg[num]['count_phot'] != 0: #只对有pair的运行这个程序
			print(p_ctg[num]['count_spec']+p_ctg[num]['count_phot'])
			img = fits.getdata(img_path+field+'_1000kms/'+field+'_f160w/'+field+'_f160w_'+str(id_h)+'.fits',0)
			hdr = fits.getheader(img_path+field+'_1000kms/'+field+'_f160w/'+field+'_f160w_'+str(id_h)+'.fits',0)
			hdr['CTYPE1'] = 'RA---TAN-SIP'
			hdr['CTYPE2'] = 'DEC--TAN-SIP'
			ra_h,dec_h = h_ctg[h_ctg['id']==id_h]['ra_candels'],h_ctg[h_ctg['id']==id_h]['dec_candels']
			x_h,y_h = WCS(hdr).wcs_world2pix(ra_h,dec_h,0)
			plt.figure(figsize=[8,8])
			plt.imshow(display.logstretch(img,a=1000),cmap='binary')
			plt.text(5,5,field+' '+str(id_h),fontsize=15,color='red')
			plt.scatter(x_h,y_h,color=flag_color[f_h-1],marker='o')
			pair_id_list = p_ctg[num]['com_id_can']
			pair_z_flag = p_ctg[num]['com_z_flag']
			for num1 in range(0,len(pair_id_list)):
				id_p = pair_id_list[num1]
				if id_p != -99:
					ra_p,dec_p = c_ctg[c_ctg['id']==id_p]['ra_1'],c_ctg[c_ctg['id']==id_p]['dec_1']
					x_p,y_p = WCS(hdr).wcs_world2pix(ra_p,dec_p,0)
					f_p = pair_z_flag[num1]
					plt.scatter(x_p,y_p,color=flag_color[f_p-1],marker='x')
			plt.savefig(img_path+'1000kms_'+field+'_'+str(id_h)+'.pdf')

hdu = fits.getdata(img_path+'pair_in_16sample.fits',1)
for idx in hdu['id']:
	for field in fields:
		if '1000kms_'+field+'_'+str(idx)+'.pdf' in os.listdir(img_path):
			copy(img_path+'1000kms_'+field+'_'+str(idx)+'.pdf',img_path+'pair_in_16sample/'+'1000kms_'+field+'_'+str(idx)+'.pdf')
from astropy.io import fits
import numpy as np

# 根据之前的粗略选择（只满足separation），在判断如果两者都有光谱红移，那么要满足500km/s(deltav = c*deltaz/(1+z))，若果是测光红移，那么就要满足delta z / (1 + z) < 0.1且小于2sigma,sigma的值GOODS场为0.031。还要满足质量比1：10.

fields = ['goodsn','goodss','egs']
path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
candels_path = '/Users/lpr/Data/lirg_project/intake/CANDELS/catalog/JFang_CANDELS_Data/'
fang_ctg = {'goodsn':'gdn','goodss':'gds','egs':'egs'}

from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from astropy.cosmology import FlatLambdaCDM as flcdm
import time

st = time.process_time()
sigma_photz = 0.031 #GOODS的photoz的sigma的中值，来自dalen+12
velocity_specz =1.e6  #1000km/s()或者500km/s(5.e5)的光谱红移的差别
for field in fields:
	com_ctg = fits.getdata(path+field+'_companion.fits',1)
	h_ctg = fits.getdata(path+field+'_Huangall_candels_van_params_zmasshuangsfr.fits',1)
	id_ctg = fits.getdata(path+field+'_id.fits',1)
	h_ctg = h_ctg[np.isin(h_ctg['id'],id_ctg['id'])]
	c_ctg = fits.getdata(candels_path+fang_ctg[field]+'_all.fits',1)
	idx_list = np.full(len(com_ctg),-99) #用来装样本的id的
	companion_list = np.full((len(com_ctg),10),-99) #用来装可能的companion的id的
	companion_z_list = np.full((len(com_ctg),10),-99.) ##用来装可能的companion的红移的
	object_z_list = np.full(len(com_ctg),-99.)
	object_z_flag_list = np.full(len(com_ctg),-99) #用来装样本的红移是否是光谱红移的
	companion_z_flag_list = np.full((len(com_ctg),10),-99.) #用来装可能的companion的红移是否是光谱红移的
	spec_pair_count_list = np.full(len(com_ctg),-99) #用来统计有多少光谱认证的pair
	phot_pair_count_list = np.full(len(com_ctg),-99) #用来统计有多少满足测光红移条件的潜在pair
	for num in range(0,len(com_ctg)):
		id_h = com_ctg[num]['id']
		f_h = com_ctg[num]['obj_z_flag']
		z_h = com_ctg[num]['object_z']
		mass_h = h_ctg[h_ctg['id']==id_h]['lmass_candels']
		print('object mass: ' + str(mass_h))
		count_spec = 0
		count_phot = 0
		for num1 in range(0,len(com_ctg[num]['com_id_can'])):
			id_com = com_ctg[num]['com_id_can'][num1]
			if id_com != -99:
				f_com = com_ctg[num]['com_z_flag'][num1]
				z_com = com_ctg[num]['com_z_can'][num1]
				mass_com = c_ctg[c_ctg['id']==id_com]['m_med']
				print('companion mass: ' + str(mass_com))
				if f_h == 1 and f_com == 1 and abs(z_h - z_com) <= velocity_specz * (1 + z_h) / 3.e8 and abs(mass_h - mass_com) <= 1: #满足红移要求和质量比
					print('satisfy mass limit')
					companion_list[num,num1] = id_com
					companion_z_list[num,num1] = z_com
					companion_z_flag_list[num,num1] = f_com
					count_spec += 1
				elif f_h == 2 or f_com == 2:
					if abs(mass_h - mass_com) <= 1 and abs(z_h - z_com) <= 0.1 * (1 + z_h) and abs(z_h - z_com) <= sigma_photz*2:
						companion_list[num,num1] = id_com
						companion_z_list[num,num1] = z_com
						companion_z_flag_list[num,num1] = f_com
						count_phot += 1
		idx_list[num] = id_h
		object_z_list[num] = z_h
		object_z_flag_list[num] = f_h
		spec_pair_count_list[num] = count_spec
		phot_pair_count_list[num] = count_phot
	# 开始用数据生成fits表格
	col0 = fits.Column(name='id',format='K',array=idx_list)
	col1 = fits.Column(name='obj_z_flag',format='K',array=object_z_flag_list)
	col2 = fits.Column(name='object_z',format='D',array=object_z_list)
	col3 = fits.Column(name='com_id_can',format='PK',array=companion_list)
	col4 = fits.Column(name='com_z_can',format='PD',array=companion_z_list)
	col5 = fits.Column(name='com_z_flag',format='PK',array=companion_z_flag_list)
	col6 = fits.Column(name='count_spec',format='K',array=spec_pair_count_list)
	col7 = fits.Column(name='count_phot',format='K',array=phot_pair_count_list)
	hdu = fits.BinTableHDU.from_columns([col0,col1,col2,col3,col4,col5,col6,col7])
	hdu.writeto('/Users/lpr/Data/lirg_project/output/catalog_radec/'+field+'_pair_1000kms.fits',overwrite=True)
	print('---------- ' + field + ' ----------')
et = time.process_time()
print('time consume: ' + str(et - st))
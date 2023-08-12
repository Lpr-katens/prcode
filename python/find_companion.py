# 这个程序是文章中第三节，用来找我们样本潜在的companion，条件是满足0.5<z<1.6，质量大于我们样本的最小质量的十分之一。
# 根据上面的条件把大表裁剪成小的表，然后对于我们的样本中的每一个，在天区上画一个方形的区域，ra-0.009：ra+0.009，dec-0.005：的车+0.005。
# 这个方形的区域就是，最后需要用for循环来找满足separation的最终小表。
# 然后，对于这个小表中的每一个源，计算它们和目标源的separation，满足30/(h kpc)的，而且满足质量比大于0.1的都算。
# 之后选什么光谱红移，测光红移之类的再说。

fields = ['goodsn','goodss','egs']#
huang_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
candels_path = '/Users/lpr/Data/lirg_project/intake/CANDELS/catalog/JFang_CANDELS_Data/'
fang_ctg = {'goodsn':'gdn','goodss':'gds','egs':'egs'}

from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from astropy.cosmology import FlatLambdaCDM as flcdm
import time

st = time.process_time()
def kpc_per_arcsec(z):
    angular_distance = flcdm(H0=70,Om0=0.3).angular_diameter_distance(z).to_value()
    arc_scale = angular_distance*np.pi*1000/(180*3600)
    return(arc_scale)

separation = 30/0.7 #0.7是h,这个参数是从newman+12中拿出来的,unit:kpc
distance = separation/kpc_per_arcsec(1) #红移1的时候，这个separation对应的角距离,大概是5个角秒
# delta RA: 0.009大概是15个角秒; delta DEC: 0.005大概是18个角秒

for field in fields:
	h_ctg = fits.getdata(huang_path+field+'_Huangall_candels_van_params_zmasshuangsfr.fits',1)
	id_ctg = fits.getdata(huang_path+field+'_id.fits',1)
	h_ctg = h_ctg[np.isin(h_ctg['id'],id_ctg['id'])]
	c_ctg = fits.getdata(candels_path+fang_ctg[field]+'_all.fits',1)
	c_ctg = c_ctg[(c_ctg['m_med'] > 6.86)&(c_ctg['zbest']>=0.5)&(c_ctg['zbest']<=1.6)&(c_ctg['class_star']<0.8)] #先限制质量比和红移来缩小数据量，我们数据中的最小质量是7.68
	idx_list = np.full(len(h_ctg),-99) #用来装样本的id的
	companion_list = np.full((len(h_ctg),10),-99) #用来装可能的companion的id的
	companion_z_list = np.full((len(h_ctg),10),-99.) ##用来装可能的companion的红移的
	object_z_list = np.full(len(h_ctg),-99.)
	object_z_flag_list = np.full(len(h_ctg),-99) #用来装样本的红移是否是光谱红移的
	companion_z_flag_list = np.full((len(h_ctg),10),-99.) #用来装可能的companion的红移是否是光谱红移的
	for num in range(0,len(h_ctg)):
		c_ctg_copy = c_ctg
		idx = h_ctg[num]['id']
		print('------' + field + '------' + str(idx) + '------')
		idx_h_can = h_ctg[num]['id_candels']
		ra,dec = h_ctg[num]['ra_candels'],h_ctg[num]['dec_candels']
		z_h = h_ctg[num]['z_used'] #huang表中的红移
		mass_h = h_ctg[num]['lmass_candels'] #huang表中的质量
		if h_ctg[num]['zspec1'] != 0. and h_ctg[num]['zspec1'] != -99.:
			z_h_flag = 1 #红移的flag：1是有光谱红移的，2是没有光谱红移的
		else:
			z_h_flag = 2
		p1 = SkyCoord(ra*u.deg,dec*u.deg,frame='fk5')
		com_idx_list = [] #临时存储companion的id
		com_z_list = [] #临时存储companion的红移
		com_z_flag_list = [] #临时存储companion的红移的flag
		c_ctg_copy = c_ctg_copy[(c_ctg_copy['ra_1'] > ra - 0.009) & (c_ctg_copy['ra_1'] < ra + 0.009) & (c_ctg_copy['dec_1'] > dec - 0.005) & (c_ctg_copy['dec_1'] < dec + 0.005)] #先把表裁成小表
		count = 0 #用来对满足separation条件的companion计数
		for num1 in range(0,len(c_ctg_copy)):
			idx_can = c_ctg_copy[num1]['id']
			z = c_ctg_copy[num1]['zbest']
			mass = c_ctg_copy[num1]['m_med']
			if c_ctg_copy[num1]['spec_z'] != 0. and c_ctg_copy[num1]['spec_z'] != -99.:
				z_c_flag = 1 #红移的flag：1是有光谱红移的，2是没有光谱红移的
			else:
				z_c_flag = 2 
			ra_can,dec_can = c_ctg_copy[num1]['ra_1'],c_ctg_copy[num1]['dec_1']
			p2 = SkyCoord(ra_can*u.deg,dec_can*u.deg,frame='fk5')
			sep = p1.separation(p2).arcsec
			if sep <= distance and idx_h_can != idx_can and abs(mass/mass_h) >= 0.1:
				print('--------' + str(idx) + '-------- candels companion' + str(idx_can) + '------')
				companion_list[num,count] = idx_can
				companion_z_list[num,count] = z
				companion_z_flag_list[num,count] = z_c_flag
				count += 1

		idx_list[num] = idx
		object_z_list[num] = z_h
		object_z_flag_list[num] = z_h_flag
	
	# 开始用数据生成fits表格
	col0 = fits.Column(name='id',format='K',array=idx_list)
	col1 = fits.Column(name='obj_z_flag',format='K',array=object_z_flag_list)
	col2 = fits.Column(name='object_z',format='D',array=object_z_list)
	col3 = fits.Column(name='com_id_can',format='PK',array=companion_list)
	col4 = fits.Column(name='com_z_can',format='PD',array=companion_z_list)
	col5 = fits.Column(name='com_z_flag',format='PK',array=companion_z_flag_list)
	hdu = fits.BinTableHDU.from_columns([col0,col1,col2,col3,col4,col5])
	hdu.writeto('/Users/lpr/Data/lirg_project/output/catalog_radec/'+field+'_companion.fits',overwrite=True)
et = time.process_time()
print('time consume: ' + str(et - st))
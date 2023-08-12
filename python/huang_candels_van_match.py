# 把huang21文章和CANDELS数据做坐标的匹配

fields = ['goodsn','goodss','egs']#
huang_path = '/Users/lpr/Data/lirg_project/intake/huang_catalog/'
candels_path = '/Users/lpr/Data/lirg_project/intake/CANDELS/catalog/JFang_CANDELS_Data/'
fang_ctg = {'goodsn':'gdn','goodss':'gds','egs':'egs'}
van_path = '/Users/lpr/Data/lirg_project/intake/3dhst/3dhst/'
van_ctg_name = {'goodsn':'goodsn','goodss':'goodss','egs':'aegis'}
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from tqdm import tqdm

# # delta RA: 0.009大概是15个角秒; delta DEC: 0.005大概是18个角秒
# for field in fields:
# 	h_ctg = fits.open(huang_path+field+'_f16_zeq1.cat')[1].data
# 	h_ctg = h_ctg[h_ctg['id']!=-1]
# 	c_ctg = fits.open(candels_path+fang_ctg[field]+'_all.fits')[1].data
# 	idx_list = np.full((len(h_ctg),2),-99) #用来装样本的id的
# 	separation_list = np.full(len(h_ctg),-99.)
# 	coordcan_list = np.full((len(h_ctg),2),-99.)
# 	for num in tqdm(range(len(h_ctg))):
# 		idx = h_ctg[num]['id']
# 		ra,dec = h_ctg[num]['ra'],h_ctg[num]['dec']
# 		p1 = SkyCoord(ra*u.deg,dec*u.deg,frame='fk5')
# 		# 先把表裁成小表
# 		c_ctg_copy = c_ctg[(c_ctg['ra_1'] > ra - 0.009) & (c_ctg['ra_1'] < ra + 0.009) & (c_ctg['dec_1'] > dec - 0.005) & (c_ctg['dec_1'] < dec + 0.005)]
# 		ra_can = c_ctg_copy['ra_1']
# 		dec_can = c_ctg_copy['dec_1']
# 		p2 = SkyCoord(ra_can*u.deg,dec_can*u.deg,frame='fk5')
# 		sep = p2.separation(p1).arcsec
# 		try:
# 			idx_can = c_ctg_copy[np.argmin(sep)]['id']
# 			idx_list[num] = [idx,idx_can]
# 			separation_list[num] = sep[np.argmin(sep)]
# 			coordcan_list[num] = [ra_can[np.argmin(sep)],dec_can[np.argmin(sep)]]
# 		except:
# 			continue
# 	# 开始用数据生成fits表格
# 	col0 = fits.Column(name='id',format='K',array=idx_list[:,0])
# 	col1 = fits.Column(name='id_candels',format='K',array=idx_list[:,1])
# 	col2 = fits.Column(name='separation_16_candels',format='D',array=separation_list)
# 	col3 = fits.Column(name='ra_candels',format='D',array=coordcan_list[:,0])
# 	col4 = fits.Column(name='dec_candels',format='D',array=coordcan_list[:,1])
# 	hdu = fits.BinTableHDU.from_columns([col0,col1,col2,col3,col4])
# 	hdu.writeto('/Users/lpr/Data/lirg_project/output/catalog_radec/'+field+'_match_candels.fits',overwrite=True)

for field in fields:
	h_ctg = fits.getdata('/Users/lpr/Data/lirg_project/output/catalog_radec/'+field+'_match_candels.fits',1)
	h_ctg = h_ctg[(h_ctg['id']!=-1)&(h_ctg['id']!=-99)]
	van_ctg = fits.getdata(van_path+van_ctg_name[field]+'/'+van_ctg_name[field]+'_3dhst.v4.1_f160wf125w.fits',1)
	idx_list = np.full((len(h_ctg),3),-99) #用来装样本的id的
	separation_list = np.full((len(h_ctg),2),-99.)
	coord_list = np.full((len(h_ctg),4),-99.)
	for num in tqdm(range(len(h_ctg))):
		idx = h_ctg[num]['id']
		idx_can = h_ctg[num]['id_candels']
		sep_16_can = h_ctg[num]['separation_16_candels']
		ra,dec = h_ctg[num]['ra_candels'],h_ctg[num]['dec_candels']
		p1 = SkyCoord(ra*u.deg,dec*u.deg,frame='fk5')
		# 先把表裁成小表
		van_ctg_copy = van_ctg[(van_ctg['ra_van'] > ra - 0.009) & (van_ctg['ra_van'] < ra + 0.009) & (van_ctg['dec_van'] > dec - 0.005) & (van_ctg['dec_van'] < dec + 0.005)]
		ra_van = van_ctg_copy['ra_van']
		dec_van = van_ctg_copy['dec_van']
		p2 = SkyCoord(ra_van*u.deg,dec_van*u.deg,frame='fk5')
		sep = p1.separation(p2).arcsec
		try:

			idx_van = van_ctg_copy[np.argmin(sep)]['id_van']
			idx_list[num] = [idx,idx_can,idx_van]
			separation_list[num] = [sep_16_can,sep[np.argmin(sep)]]
			coord_list[num] = [ra,dec,ra_van[np.argmin(sep)],dec_van[np.argmin(sep)]]
		except:
			continue
	# 开始用数据生成fits表格
	col0 = fits.Column(name='id',format='K',array=idx_list[:,0])
	col1 = fits.Column(name='id_candels',format='K',array=idx_list[:,1])
	col2 = fits.Column(name='id_van',format='K',array=idx_list[:,2])
	col3 = fits.Column(name='separation_16_candels',format='D',array=separation_list[:,0])
	col4 = fits.Column(name='separation_candels_van',format='D',array=separation_list[:,1])
	col5 = fits.Column(name='ra_candels',format='D',array=coord_list[:,0])
	col6  = fits.Column(name='dec_candels',format='D',array=coord_list[:,1])
	col7 = fits.Column(name='ra_van',format='D',array=coord_list[:,2])
	col8  = fits.Column(name='dec_van',format='D',array=coord_list[:,3])
	hdu = fits.BinTableHDU.from_columns([col0,col1,col2,col3,col4,col5,col6,col7,col8])
	hdu.writeto('/Users/lpr/Data/lirg_project/output/catalog_radec/'+field+'_match_candels_van.fits',overwrite=True)	
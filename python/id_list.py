from astropy.io import fits
import numpy as np
from astropy.table import join

fields = ['goodsn','goodss','egs']#
# for field in fields:
# 	hdu = fits.getdata('/Users/lpr/Data/lirg_project/output/catalog_radec/'+field+'_match_candels_van.fits',1)
# 	huang = fits.getdata('/Users/lpr/Data/lirg_project/intake/huang_catalog/'+field+'_Huang_all.fits',1)
# 	huang = huang[np.isin(huang['id'],hdu['id'])]
# 	data = join(hdu,huang,keys_left='id',keys_right='ID')
# 	data.write('/Users/lpr/Data/lirg_project/output/catalog_radec/'+field+'_Huangall_candels_van.fits',overwrite=True)
# van_path = '/Users/lpr/Data/lirg_project/intake/3dhst/3dhst/'
# van_ctg_name = {'goodsn':'goodsn','goodss':'goodss','egs':'aegis'}
# for field in fields:
# 	hdu = fits.getdata('/Users/lpr/Data/lirg_project/output/catalog_radec/'+field+'_Huangall_candels_van.fits',1)
# 	van = fits.getdata(van_path+van_ctg_name[field]+'/'+van_ctg_name[field]+'_3dhst.v4.1_f160wf125w.fits',1)
# 	van = van[np.isin(van['id_van'],hdu['id_van'])]
# 	data = join(hdu,van,keys_left='id_van',keys_right='ID_van')
# 	data.write('/Users/lpr/Data/lirg_project/output/catalog_radec/'+field+'_Huangall_candels_van_params.fits',overwrite=True)
del_id = {'egs':[13011807,12101077,-1,13025385,13011815],'goodss':[2324,7748,16929,17215,15278,16731,16003,14856,14845,13663,13557,13121,12953,12814,12808,12467,11380,12237],'goodsn':[36821,31002,26981,11380,12237,19196,34172,42252,13027]}
ctg_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
for field in fields:
	ctg = fits.getdata(ctg_path+field+'_Huangall_candels_van_params_modifyz.fits',1)
	idx_list = []
	for num in range(0,len(ctg)):
		idx = ctg[num]['id']
		separation = ctg[num]['Separation_16_candels']
		z = ctg[num]['z_used']
		if idx not in del_id[field] and separation <= 1. and z <= 1.3 and z >= 0.8:
			idx_list.append(idx)
	idx_list = np.array(idx_list)
	col = fits.Column(name='id',array=idx_list,format='K')
	hdu = fits.BinTableHDU.from_columns([col])
	hdu.writeto(ctg_path+field+'_id.fits',overwrite=True)
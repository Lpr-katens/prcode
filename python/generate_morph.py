# 生成fits表格
from astropy.io import fits
import numpy as np

fields = ['goodsn','goodss','egs']
root_path = '/Users/lpr/Data/lirg_project/output/'
for field in fields:
	morph_ctg = fits.getdata(root_path+'catalog/'+field+"_Huangall_morph.fits",1)
	id_ctg = fits.getdata(root_path+'catalog_radec/'+field+'_id.fits',1)
	info_list = np.full([len(id_ctg),2],-99)
	for num in range(0,len(id_ctg)):
		idx = id_ctg[num]['id']
		# spiral/disky=1,irregular=2,compact=3,two components=4
		morph = morph_ctg[morph_ctg['id']==idx]['morph'][0]
		info_list[num] = [idx,morph]
	col1 = fits.Column(name='id',array=info_list[:,0],format='K')
	col2 = fits.Column(name='morph',array=info_list[:,1],format='K')
	hdu = fits.BinTableHDU.from_columns([col1,col2])
	hdu.writeto(root_path+'catalog_radec/'+field+'_morph.fits',overwrite=True)
	print('                       ***\n'
		'                      **\n'
		'                     **\n'
		'                    **\n'
		'                   **\n'
		'                  **\n'
		'                 **\n'
		'                **\n'
		'               **\n'
		' ***          **\n'
		'   **        **\n'
		'    **      **\n'
		'     **    **\n'
		'      **  **\n'
		'       ****\n'
		'        **\n')
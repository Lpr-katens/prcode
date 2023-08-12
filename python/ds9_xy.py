from astropy.io import fits
from astropy.wcs import WCS
ctg = fits.getdata('/Users/lpr/Data/lirg_project/intake/CANDELS/catalog/JFang_CANDELS_Data/egs_all.fits',1)
img,hdr = fits.getdata('/Users/lpr/Data/lirg_project/intake/CANDELS/egs/aegis_3dhst.v4.0.F160W_orig_sci.fits',0,header=True)
idx_list = ctg['id']
ra,dec = ctg['ra_1'],ctg['dec_1']
x,y = WCS(hdr).wcs_world2pix(ra,dec,0)
col1=fits.Column(name='id',array=idx_list,format='K')
col2=fits.Column(name='x',array=x,format='D')
col3=fits.Column(name='y',array=y,format='D')
hdu=fits.BinTableHDU.from_columns([col1,col2,col3])
outputpath = input('Input the path of the output catalog: ' + '/')
outputname = input('Input the name of the output catalog: ')
hdu.writeto(outputpath+'/'+outputname,overwrite=True)
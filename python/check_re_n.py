from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

catalog = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/goodsn_all.fits')[1].data
image_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/singleGALAXYimage/gdn_f160w'
export_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/re_n_check/1'
# circle function
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
radius = 41 # 2.5"
x = radius * np.cos( angle )+50
y = radius * np.sin( angle )+50

idx = np.array([13918,9988,2929,24222,17749,26091,10437,26797,18654,27048,1323,11519,9094,2027,6403,17020,16884,1308,12288,22407,11010,14408,26188,38311,21468,15742,34371,35332,34537,27302,25645,37446,33272,34501,24080,18711,26926])
for num in range(0,len(catalog)):
	galaxyid = catalog[num]['ID_Huang']
	redshift = catalog[num]['REDSHIFT']
	re = catalog[num]['re_f160w']
	dre = catalog[num]['dre_f160w']
	n = catalog[num]['n_f160w']
	dn = catalog[num]['dn_f160w']
	flag = catalog[num]['f']
	if galaxyid in idx:
		image = np.log10(fits.open(image_path+'/gdn_f160w_'+str(galaxyid)+'.fits')[0].data)
		# image = fits.open(image_path+'/gdn_f160w_'+str(galaxyid)+'.fits')[0].data
		plt.imshow(image,interpolation=None,origin='lower',cmap='Greys')
		plt.text(0,95,'z='+str(np.around(redshift,2)),fontsize=15,color='red')
		plt.text(70,95,'re='+str(np.around(re,2)),fontsize=15,color='red')
		plt.text(70,85,'dre='+str(np.around(dre,2)),fontsize=15,color='red')
		plt.text(0,15,'n='+str(np.around(n,2)),fontsize=15,color='red')
		plt.text(0,5,'dn='+str(np.around(dn,2)),fontsize=15,color='red')
		plt.text(70,5,'flag: '+str(flag),fontsize=15,color='red')
		plt.plot(x,y,color='red')
		plt.savefig(export_path+'/'+str(galaxyid)+'_linear.eps')
		plt.close()
		print(str(galaxyid)+' done')
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

ctg_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
img_path = '/Users/lpr/Data/lirg_project/output/images/'
merger_path = '/Users/lpr/Data/lirg_project/output/gini_m20_segmap/'
fields = ['goodsn']#,'goodss','egs'
pair_flag = {True:'red',False:'blue'}
x = np.linspace(-3,0,100)
y1 = -0.14*x+0.33
y2 = 0.14*x+0.8
count_specpair = []

# regular_list = {'goodsn':[15112,13404,15995,15722,16115,13918,10534,24051,9988,11762,10685,13054,5383,14258,9457,22595,6981,5788,5271,2929,5075,20279,29082,10738,16441,14130,17749,11003,6248,24758,26797,13814,5297,29349,14780,30540,12762,16957,18654,27912,6513,3668,5858,1323,17308,1209,26506,1111,28294,32419,32851,27324,32014,34293,1114,12538,33568,12288,22758,3631,35739,29451,34043,23486,30995,13012,26118,13130,10645,26024,2142,28448,14887,25856,36636,32849,12336,31462,32411,9720,24141,14480,10776,23076,15742,14500,34371,29304,31417,38001,15758,25587,23844,19436,6796,25645,35985,27302,14565,23316,24082,28688,30235,23971,24705,30322,26926,33021,18711,32963,33322,15477,29835,24080,23226,16767,17565,33272,34501,24035,30265,13906],'goodss':[1928,2767,3452,4490,3919,26981,3631,6721,6825,8519,8326,8049,7223,10680,11793,11657,11957,12040,12993,12780,12353,13193,12997,13326,13342,14707],'egs':[13010519,13101979,13017707,13033913,13026220,13018656,13003614,13026194,13042381,12024267,13010765,12101109,13004254,13034520,13034076,13011067,13042029,13004240,13026626,12028644,12028644,12028094,12028092,13011748,13026391,13034272,12028563,13034019,13026825,13019229,13026703,13026797,12028490,13035131,13026418,13034617,13011723]}
for field in fields:
	ginim20_ctg = fits.getdata(ctg_path+field+'_GiniM20_10arcseccutout.fits',1)
	h_ctg = fits.getdata(ctg_path+field+'_Huangall_candels_radec_van_zmasssfr8sfrhuangsfr.fits',1)
	# pair_ctg = fits.getdata(ctg_path+field+'_pair_1000kms.fits',1)
	plt.figure(figsize=[5,5])
	plt.xlabel('$M_{20}$',fontsize=15)
	plt.ylabel('Gini',fontsize=15)
	plt.plot(x,y1,color='black')
	plt.plot(x,y2,color='black')
	for num in range(0,len(ginim20_ctg)):
		idx = ginim20_ctg[num]['id']
		if idx != -99:
			gini = ginim20_ctg[num]['gini_coeff_statmorph']
			m20 = ginim20_ctg[num]['moment_20_statmorph']
			# count_spec_pair = pair_ctg[pair_ctg['id']==idx]['count_spec'][0]
			# sersic = h_ctg[h_ctg['id']==idx]['n_f160w'][0]
			sfr = h_ctg[h_ctg['id']==idx]['sfr_huang'][0]
			# plt.scatter(m20,gini,color=pair_flag[count_spec_pair>0],s=5)
			# plt.scatter(m20,gini,c=sersic,cmap='jet',vmin=0.2,vmax=4,s=5)
			plt.scatter(m20,gini,c=sfr,cmap='jet',vmin=0,vmax=100,s=5)
			# if gini<-0.14*m20+0.33 and gini>0.14*m20+0.8:#pair_ctg[pair_ctg['id']==idx]['count_spec'][0]>0 and 
			# 	count_specpair.append(idx)
	plt.xlim(0,-3)
	plt.ylim(0.3,0.8)
	# plt.savefig(img_path+field+'_gini_m20_regulars_13petro_3detectsigma_03eta.png')
	plt.savefig(img_path+field+'_gini_m20_sfrcoded.png')
# print(len(count_specpair)/len(regular_list['goodsn']))
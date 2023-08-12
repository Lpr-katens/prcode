# import os

# path = '/Users/lpr/Data/galfit/LIRGproject/goodss/'#SExtractor/goodsn_sex/'
# temp = os.listdir(path)
# gala = []
# for file in temp:
# 	if file.endswith('.DS_Store') or file.endswith('.fits'):
# 		print(file)
# 	else:
# 		gala.append(file)
# remove_list=['residual.fits','galfit.01','fit.log','galfit.02','galfit.03','subcomps.fits','galfit2.feedme','galfit3.feedme','galfit3.CONSTRAINTS']
# #remove_list=['subcomps.fits','galfit2.feedme']
# for num in range(0,len(gala)):
# 	temp = os.listdir(path+'/'+gala[num])
# 	for file2 in temp:
# 		if file2 in remove_list:# if file2 == 'fit.log' or file2 == 'galfit.01' or file2 == 'residual.fits' or file2 == 'galfit.02':
# 			os.remove(path+'/'+gala[num]+'/'+file2)
# 			print(path+'/'+gala[num]+'/'+file2+' has been removed')
# import os
# from shutil import copy
# fields_list = ['goodsn','goodss','egs']
# img_path = '/Users/lpr/Data/lirg_project/output/'
# for field in fields_list:
# 	files = os.listdir(img_path+field+'/'+field+'_f160w')
# 	os.mkdir('/Users/lpr/Data/SExtractor/'+field+'_f160w')
# 	for file in files:
# 		if file.endswith('.fits') and not file.endswith('_wht.fits') and not file.endswith('_photutils.fits'):
# 			copy(img_path+field+'/'+field+'_f160w/'+file,'/Users/lpr/Data/SExtractor/'+field+'_f160w/'+file)
# 			print(file+' done')

# import os
# path = '/Users/lpr/Data/SAMCS/GALFIT/f1800w/'#SExtractor/goodsn_sex/'
# temp = os.listdir(path)
# gala = []
# for file in temp:
# 	if file.endswith('.DS_Store') or file.endswith('.fits'):
# 		print(file)
# 	else:
# 		gala.append(file)
# remove_list=['residual.fits','galfit.01','fit.log','galfit.02','galfit.03','subcomps.fits','galfit2.feedme','galfit3.feedme','galfit3.CONSTRAINTS']
# #remove_list=['subcomps.fits','galfit2.feedme']
# for num in range(0,len(gala)):
# 	temp = os.listdir(path+'/'+gala[num])
# 	for file2 in temp:
# 		if file2 in remove_list:# if file2 == 'fit.log' or file2 == 'galfit.01' or file2 == 'residual.fits' or file2 == 'galfit.02':
# 			os.remove(path+'/'+gala[num]+'/'+file2)
# 			print(path+'/'+gala[num]+'/'+file2+' has been removed')

import os
bands = ['f115w','f150w','f200w','f277w','f356w','f444w']
remove_list=['residual.fits','galfit.01','fit.log','galfit.02','galfit.03','galfit.04']
path = '/Users/lpr/Data/lirg_project/output/egs/jwst_candels/galfit/'
for band in bands:
	temp = os.listdir(path+band)
	gala = []
	for file in temp:
		if file.endswith('.DS_Store') or file.endswith('.fits'):
			print(file)
		else:
			gala.append(file)	
	for num in range(0,len(gala)):
		temp = os.listdir(path+band+'/'+gala[num])
		for file2 in temp:
			if file2 in remove_list:
				os.remove(path+band+'/'+gala[num]+'/'+file2)
				print(path+band+'/'+gala[num]+'/'+file2+' has been removed')
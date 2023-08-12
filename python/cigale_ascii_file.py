# ---------------------------------------------------------
# 用光学到红外的所有测光数据，来对lirgz~1样本的所有AGN做SED拟合，看
# 是否AGN对远红外的贡献大到影响论文的结果。
# ---------------------------------------------------------
from astropy.io import fits, ascii
from astropy.table import Table as tb
import numpy as np

fields_list = ['goodsn','goodss','egs']#
candels_fields = {'goodsn':'gdn','goodss':'gds','egs':'egs'}
cata_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
cata_suffix = '_Huangall_candels_radec_van_modifyz.fits'
candels_path = '/Users/lpr/Data/lirg_project/intake/CANDELS/catalog/JFang_CANDELS_Data/'
# 每10个星系一组，防止程序崩溃
for field in fields_list:
	match_catalog = fits.getdata(cata_path+field+cata_suffix,1)
	id_ctg = fits.getdata(cata_path+field+'_id.fits',1)
	candels_catalog = fits.getdata(candels_path+candels_fields[field]+'_all.fits',1)
	match_catalog = match_catalog[np.isin(match_catalog['id'],id_ctg['id'])]
	herschel_ctg = fits.getdata('/Users/lpr/Data/lirg_project/intake/huang_catalog/'+field+'_f16_zeq1.cat')
	# 定义列
	idx_list = np.full(len(match_catalog),-999,dtype=int)
	redshift_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f435_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f435err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f606_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f606err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f775_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f775err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f814_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f814err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f105_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f105err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f125_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f125err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f140_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f140err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f160_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	f160err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	ks_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	kserr_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	ch1_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	ch1err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	ch2_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	ch2err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	ch3_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	ch3err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	ch4_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	ch4err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	cfht_u_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	cfht_uerr_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	cfht_g_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	cfht_gerr_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	cfht_r_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	cfht_rerr_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	cfht_i_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	cfht_ierr_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	cfht_z_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	cfht_zerr_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	irs_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	irserr_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	mips_24_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	mips_24err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# herschel_70_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# herschel_70err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# herschel_100_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# herschel_100err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# herschel_160_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# herschel_160err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# herschel_250_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# herschel_250err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# herschel_350_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# herschel_350err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# herschel_500_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# herschel_500err_list = np.full(len(match_catalog),-999.,dtype=np.float32)
	# 添加数据进数组
	column_list = candels_catalog.names
	for num1 in range(0,len(match_catalog)):
		idx = match_catalog[num1]['ID']
		idx_candels = match_catalog[num1]['ID_CANDELS']
		redshift = match_catalog[num1]['z_used']
		idx_list[num1] = idx
		redshift_list[num1] = redshift
		# gdn and gds field flux should be in capital FLUX
		if 'CFHT_u_flux' in column_list:
			cfht_u = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['CFHT_u_flux']
			cfht_uerr = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['CFHT_u_fluxerr']
			cfht_u_list[num1] = cfht_u
			cfht_uerr_list[num1] = cfht_uerr
		if 'CFHT_g_flux' in column_list:
			cfht_g = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['CFHT_g_flux']
			cfht_gerr = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['CFHT_g_fluxerr']
			cfht_g_list[num1] = cfht_g
			cfht_gerr_list[num1] = cfht_gerr
		if 'CFHT_i_flux' in column_list:
			cfht_i = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['CFHT_i_flux']
			cfht_ierr = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['CFHT_i_fluxerr']
			cfht_i_list[num1] = cfht_i
			cfht_ierr_list[num1] = cfht_ierr
		if 'CFHT_z_flux' in column_list:
			cfht_z = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['CFHT_z_flux']
			cfht_zerr = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['CFHT_z_fluxerr']
			cfht_z_list[num1] = cfht_z
			cfht_zerr_list[num1] = cfht_zerr
		if 'ACS_F435W_FLUX' in column_list:
			f435 = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['ACS_F435W_FLUX']
			f435err = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['ACS_F435W_FLUXERR']
			f435_list[num1] = f435
			f435err_list[num1] = f435err
		if 'ACS_F606W_FLUX' in column_list:
			f606 = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['ACS_F606W_FLUX']
			f606err = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['ACS_F606W_FLUXERR']
			f606_list[num1] = f606
			f606err_list[num1] = f606err
		if 'ACS_F775W_FLUX' in column_list:
			f775 = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['ACS_F775W_FLUX']
			f775err = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['ACS_F775W_FLUXERR']
			f775_list[num1] = f775
			f775err_list[num1] = f775err
		if 'ACS_F814W_FLUX' in column_list:
			f814 = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['ACS_F814W_FLUX']
			f814err = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['ACS_F814W_FLUXERR']
			f814_list[num1] = f814
			f814err_list[num1] = f814err
		if 'WFC3_F105W_FLUX' in column_list:
			f105 = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['WFC3_F105W_FLUX']
			f105err = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['WFC3_F105W_FLUXERR']
			f105_list[num1] = f105
			f105err_list[num1] = f105err
		if 'WFC3_F125W_FLUX' in column_list:
			f125 = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['WFC3_F125W_FLUX']
			f125err = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['WFC3_F125W_FLUXERR']
			f125_list[num1] = f125
			f125err_list[num1] = f125err
		if 'WFC3_F140W_FLUX' in column_list:
			f140 = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['WFC3_F140W_FLUX']
			f140err = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['WFC3_F140W_FLUXERR']
			f140_list[num1] = f140
			f140err_list[num1] = f140err
		if 'WFC3_F160W_FLUX' in column_list:
			f160 = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['WFC3_F160W_FLUX']
			f160err = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['WFC3_F160W_FLUXERR']
			f160_list[num1] = f160
			f160err_list[num1] = f160err
		if 'CFHT_Ks_FLUX' in column_list:
			ks = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['CFHT_Ks_FLUX']
			kserr = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['CFHT_Ks_FLUXERR']
			ks_list[num1] = ks
			kserr_list[num1] = kserr
		if 'IRAC_CH1_SCANDELS_FLUX' in column_list:
			ch1 = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['IRAC_CH1_SCANDELS_FLUX']
			ch1err = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['IRAC_CH1_SCANDELS_FLUXERR']
			ch1_list[num1] = ch1
			ch1err_list[num1] = ch1err
		if 'IRAC_CH2_SCANDELS_FLUX' in column_list:
			ch2 = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['IRAC_CH2_SCANDELS_FLUX']
			ch2err = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['IRAC_CH2_SCANDELS_FLUXERR']
			ch2_list[num1] = ch2
			ch2err_list[num1] = ch2err
		if 'IRAC_CH3_FLUX' in column_list:
			ch3 = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['IRAC_CH3_FLUX']
			ch3err = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['IRAC_CH3_FLUXERR']
			ch3_list[num1] = ch3
			ch3err_list[num1] = ch3err
		if 'IRAC_CH4_FLUX' in column_list:
			ch4 = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['IRAC_CH4_FLUX']
			ch4err = candels_catalog[np.where(candels_catalog['ID']==idx_candels)]['IRAC_CH4_FLUXERR']	
			ch4_list[num1] = ch4
			ch4err_list[num1] = ch4err
		# 黄老师表格中的IRS,MIPS,Herschel数据 
			irs_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['F16'][0]
			irserr_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['E16']
			mips_24_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['F24'][0]
			mips_24err_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['E24']
			# herschel_70_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['F70'][0]
			# herschel_70err_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['E70']
			# herschel_100_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['F100'][0]
			# herschel_100err_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['E100']
			# herschel_160_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['F160'][0]
			# herschel_160err_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['E160']
			# herschel_250_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['F250'][0]
			# herschel_250err_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['E250']
			# herschel_350_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['F350'][0]
			# herschel_350err_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['E350']
			# herschel_500_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['F500'][0]
			# herschel_500err_list[num1] = herschel_ctg[herschel_ctg['ID']==idx]['E500']
	# 生成表格
	data = tb()
	data['id'] = idx_list
	data['redshift'] = redshift_list
	if 'CFHT_u_flux' in column_list:
		data['cfht.megacam.u'] = cfht_u_list*10**-3
		data['cfht.megacam.u_err'] = cfht_uerr_list*10**-3
	if 'CFHT_g_flux' in column_list:
		data['cfht.megacam.g'] = cfht_g_list*10**-3
		data['cfht.megacam.g_err'] = cfht_gerr_list*10**-3
	if 'CFHT_r_flux' in column_list:
		data['cfht.megacam.r'] = cfht_r_list*10**-3
		data['cfht.megacam.r_err'] = cfht_rerr_list*10**-3
	if 'CFHT_z_flux' in column_list:
		data['cfht.megacam.z'] = cfht_z_list*10**-3
		data['cfht.megacam.z_err'] = cfht_zerr_list*10**-3
	if 'ACS_F435W_FLUX' in column_list:
		data['hst.wfc.F435W'] = f435_list*10**-3
		data['hst.wfc.F435W_err'] = f435err_list*10**-3
	if 'ACS_F606W_FLUX' in column_list:
		data['hst.wfc.F606W'] = f606_list*10**-3
		data['hst.wfc.F606W_err'] = f606err_list*10**-3
	if 'ACS_F775W_FLUX' in column_list:
		data['hst.wfc.F775W'] = f775_list*10**-3
		data['hst.wfc.F775W_err'] = f775err_list*10**-3
	if 'ACS_F814W_FLUX' in column_list:
		data['hst.wfc.F814W'] = f814_list*10**-3
		data['hst.wfc.F814W_err'] = f814err_list*10**-3
	if 'WFC3_F105W_FLUX' in column_list:
		data['HST-WFC3_IR.F105W'] = f105_list*10**-3
		data['HST-WFC3_IR.F105W_err'] = f105err_list*10**-3
	if 'WFC3_F125W_FLUX' in column_list:
		data['hst.wfc3.F125W'] = f125_list*10**-3
		data['hst.wfc3.F125W_err'] = f125err_list*10**-3
	if 'WFC3_F140W_FLUX' in column_list:
		data['hst.wfc3.F140W'] = f140_list*10**-3
		data['hst.wfc3.F140W_err'] = f140err_list*10**-3
	if 'WFC3_F160W_FLUX' in column_list:
		data['hst.wfc3.F160W'] = f160_list*10**-3
		data['hst.wfc3.F160W_err'] = f160err_list*10**-3
	if 'CFHT_Ks_FLUX' in column_list:
		data['cfht.wircam.Ks'] = ks_list*10**-3
		data['cfht.wircam.Ks_err'] = kserr_list*10**-3
	if 'IRAC_CH1_SCANDELS_FLUX' in column_list:
		data['spitzer.irac.ch1'] = ch1_list*10**-3
		data['spitzer.irac.ch1_err'] = ch1err_list*10**-3
	if 'IRAC_CH2_SCANDELS_FLUX' in column_list:
		data['spitzer.irac.ch2'] = ch2_list*10**-3
		data['spitzer.irac.ch2_err'] = ch2err_list*10**-3
	if 'IRAC_CH3_FLUX' in column_list:
		data['spitzer.irac.ch3'] = ch3_list*10**-3
		data['spitzer.irac.ch3_err'] = ch3err_list*10**-3
	if 'IRAC_CH4_FLUX' in column_list:
		data['spitzer.irac.ch4'] = ch4_list*10**-3
		data['spitzer.irac.ch4_err'] = ch4err_list*10**-3
	irs_list = np.where(irs_list!=-99.,irs_list*10**-3,irs_list)
	data['spitzer.irs.Blue'] = irs_list
	irserr_list = np.where(irserr_list!=-99.,irserr_list*10**-3,irserr_list)
	data['spitzer.irs.Blue_err'] = irserr_list
	mips_24_list = np.where(mips_24_list!=-99.,mips_24_list*10**-3,mips_24_list)
	data['spitzer.mips.24'] = mips_24_list
	mips_24err_list = np.where(mips_24err_list!=-99.,mips_24err_list*10**-3,mips_24err_list)
	data['spitzer.mips.24_err'] = mips_24err_list
	# herschel_70_list = np.where(herschel_70_list!=-99.,herschel_70_list*10**-3,herschel_70_list)
	# data['herschel.pacs.70'] = herschel_70_list
	# herschel_70err_list = np.where(herschel_70err_list!=-99.,herschel_70err_list*10**-3,herschel_70err_list)
	# data['herschel.pacs.70_err'] = herschel_70err_list
	# herschel_100_list = np.where(herschel_100_list!=-99000.0,herschel_100_list*10**-3,herschel_100_list)
	# data['herschel.pacs.100'] = herschel_100_list
	# herschel_100err_list = np.where(herschel_100err_list!=-99000.0,herschel_100err_list*10**-3,herschel_100err_list)
	# data['herschel.pacs.100_err'] = herschel_100err_list
	# herschel_160_list = np.where(herschel_160_list!=-99000.0,herschel_160_list*10**-3,herschel_160_list)
	# data['herschel.pacs.160'] = herschel_160_list
	# herschel_160err_list = np.where(herschel_160err_list!=-99000.0,herschel_160err_list*10**-3,herschel_160err_list)
	# data['herschel.pacs.160_err'] = herschel_160err_list
	# herschel_250_list = np.where(herschel_250_list!=-99000.0,herschel_250_list*10**-3,herschel_250_list)
	# data['herschel.spire.250'] = herschel_250_list
	# herschel_250err_list = np.where(herschel_250err_list!=-99000.0,herschel_250err_list*10**-3,herschel_250err_list)
	# data['herschel.spire.250_err'] = herschel_250err_list
	# herschel_350_list = np.where(herschel_350_list!=-99000.0,herschel_350_list*10**-3,herschel_350_list)
	# data['herschel.spire.350'] = herschel_350_list
	# herschel_350err_list = np.where(herschel_350err_list!=-99000.0,herschel_350err_list*10**-3,herschel_350err_list)
	# data['herschel.spire.350_err'] = herschel_350err_list
	# herschel_500_list = np.where(herschel_500_list!=-99000.0,herschel_500_list*10**-3,herschel_500_list)
	# data['herschel.spire.500'] = herschel_500_list
	# herschel_500err_list = np.where(herschel_500err_list!=-99000.0,herschel_500err_list*10**-3,herschel_500err_list)
	# data['herschel.spire.500_err'] = herschel_500err_list
	# headline='# id redshift hst.wfc.F435W hst.wfc.F435W_err hst.wfc.F606W hst.wfc.F606W_err hst.wfc.F775W hst.wfc.F775W_err hst.wfc.F814W hst.wfc.F814W_err HST-WFC3_IR.F105W HST-WFC3_IR.F105W_err hst.wfc3.F125W hst.wfc3.F125W_err hst.wfc3.F140W hst.wfc3.F140W_err hst.wfc3.F160W hst.wfc3.F160W_err cfht.wircam.Ks cfht.wircam.Ks_err spitzer.irac.ch1 spitzer.irac.ch1_err spitzer.irac.ch2 spitzer.irac.ch2_err spitzer.irac.ch3 spitzer.irac.ch3_err spitzer.irac.ch4 spitzer.irac.ch4_err spitzer.irs.Blue.dat spitzer.irs.Blue.dat_err spitzer.mips.24 spitzer.mips.24_err herschel.pacs.70 herschel.pacs.70_err herschel.pacs.100 herschel.pacs.100_err herschel.pacs.160 herschel.pacs.160_err herschel.spire.250 herschel.spire.250_err herschel.spire.350 herschel.spire.350_err herschel.spire.500 herschel.spire.500_err'
	data.write('/Users/lpr/Data/lirg_project/output/cigale/'+field+'/'+field+'_photometry.txt',format='ascii',overwrite=True)
	print('------------- '+field+' done -------------')
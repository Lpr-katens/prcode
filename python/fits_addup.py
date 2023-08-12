# =============================================
#combine two fits binary tables from columns
# =============================================
from astropy.io import fits
import numpy as np
# from astropy.table import Table, hstack
fiels = ['goodsn','goodss','egs']
hdu1 = fits.open('/Users/lpr/Data/lirg_project/intake/huang_catalog/goodsn_f16_sfr_agn_20201217.fits')[1]
hdu2 = fits.open('/Users/lpr/Data/lirg_project/intake/huang_catalog/goodsn_f16_zeq1_candels.fits')[1]
hdu3 = fits.open('/Users/lpr/Data/lirg_project/intake/huang_catalog/goodsn_f16_zeq1_radec.cat')[1]
hdu4 = fits.open('/Users/lpr/Data/lirg_project/intake/huang_catalog/goodsn_fir_20200902.fits')[1]
hdu5 = fits.open('/Users/lpr/Data/lirg_project/intake/huang_catalog/goodsn_newfit_20200808.fits')[1]
hdu6 = fits.open('/Users/lpr/Data/lirg_project/intake/huang_catalog/goodsn_redshift.fits')[1]
hdu7 = fits.open('/Users/lpr/Data/lirg_project/intake/huang_catalog/goodsn_zeq1_x_20201120.fits')[1]
hdu = fits.BinTableHDU.from_columns(hdu1.data.columns+hdu2.data.columns+hdu3.data.columns+hdu4.data.columns+hdu5.data.columns+hdu6.data.columns+hdu7.data.columns)
hdu.writeto('/Users/lpr/Data/lirg_project/intake/huang_catalog/goodsn_Huang_all.fits',overwrite=True)
# =============================================
# combine two fits binary tables from rows
# =============================================
# from astropy.io import fits
# t1 = fits.open('')
# t2 = fits.open('')
# t3 = fits.open('')
# nrows1 = t1[1].data.shape[0]
# nrows2 = t2[1].data.shape[0]
# nrows3 = t3[1].data.shape[0]
# nrows = nrows1 + nrows2 + nrows3
# hdu = fits.BinTableHDU.from_columns(t1[1].columns, nrows=nrows)
# for colname in t1[1].columns.names:
#     hdu.data[colname][nrows1:] = t2[1].data[colname]
# for colname in t1[1].columns.names:
# 	hdu.data[colname][nrows1+nrows2:] = t3[1].data[colnam]
# hdu.writeto(' /Users/lpr/Data/fits/expdata/HST/RGBstamp/field3_combine.fits',overwrite=True)
# from astropy.io import fits
# import numpy as np

# hdu = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/tem_re_sersic.fits')[1].data
# class_candidate = {}
# col_tmp = np.zeros(len(hdu))
# col_ki = np.zeros(len(hdu))
# for num in range(0,len(hdu)):
#     col_tmp[num] = hdu[num]['TMP_CLASS'][0]
#     col_ki[num] = hdu[num]['KI'][0]
# # for num in range(0,len(hdu)):
# #     class_temp = [hdu[num]['TMP_CLASS'][0]]
# #     ki_temp = [hdu[num]['KI'][0]]
# # #     class_candidate[hdu[num]['galaxy_id']] = hdu[num]['TMP_CLASS'][0]
# #     for num2 in range(1,10):
# #         if np.abs(hdu[num]['KI'][num2]-hdu[num]['KI'][0]) <= 0.5: # if delta(ki)<0.5, take both classes as possible class
# # #             class_candidate[hdu[num]['galaxy_id']] = hdu[num]['TMP_CLASS'][num2]
# #             class_temp.append(hdu[num]['TMP_CLASS'][num2])
# #             ki_temp.append(hdu[num]['KI'][num2])
# #         else:
# #             break
# #     class_candidate[hdu[num]['galaxy_id']] = class_temp,ki_temp
# #     col_tmp[num] = np.array(class_temp)
# #     col_ki[num] = np.array(ki_temp)
# # # col_tmp = np.array(col_tmp)
# # # col_ki = np.array(col_ki)
# c1 = fits.Column(name='possible_class', array=col_tmp,format='D')
# c2 = fits.Column(name='possible_class_ki', array=col_ki,format='D')
# orig_cols = hdu.columns
# new_cols = fits.ColDefs([c1,c2])
# hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
# hdu.writeto('/Users/lpr/Data/fits/expdata/HST/goodsn_all/tem_re_sersic2.fits')
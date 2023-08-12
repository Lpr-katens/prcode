#find the best match between gds_l4p5ex and gds_f160f125.fits
from astropy.io import fits
import numpy as np
hdu1 = fits.open('/Users/lpr/Data/fits/pridata/LIRGs_project/egs_f16_sfr_agn_20201217.fits')[1].data
hdu2 = fits.open('/Users/lpr/Data/fits/pridata/LIRGs_project/egs_f16_zeq1_candels.fits')[1].data
hdu3 = fits.open('/Users/lpr/Data/fits/pridata/LIRGs_project/egs_newfit_20200808.fits')[1].data
# hdu2 = fits.open('/Users/lpr/Data/fits/pridata/LIRGs_project/JFang_CANDELS_Data/egs_all.fits')[1].data
# print('===========ready to find matches===========')
# ra = np.zeros(len(hdu1))
# dec = np.zeros(len(hdu1))
# re_160 = np.zeros(len(hdu1))
# dre_160 = np.zeros(len(hdu1))
# n_160 = np.zeros(len(hdu1))
# dn_160 = np.zeros(len(hdu1))
# q_160 = np.zeros(len(hdu1))
# dq_160 = np.zeros(len(hdu1))
# re_125 = np.zeros(len(hdu1))
# dre_125 = np.zeros(len(hdu1))
# n_125 = np.zeros(len(hdu1))
# dn_125 = np.zeros(len(hdu1))
# q_125 = np.zeros(len(hdu1))
# dq_125 = np.zeros(len(hdu1))
# delta_distance = np.zeros(len(hdu1))
# for num1 in range(0,len(hdu1)):
#     distance = 0
#     index = 0
#     for num2 in range(0,len(hdu2)):
#         if num2 == 0:
#             distance = np.sqrt((hdu1[num1]['RA']-hdu2[num2]['RA_1'])**2)
#             index = num2
#         else:
#             if np.sqrt((hdu1[num1]['RA']-hdu2[num2]['RA_1'])**2) < distance:
#                 distance = np.sqrt((hdu1[num1]['RA']-hdu2[num2]['RA_1'])**2)
#                 index = num2
#     ra[num1] = hdu2[index]['RA_1']
#     dec[num1] = hdu2[index]['DEC_1']
#     re_160[num1] = hdu2[index]['re_f160w']
#     dre_160[num1] = hdu2[index]['dre_f160w']
#     n_160[num1] = hdu2[index]['n_f160w']
#     dn_160[num1] = hdu2[index]['dn_f160w']
#     q_160[num1] = hdu2[index]['q_f160w']
#     dq_160[num1] = hdu2[index]['dq_f160w']
#     re_125[num1] = hdu2[index]['re_f125w']
#     dre_125[num1] = hdu2[index]['dre_f125w']
#     n_125[num1] = hdu2[index]['n_f125w']
#     dn_125[num1] = hdu2[index]['dn_f125w']
#     q_125[num1] = hdu2[index]['q_f125w']
#     dq_125[num1] = hdu2[index]['dq_f125w']
#     delta_distance[num1] = distance
# print('===========match done===========')
# col1 = fits.Column(name='ra', array=ra, format='D')
# col2 = fits.Column(name='dec', array=dec, format='D')
# col3 = fits.Column(name='re_f160w', array=re_160, format='D')
# col4 = fits.Column(name='dre_f160w', array=dre_160, format='D')
# col5 = fits.Column(name='n_f160w', array=n_160, format='D')
# col6 = fits.Column(name='dn_f160w', array=dn_160, format='D')
# col7 = fits.Column(name='q_f160w', array=q_160, format='D')
# col8 = fits.Column(name='dq_f160w', array=dq_160, format='D')
# col9 = fits.Column(name='re_f125w', array=re_125, format='D')
# col10 = fits.Column(name='dre_f125w', array=dre_125, format='D')
# col11 = fits.Column(name='n_f125w', array=n_125, format='D')
# col12 = fits.Column(name='dn_f125w', array=dn_125, format='D')
# col13 = fits.Column(name='q_f125w', array=q_125, format='D')
# col14 = fits.Column(name='dq_f125w', array=dq_125, format='D')
# col15 = fits.Column(name='distance',array=delta_distance,format='D')
# table = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15])
# table.writeto('/Users/lpr/Data/fits/expdata/HST/egs_all/egsTEMPtable.fits')
# print('===========table is wrtieto egsTEMPtable.fits=========== ')
# hdu3 = fits.open('/Users/lpr/Data/fits/expdata/HST/egs_all/egsTEMPtable.fits')[1]
hdu = fits.BinTableHDU.from_columns(hdu1.columns+hdu2.columns+hdu3.columns)
hdu.writeto('/Users/lpr/Data/fits/pridata/LIRGs_project/egs_Huang_all.fits')
print('===========table is wrtieto gds_Huang_all.fits=========== ')
#画gdn的sersic index和effective radius——galaxy type的直方图
from scipy.integrate import quad
import matplotlib.pyplot as plt
def dc(x,lambda0,m0):
    hz=np.sqrt(lambda0+m0*(1+x)**3)
    return 1/hz
H0=70 # km/s/Mpc, nowadays commonly used
Lambda0=0.7
M0=0.3
def kpc_per_arcsec(z):
    angular_distance=((3*10**5)/H0)*quad(dc,0,z,args=(Lambda0,M0))[0]/(1+z)
    arc_scale = angular_distance*np.pi*1000/(180*3600)
    return(arc_scale)
from astropy.io import fits
import numpy as np
hdu1 = fits.open('/Users/lpr/Data/fits/expdata/HST/goodss_all/gds_Huang_van.fits')[1].data
hdu2 = fits.open('/Users/lpr/Data/fits/expdata/HST/goodsn_all/gdn_Huang_van.fits')[1].data
hdu3 = fits.open('/Users/lpr/Data/fits/expdata/HST/egs_all/egs_Huang_van.fits')[1].data
# Huang(2021) selected galaxies with 0.8<z<1.3
hdu1 = hdu1[np.where((hdu1['offset']<1e-4)&(hdu1['REDSHIFT']>0.8)&(hdu1['REDSHIFT']<1.3))]
hdu2 = hdu2[np.where((hdu2['offset']<1e-4)&(hdu2['REDSHIFT']>0.8)&(hdu2['REDSHIFT']<1.3))]
hdu3 = hdu3[np.where((hdu3['offset']<1e-4)&(hdu3['REDSHIFT']>0.8)&(hdu3['REDSHIFT']<1.3))]
hdu1 = hdu1[np.where((hdu1['f_f160w']==1)|(hdu1['f_f160w']==0))]#
hdu2 = hdu2[np.where((hdu2['f_f160w']==1)|(hdu2['f_f160w']==0))]#
hdu3 = hdu3[np.where((hdu3['f_f160w']==1)|(hdu3['f_f160w']==0))]#
for num1 in range(0,len(hdu1)):
    hdu1[num1]['re_f160w'] = hdu1[num1]['re_f160w']*kpc_per_arcsec(hdu1[num1]['REDSHIFT'])
for num2 in range(0,len(hdu2)):
    hdu2[num2]['re_f160w'] = hdu2[num2]['re_f160w']*kpc_per_arcsec(hdu2[num2]['REDSHIFT'])
for num3 in range(0,len(hdu3)):
    hdu3[num3]['re_f160w'] = hdu3[num3]['re_f160w']*kpc_per_arcsec(hdu3[num3]['REDSHIFT'])
agn1 = hdu1[np.where(hdu1['TMP_CLASS'][:,0] == 1)]
agn2 = hdu2[np.where(hdu2['TMP_CLASS'][:,0] == 1)]
agn3 = hdu3[np.where(hdu3['TMP_CLASS'][:,0] == 1)]
#agn2 = agn1[np.where(agn1['L4P5EX']/agn1['E4P5EX']>3)]
composite1 = hdu1[np.where(hdu1['TMP_CLASS'][:,0] == 2)]
composite2 = hdu2[np.where(hdu2['TMP_CLASS'][:,0] == 2)]
composite3 = hdu3[np.where(hdu3['TMP_CLASS'][:,0] == 2)]
sf1 = hdu1[np.where(hdu1['TMP_CLASS'][:,0] == 3)]
sf2 = hdu2[np.where(hdu2['TMP_CLASS'][:,0] == 3)]
sf3 = hdu3[np.where(hdu3['TMP_CLASS'][:,0] == 3)]
quiescent1 = hdu1[np.where(hdu1['TMP_CLASS'][:,0] == 4)]
quiescent2 = hdu2[np.where(hdu2['TMP_CLASS'][:,0] == 4)]
quiescent3 = hdu3[np.where(hdu3['TMP_CLASS'][:,0] == 4)]
blue_compact1 = hdu1[np.where(hdu1['TMP_CLASS'][:,0] == 5)]
blue_compact2 = hdu2[np.where(hdu2['TMP_CLASS'][:,0] == 5)]
blue_compact3 = hdu3[np.where(hdu3['TMP_CLASS'][:,0] == 5)]
agn_n = np.concatenate((agn1['n_f160w'],agn2['n_f160w'],agn3['n_f160w']),axis=0)
composite_n = np.concatenate((composite1['n_f160w'],composite2['n_f160w'],composite3['n_f160w']),axis=0)
sf_n = np.concatenate((sf1['n_f160w'],sf2['n_f160w'],sf3['n_f160w']),axis=0)
blue_compact_n = np.concatenate((blue_compact1['n_f160w'],blue_compact2['n_f160w'],blue_compact3['n_f160w']),axis=0)
# agn_n = np.concatenate((agn1['re_f160w'],agn2['re_f160w'],agn3['re_f160w']),axis=0)
# composite_n = np.concatenate((composite1['re_f160w'],composite2['re_f160w'],composite3['re_f160w']),axis=0)
# sf_n = np.concatenate((sf1['re_f160w'],sf2['re_f160w'],sf3['re_f160w']),axis=0)
# blue_compact_n = np.concatenate((blue_compact1['re_f160w'],blue_compact2['re_f160w'],blue_compact3['re_f160w']),axis=0)
# #-------------------------------------------------------sersic index---------------------------------------------------
# plt.figure(figsize=(13,13))
# ax1 = plt.subplot(311)#=================
# plt.hist(agn_n,label='AGN',color='red',bins=15,range=[0,8],align='right',width=0.4)#re:range to 15
# plt.legend(fontsize=20)
# plt.text(7,16,'ALL: '+str(len(hdu1)+len(hdu2)+len(hdu3)),fontsize=20)
# plt.setp(ax1.get_xticklabels(), visible=False)
# plt.title('sersic index fraction distribution(%)',fontsize=20)
# plt.xlim(0,8.5)
# plt.yticks(ticks=np.arange(0,16,5),labels=np.around(100*np.arange(0,16,5)/(len(hdu)+len(hdu2)+len(hdu3)),1),fontsize=20)
# ax2 = plt.subplot(312,sharex=ax1)#=================
# plt.hist(composite_n,label='composite',color='orange',bins=15,range=[0,8],align='right',width=0.4)
# plt.legend(fontsize=20)
# plt.setp(ax2.get_xticklabels(), visible=False)
# plt.yticks(ticks=np.arange(0,65,32),labels=np.around(100*np.arange(0,65,32)/(len(hdu)+len(hdu2)+len(hdu3)),1),fontsize=20)
# ax3 = plt.subplot(313,sharex=ax1)#=================
# plt.hist(sf_n,label='star forming',color='green',bins=15,range=[0,8],align='right',width=0.4)
# plt.legend(fontsize=20)
# plt.yticks(ticks=np.arange(0,71,35),labels=np.around(100*np.arange(0,71,35)/(len(hdu)+len(hdu2)+len(hdu3)),1),fontsize=20)
# plt.xlabel('n',fontsize=20)
# plt.xticks(ticks=np.arange(0,8.5,1),fontsize=20)
# #plt.savefig('/Users/lpr/Desktop/yixing-meeting/n_freq-5class_exc.eps')
# plt.savefig('/Users/lpr/Data/fits/expdata/HST/egs_all/plot/n_freq-5class_exc.eps')
# #-------------------------------------------------------re---------------------------------------------------
# plt.figure(figsize=(13,13))
# ax1 = plt.subplot(311)#====================
# plt.hist(agn_n,label='AGN',color='red',bins=15,range=[0,15],align='right',width=0.7)
# plt.legend(fontsize=20)
# plt.plot([4,4],[0,23],linewidth=3)
# plt.text(4.5,12,'72.72%',fontsize=20)
# plt.text(4.5,15,'$r_e$<4',fontsize=20)
# # plt.text(12,22,'ALL: '+str(len(hdu1)+len(hdu2)+len(hdu3)),fontsize=20)
# plt.setp(ax1.get_xticklabels(), visible=False)
# # plt.title('$r_e$ fraction distribution(%)',fontsize=20)
# plt.xlim(0,15)
# plt.yticks(ticks=np.arange(0,23,7),labels=np.around(100*np.arange(0,23,7)/(len(hdu)+len(hdu2)+len(hdu3)),1),fontsize=20)
# ax2 = plt.subplot(312,sharex=ax1)#====================
# plt.hist(composite_n,label='composite',color='orange',bins=15,range=[0,15],align='right',width=0.7)
# plt.legend(fontsize=20)
# plt.plot([4,4],[0,47],linewidth=3)
# plt.text(4.5,25,'63.58%',fontsize=20)
# plt.text(4.5,31,'$r_e$<4',fontsize=20)
# plt.setp(ax2.get_xticklabels(), visible=False)
# plt.yticks(ticks=np.arange(0,47,20),labels=np.around(100*np.arange(0,47,20)/(len(hdu)+len(hdu2)+len(hdu3)),1),fontsize=20)
# ax3 = plt.subplot(313,sharex=ax1)#====================
# plt.hist(sf_n,label='star forming',color='green',bins=15,range=[0,15],align='right',width=0.7)
# plt.legend(fontsize=20)
# plt.plot([4,4],[0,62],linewidth=3)
# plt.text(4.5,40,'59.19%',fontsize=20)
# plt.text(4.5,48,'$r_e$<4',fontsize=20)
# plt.yticks(ticks=np.arange(0,62,30),labels=np.around(100*np.arange(0,62,30)/(len(hdu)+len(hdu2)+len(hdu3)),1),fontsize=20)
# plt.xlabel('$r_e$ (kpc)',fontsize=30)
# plt.xticks(ticks=np.arange(0,16,3),fontsize=20)
# #plt.savefig('/Users/lpr/Desktop/yixing-meeting/re_freq-5class_exc_text.eps')
# #plt.savefig('/Users/lpr/Data/fits/expdata/HST/egs_all/plot/re_freq-5class_exc.eps')
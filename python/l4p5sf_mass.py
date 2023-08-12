import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy.stats import linregress
hdu_gdn=fits.open('/Users/lpr/Data/lirg_project/output/catalog/goodsn_Huangall_candels_radec_van_zmasssfr8sfrhuangsfr.fits')[1].data
hdu_gds=fits.open('/Users/lpr/Data/lirg_project/output/catalog/goodss_Huangall_candels_radec_van_zmasssfr8sfrhuangsfr.fits')[1].data
hdu_egs=fits.open('/Users/lpr/Data/lirg_project/output/catalog/egs_Huangall_candels_radec_van_zmasssfr8sfrhuangsfr.fits')[1].data
hdu_gdn=hdu_gdn[np.where((hdu_gdn['z_used']>=0.8)&(hdu_gdn['z_used']<=1.3)&(hdu_gdn['id']!=-1)&(hdu_gdn['separation_candels_16']<1))]
hdu_gds=hdu_gds[np.where((hdu_gds['z_used']>=0.8)&(hdu_gds['z_used']<=1.3)&(hdu_gds['id']!=-1)&(hdu_gds['separation_candels_16']<1))]
hdu_egs=hdu_egs[np.where((hdu_egs['z_used']>=0.8)&(hdu_egs['z_used']<=1.3)&(hdu_egs['id']!=-1)&(hdu_egs['separation_candels_16']<1))]
# ensure the star forming 4.5um luminosity vs. mass relation do not contaminated by AGN
nonagn_gdn=hdu_gdn[np.where((hdu_gdn['L4P5EX']/hdu_gdn['E4P5EX']<3)&(hdu_gdn['TMP_CLASS'][:,0]!=1))]
nonagn_gds=hdu_gds[np.where((hdu_gds['L4P5EX']/hdu_gds['E4P5EX']<3)&(hdu_gds['TMP_CLASS'][:,0]!=1))]
nonagn_egs=hdu_egs[np.where((hdu_egs['L4P5EX']/hdu_egs['E4P5EX']<3)&(hdu_egs['TMP_CLASS'][:,0]!=1))]
agn_gdn=hdu_gdn[np.where(hdu_gdn['TMP_CLASS'][:,0]==1)]
agn_gds=hdu_gds[np.where(hdu_gds['TMP_CLASS'][:,0]==1)]
agn_egs=hdu_egs[np.where(hdu_egs['TMP_CLASS'][:,0]==1)]
# calculate correlation
# initialize a linear fitter
fit = fitting.LinearLSQFitter()
# initialize the outlier removal fitter
or_fit = fitting.FittingWithOutlierRemoval(fit,sigma_clip,niter=3,sigma=3.0)
# initialize a linear model
line_init = models.Linear1D()
# fit the data with the fitter
mass=np.concatenate((nonagn_gdn['LMASS_candels'],nonagn_gds['LMASS_candels'],nonagn_egs['LMASS_candels']),axis=0)
log_l4p5sf=np.concatenate((np.log10((10.**(nonagn_gdn['L4P5'][:,0])-nonagn_gdn['L4P5EX'])),np.log10((10.**(nonagn_gds['L4P5'][:,0])-nonagn_gds['L4P5EX'])),np.log10((10.**(nonagn_egs['L4P5'][:,0])-nonagn_egs['L4P5EX']))),axis=0)
log_l4p5sferr=np.concatenate(((np.sqrt((10.**nonagn_gdn['L4P5'][:,0]*np.log(10)*nonagn_gdn['E4P5'][:,0])**2+(nonagn_gdn['E4P5EX'])**2))/((10.**nonagn_gdn['L4P5'][:,0]-nonagn_gdn['L4P5EX'])*np.log(10)),(np.sqrt((10.**nonagn_gds['L4P5'][:,0]*np.log(10)*nonagn_gds['E4P5'][:,0])**2+(nonagn_gds['E4P5EX'])**2))/((10.**nonagn_gds['L4P5'][:,0]-nonagn_gds['L4P5EX'])*np.log(10)),(np.sqrt((10.**nonagn_egs['L4P5'][:,0]*np.log(10)*nonagn_egs['E4P5'][:,0])**2+(nonagn_egs['E4P5EX'])**2))/((10.**nonagn_egs['L4P5'][:,0]-nonagn_egs['L4P5EX'])*np.log(10))),axis=0)
fitted_line,mask = or_fit(line_init,log_l4p5sf,mass,weights=1.0/log_l4p5sferr)
# fitted_line,mask = or_fit(line_init,mass,log_l4p5sf,weights=1.0/log_l4p5sferr)
filtered_data = np.ma.masked_array(log_l4p5sf,mask=mask)
# plot
plt.figure()
plt.errorbar(mass,log_l4p5sf,yerr=log_l4p5sferr, fmt="ko", fillstyle="none", label="Clipped Data")
plt.plot(mass, filtered_data, "ko", label="Fitted Data")
plt.plot(mass, fitted_line(mass), 'k-', label='Fitted Model')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()
# use scipy.linregress to do fitting without errors
# result=linregress(mass,log_l4p5sf)
# x=np.linspace(np.min(mass),np.max(mass),50)
# y=result.intercept+result.slope*x
# plt.figure()
# plt.errorbar(mass,log_l4p5sf,yerr=log_l4p5sferr,fmt="ko",fillstyle="none",label="Data")
# plt.plot(x,y,"-",label="Fitting result")
# plt.text(10.2,9,'y='+str(round(result.intercept,2))+'+'+str(round(result.slope,2))+'x',color='red')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.legend()
# plt.show()
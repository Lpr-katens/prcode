# 根据EGS的图像的3sigma的面亮度来做星系的segmentation
# 大体思路是这样的：
# 1.先根据EGS F160W的图像，得到背景的高斯分布的均值。
# 2.然后对于GOODS场和EGS场所有的目标源，先做SExtractor，得到附近的污染源的segmentation然后mask掉。
# 	这个时候，即使对于污染源的扣除不是特别干净，后边应用3sigma的时候应该也没有太大的影响（待定）。
# 3.把第一步得到的背景的值拿来，应用在GOODS场和EGS场的所有目标源的cutout。
# 	只要是大于3sigma的像素都保留，然后先保存在一个文件夹里检查一下。
# 4.上一步得到的segmentation，用photutils.isophote.ellipse来拟合星系的椭圆。
# 	得到星系的半长轴作为星系的size
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from photutils.isophote import Ellipse
from photutils.morphology import data_properties
from photutils.aperture import EllipticalAperture as ea
import numpy as np
from lpr.image.display import logstretch
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from photutils.isophote import EllipseGeometry #用来展示星系的geometry

image_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
export_path = '/Users/lpr/Data/lirg_project/output/galaxy_size/'
segmap_path = '/Users/lpr/Data/lirg_project/output/van_re_half_light_radius/segmentation/'
re_path = '/Users/lpr/Data/lirg_project/output/van_re_half_light_radius/'
fields = ['goodss']#'goodsn',,'egs'
band_used = 'f160w'

# 用来拟合背景
def gaussian(x,mu,sigma,param):
    return param*np.exp(-((x-mu)/sigma)**2)/(sigma*np.sqrt(2*np.pi))
# 先分析EGS F160W图像的背景
egs_img = fits.getdata('/Users/lpr/Data/lirg_project/intake/CANDELS/egs/egs_all_wfc3_ir_f160w_060mas_v1.1_drz.fits',0)
noise = egs_img[(egs_img>-0.01)&(egs_img<0.01)&(egs_img!=0)]
ydata,t = np.histogram(noise,bins=100,range=[np.min(noise),np.max(noise)])
xdata = t[1:]
popt,pcov = curve_fit(gaussian,xdata,ydata)
bkg,bkg_sigma = popt[0],popt[1] #bkg,bkg_sigma 就是之后要用的背景的值,bkg+3*bkg_sigma

for field in fields:
	root_path = image_path+field+'/'+field+'_'+band_used
	id_ctg = fits.getdata(image_path+field+'_id.fits',1)
	h_ctg = fits.getdata(image_path+field+'_Huangall_candels_radec_van.fits',1)
	re_ctg = fits.getdata(re_path+field+'_re_hlr.fits',1)
	idx_list = np.full(len(id_ctg),-99)
	size_list = np.full(len(id_ctg),-99.)
	for num in range(0,len(id_ctg)):
		idx = id_ctg[num]['id']
		if idx!=32411:#!=32411
			print('ID:'+str(idx))
			print('         *              *\n'
'        %$@            !@#\n'
'       &%*&%          ^$#!@\n'
'      <:{!+@"|&^*&%$#!@%"}#)\n'
'     #!@%$**************&%$@~\n'
'    &++&@%***  code  ***%&%*&%\n'
'    @%*#!&**   begin  **"?)^&%\n'
'     %*#!&**   now    **?%*@#\n'
'      ++&@**************@~#^\n'
'          #!*&&^$#!@%*()\n'
'            $%(_+(%@?^\n'
'         )!(#*&++&@%%!*_#+\n'
'      !@#$)(*&^$#!@%*()*&%$@\n'
'    (*&^$#!@%*()*&%$@?#^%*(*&*\n'
'  #$)(*&^$#<{}"?)^&%$@?#^%*(*&*)\n'
' !@#$)(*&^$#!@%*&%*&%$@~#^%*(*&*)\n'
'  )(*&^$#!@%*$#*&^$#!@%*()*&%$@~\n'
'   $)(*&^$#!@%*!@*&%$@|#^%&%$@~\n'
'    (*&^$#!@%*@|_>"}#)$(*()*&%          **   **\n'
'     (*&^$#!@%*#!*&&^$#!@%*()       **          **\n'
'      )*&%$@|%*@%*&&^$#!@%*(      **      **     **\n'
'       @|#^%&%$>:}?%*@#!@%*     **       **      **\n'
'        @%*@#!@%*!@%*$#*&@*     **         **   **\n'
'        !@%*()*&:>%#!~>:"& ** **\n'
'        <:{!+@"|&^*&%$#!@%   *\n'
'       *$#!+@"*&  &%$@~<{}"')
			ra,dec = h_ctg[h_ctg['id']==idx]['ra_candels'],h_ctg[h_ctg['id']==idx]['dec_candels']
			img = fits.getdata(root_path+'/'+field+'_'+band_used+'_'+str(idx)+'.fits',0)
			hdr = fits.getheader(root_path+'/'+field+'_'+band_used+'_'+str(idx)+'.fits',0)
			hdr['CTYPE1'] = 'RA---TAN-SIP'
			hdr['CTYPE2'] = 'DEC--TAN-SIP'
			x,y = WCS(hdr).wcs_world2pix(ra,dec,0)
			segmap = fits.getdata(segmap_path+str(idx)+'.fits',0)
			# 对于是污染源的位置，赋值为-99，不参与后续的计算
			img[(segmap!=segmap[int(x),int(y)])&(segmap!=0)]=-99
			# 对于低于背景3sigma的像素，也赋值为-99，不参与后续的计算
			img[img<=bkg+3*bkg_sigma]=-99
			# 因为会有一些离散的像素，它们也高于背景3sigma，所以用一下statmorph中用的平滑segmentation的方法
			# 也就是scipy.ndimage.uniform_filter
			image = np.copy(img)
			image = np.where(image==-99,0,1) #把不属于目标星系的像素赋值0，其他地方赋值为1
			image_float = ndi.uniform_filter(np.float64(image),size=5)
			image = np.int64(image_float>0.5)
			img[image==0]=-99
			# # 如果旁边污染源对目标星系影响太大，那么画一个框，只拟合框里的
			# temp_img = np.copy(img)
			# img = temp_img[10:91,10:91]
			# 接下来正常做椭圆的拟合
			plt.figure(figsize=[5,5])
			plt.imshow(logstretch(img,a=500),cmap='binary')
			plt.savefig(export_path+field+'_'+str(idx)+'.png')
			plt.close()
			hdu = fits.HDUList([fits.PrimaryHDU(img)])
			hdu.writeto(export_path+field+'_'+str(idx)+'.fits',overwrite=True)
			cat = data_properties(img)
			position = (cat.xcentroid, cat.ycentroid)
			r = 3.0  # approximate isophotal extent
			a = cat.semimajor_sigma.value*r
			b = cat.semiminor_sigma.value*r
			theta = cat.orientation.to(u.rad).value
			apertures = ea(position, a, b, theta=theta)
			plt.figure(figsize=[5,5])
			plt.imshow(img, origin='lower', cmap='viridis',interpolation='nearest')
			# apertures.plot(color='#d62728')
			# # 画geometry
			# geometry = EllipseGeometry(x0=cat.xcentroid,y0=cat.ycentroid,sma=17.364374154231625,eps=(1-b/a),pa=cat.orientation.to(u.rad).value)
			# aper = ea((geometry.x0,geometry.y0),geometry.sma,geometry.sma*(1-geometry.eps),geometry.pa+10)
			# aper.plot(color='#d62728')
			plt.savefig(export_path+field+'_'+str(idx)+'_elli.png')
			plt.close()
			idx_list[num] = idx
			size_list[num] = a
			print(idx,a)
	col1 = fits.Column(name='id',array=idx_list,format='K')
	col2 = fits.Column(name='sma',array=size_list,format='D')
	hdu = fits.BinTableHDU.from_columns([col1,col2])
	hdu.writeto(export_path+field+'_sma.fits',overwrite=True)
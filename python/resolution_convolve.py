# ================
# improved version of image_convolution
# ================
# f606 resolution: 0.124"(2.1pixel), f160 resolution: 0.184"(3.1pixel), goodsn pixel scale: 0.059"/pix
import os
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.convolution import convolve_fft
import numpy as np
from photutils.psf import create_matching_kernel,HanningWindow
# smooth image with a gauussian kernel to make different wavelength resolution consistent
# convolution need to operate in same thing, such as pixel or arcsec
# first time, i try to convolve pixel(image) with arcsec, these two
# parameters are not same thing, so the result is wrong
# In gaussian PSF, FWHM=2*sqrt(2*np.log(2))*sigma, from Aniano,G.(2011)
# x_stddev: x standard deviation
# psf = {'f275w':1.24,'f435w':1.68 ,'f606w':1.9,'f775w':1.9,'f814w':1.9,'f850l':1.9,'f098m':2.53,'f105w':2.88,'f110w':2.96,'f125w':3.21,'f127m':3.21,'f139m':3.21,'f140w':3.21,'f153m':3.21,'f160w':3.22}
window = HanningWindow()
fields_list = ['goodsn','goodss','egs']
original_path = '/Users/lpr/Data/lirg_project/output/'
expath = '/Users/lpr/Data/lirg_project/output/'
psf_list = '/Users/lpr/Data/lirg_project/intake/CANDELS/'
for field in fields_list:
    if field == 'goodsn':# if field != 'egs':
        psf_160 = fits.open(psf_list+field+'_3dhst_v4.0_wfc3_psf/'+field+'_3dhst.v4.0.F160W_psf.fits')[0].data
    #     print(psf_160.shape)
    # else:
    #     psf_160 = fits.open(psf_list+'aegis_3dhst_v4.0_wfc3_psf/aegis_3dhst.v4.0.F160W_psf.fits')[0].data
    #     print(psf_160.shape)
    for bands_dir in os.listdir(original_path+field):
        band = bands_dir[bands_dir.find('_')+1:]
        print('band: '+band)
        if band=='f814w':
            if field == 'goodsn':
                psf_name = field+'_3dhst.v4.0.'+band.upper()+'_psf.fits'
                print('psf_name: '+psf_name)
                if psf_name in os.listdir(psf_list+field+'_3dhst_v4.0_acs_psf/'):
                    psf_band = fits.open(psf_list+field+'_3dhst_v4.0_acs_psf/'+psf_name)[0].data
                    print(psf_band.shape)
                    psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
                    img_list = os.listdir(original_path+field+'/'+bands_dir)
                    for num1 in range(0,len(img_list)):
                        if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
                            img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
                            astropy_conv = convolve(img.data,psf_kernel)
                            print(astropy_conv.shape)
                            img.data = astropy_conv
                            img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
                            print(img_list[num1]+'_is done')
                elif psf_name[:-9]+'cand_psf.fits' in os.listdir(psf_list+field+'_3dhst_v4.0_acs_psf/'):
                    psf_band = fits.open(psf_list+field+'_3dhst_v4.0_acs_psf/'+psf_name[:-9]+'cand_psf.fits')[0].data
                    print(psf_band.shape)
                    psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
                    img_list = os.listdir(original_path+field+'/'+bands_dir)
                    for num1 in range(0,len(img_list)):
                        if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
                            img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
                            astropy_conv = convolve(img.data,psf_kernel)
                            img.data = astropy_conv
                            img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
                            print(img_list[num1]+'_is done')
        # if band!='f125w' and band!='f140w' and band!='f160w' and band!='f850l':
        #     if field != 'egs':
        #         psf_name = field+'_3dhst.v4.0.'+band.upper()+'_psf.fits'
        #         print('psf_name: '+psf_name)
        #         if psf_name in os.listdir(psf_list+field+'_3dhst_v4.0_acs_psf/'):
        #             psf_band = fits.open(psf_list+field+'_3dhst_v4.0_acs_psf/'+psf_name)[0].data
        #             print(psf_band.shape)
        #             psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
        #             img_list = os.listdir(original_path+field+'/'+bands_dir)
        #             for num1 in range(0,len(img_list)):
        #                 if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
        #                     img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
        #                     astropy_conv = convolve(img.data,psf_kernel)
        #                     print(astropy_conv.shape)
        #                     img.data = astropy_conv
        #                     img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
        #                     print(img_list[num1]+'_is done')
        #         elif psf_name[:-9]+'cand_psf.fits' in os.listdir(psf_list+field+'_3dhst_v4.0_acs_psf/'):
        #             psf_band = fits.open(psf_list+field+'_3dhst_v4.0_acs_psf/'+psf_name[:-9]+'cand_psf.fits')[0].data
        #             print(psf_band.shape)
        #             psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
        #             img_list = os.listdir(original_path+field+'/'+bands_dir)
        #             for num1 in range(0,len(img_list)):
        #                 if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
        #                     img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
        #                     astropy_conv = convolve(img.data,psf_kernel)
        #                     img.data = astropy_conv
        #                     img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
        #                     print(img_list[num1]+'_is done')
        #     else:
        #         psf_name = 'aegis_3dhst.v4.0.'+band.upper()+'_psf.fits'
        #         if psf_name in os.listdir(psf_list+'aegis_3dhst_v4.0_acs_psf/'):
        #             psf_band = fits.open(psf_list+'aegis_3dhst_v4.0_acs_psf/'+psf_name)[0].data
        #             psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
        #             img_list = os.listdir(original_path+field+'/'+bands_dir)
        #             for num1 in range(0,len(img_list)):
        #                 if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
        #                     img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
        #                     astropy_conv = convolve(img.data,psf_kernel)
        #                     img.data = astropy_conv
        #                     img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
        #                     print(img_list[num1]+'_is done')
        #         elif psf_name[:-9]+'cand_psf.fits' in os.listdir(psf_list+'aegis_3dhst_v4.0_acs_psf/'):
        #             psf_band = fits.open(psf_list+'aegis_3dhst_v4.0_acs_psf/'+psf_name[:-9]+'cand_psf.fits')[0].data
        #             psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
        #             img_list = os.listdir(original_path+field+'/'+bands_dir)
        #             for num1 in range(0,len(img_list)):
        #                 if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
        #                     img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
        #                     astropy_conv = convolve(img.data,psf_kernel)
        #                     img.data = astropy_conv
        #                     img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
        #                     print(img_list[num1]+'_is done')
        # elif band=='f850l':
        #     if field != 'egs':
        #         psf_name = field+'_3dhst.v4.0.F850LP_psf.fits'
        #         print('psf_name: '+psf_name)
        #         if psf_name in os.listdir(psf_list+field+'_3dhst_v4.0_acs_psf/'):
        #             psf_band = fits.open(psf_list+field+'_3dhst_v4.0_acs_psf/'+psf_name)[0].data
        #             psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
        #             img_list = os.listdir(original_path+field+'/'+bands_dir)
        #             for num1 in range(0,len(img_list)):
        #                 if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
        #                     img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
        #                     astropy_conv = convolve(img.data,psf_kernel)
        #                     print(astropy_conv.shape)
        #                     img.data = astropy_conv
        #                     img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
        #                     print(img_list[num1]+'_is done')
        #         elif psf_name[:-9]+'cand_psf.fits' in os.listdir(psf_list+field+'_3dhst_v4.0_acs_psf/'):
        #             psf_band = fits.open(psf_list+field+'_3dhst_v4.0_acs_psf/'+psf_name[:-9]+'cand_psf.fits')[0].data
        #             print(psf_band.shape)
        #             psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
        #             img_list = os.listdir(original_path+field+'/'+bands_dir)
        #             for num1 in range(0,len(img_list)):
        #                 if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
        #                     img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
        #                     astropy_conv = convolve(img.data,psf_kernel)
        #                     img.data = astropy_conv
        #                     img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
        #                     print(img_list[num1]+'_is done')
        #     else:
        #         psf_name = 'aegis_3dhst.v4.0.F850LP_psf.fits'
        #         if psf_name in os.listdir(psf_list+'aegis_3dhst_v4.0_acs_psf/'):
        #             psf_band = fits.open(psf_list+'aegis_3dhst_v4.0_acs_psf/'+psf_name)[0].data
        #             psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
        #             img_list = os.listdir(original_path+field+'/'+bands_dir)
        #             for num1 in range(0,len(img_list)):
        #                 if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
        #                     img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
        #                     astropy_conv = convolve(img.data,psf_kernel)
        #                     img.data = astropy_conv
        #                     img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
        #                     print(img_list[num1]+'_is done')
        #         elif psf_name[:-9]+'cand_psf.fits' in os.listdir(psf_list+'aegis_3dhst_v4.0_acs_psf/'):
        #             psf_band = fits.open(psf_list+'aegis_3dhst_v4.0_acs_psf/'+psf_name[:-9]+'cand_psf.fits')[0].data
        #             psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
        #             img_list = os.listdir(original_path+field+'/'+bands_dir)
        #             for num1 in range(0,len(img_list)):
        #                 if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
        #                     img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
        #                     astropy_conv = convolve(img.data,psf_kernel)
        #                     img.data = astropy_conv
        #                     img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
        #                     print(img_list[num1]+'_is done')
        # elif band=='f125w' or band=='f140w' or band=='f160w':
        #     if field != 'egs':
        #         psf_name = field+'_3dhst.v4.0.'+band.upper()+'_psf.fits'
        #         print('psf_name: '+psf_name)
        #         if psf_name in os.listdir(psf_list+field+'_3dhst_v4.0_wfc3_psf/'):
        #             psf_band = fits.open(psf_list+field+'_3dhst_v4.0_wfc3_psf/'+psf_name)[0].data
        #             print(psf_band.shape)
        #             psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
        #             img_list = os.listdir(original_path+field+'/'+bands_dir)
        #             for num1 in range(0,len(img_list)):
        #                 if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
        #                     img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
        #                     astropy_conv = convolve(img.data,psf_kernel)
        #                     print(astropy_conv.shape)
        #                     img.data = astropy_conv
        #                     img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
        #                     print(img_list[num1]+'_is done')
        #         elif psf_name[:-9]+'cand_psf.fits' in os.listdir(psf_list+field+'_3dhst_v4.0_wfc3_psf/'):
        #             psf_band = fits.open(psf_list+field+'_3dhst_v4.0_wfc3_psf/'+psf_name[:-9]+'cand_psf.fits')[0].data
        #             psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
        #             img_list = os.listdir(original_path+field+'/'+bands_dir)
        #             for num1 in range(0,len(img_list)):
        #                 if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
        #                     img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
        #                     astropy_conv = convolve(img.data,psf_kernel)
        #                     img.data = astropy_conv
        #                     img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
        #                     print(img_list[num1]+'_is done')
        #     else:
        #         psf_name = 'aegis_3dhst.v4.0.'+band.upper()+'_psf.fits'
        #         if psf_name in os.listdir(psf_list+'aegis_3dhst_v4.0_wfc3_psf/'):
        #             psf_band = fits.open(psf_list+'aegis_3dhst_v4.0_wfc3_psf/'+psf_name)[0].data
        #             psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
        #             img_list = os.listdir(original_path+field+'/'+bands_dir)
        #             for num1 in range(0,len(img_list)):
        #                 if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
        #                     img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
        #                     astropy_conv = convolve(img.data,psf_kernel)
        #                     img.data = astropy_conv
        #                     img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
        #                     print(img_list[num1]+'_is done')
        #         elif psf_name[:-9]+'cand_psf.fits' in os.listdir(psf_list+'aegis_3dhst_v4.0_wfc3_psf/'):
        #             psf_band = fits.open(psf_list+'aegis_3dhst_v4.0_wfc3_psf/'+psf_name[:-9]+'cand_psf.fits')[0].data
        #             psf_kernel = create_matching_kernel(psf_band,psf_160,window=window)
        #             img_list = os.listdir(original_path+field+'/'+bands_dir)
        #             for num1 in range(0,len(img_list)):
        #                 if img_list[num1].endswith('.fits') and not img_list[num1].endswith('_wht.fits') and not img_list[num1].endswith('_photutils.fits'):
        #                     img = fits.open(original_path+field+'/'+bands_dir+'/'+img_list[num1])[0]
        #                     astropy_conv = convolve(img.data,psf_kernel)
        #                     img.data = astropy_conv
        #                     img.writeto(expath+field+'/'+bands_dir+'/'+img_list[num1][:-5]+'_photutils.fits',overwrite=True)
        #                     print(img_list[num1]+'_is done')


# # os.mkdir(expath)
# psf_105 = (fits.open('/Users/lpr/Data/lirg_project/output/CANDELS/PSFSTD_WFC3IR_F105W.fits')[0].data)[4]
# psf_606 = fits.open('/Users/lpr/Data/lirg_project/output/CANDELS/goodsn_3dhst.v4.0.F606W_psf.fits')[0].data
# psf_850 = fits.open('/Users/lpr/Data/lirg_project/output/CANDELS/goodsn_3dhst.v4.0.F850LP_psf.fits')[0].data
# psf_125 = fits.open('/Users/lpr/Data/lirg_project/output/CANDELS/goodsn_3dhst.v4.0.F125W_psf.fits')[0].data
# psf_140 = fits.open('/Users/lpr/Data/lirg_project/output/CANDELS/goodsn_3dhst.v4.0.F140W_psf.fits')[0].data
# psf_160 = fits.open('/Users/lpr/Data/lirg_project/output/CANDELS/goodsn_3dhst.v4.0.F160W_psf.fits')[0].data
# # psf_list = {'goodsn_f606w':psf_606,'goodsn_f850l':psf_850,'goodsn_f125w':psf_125,'goodsn_f140w':psf_140}
# # for band in band_list:
# band = 'f105w'
# files = os.listdir(original_path+'goodsn_'+band) # images
# img_list = []
# for items in files:
#     if items.endswith('.fits'):
#         img_list.append(items)
# psf_kernel = create_matching_kernel(psf_105[50-34:50+35,50-34:50+35],psf_160,window=window)# psf_kernel = create_matching_kernel(psf_list[band],psf_160,window=window)
# os.mkdir(expath+'goodsn_'+band+'_photutils/')
# for num in range(0,len(img_list)):
#     # band = filelist[num][filelist[num].index('_f')+1:]
    
#     # sigma = np.sqrt(psf['f125w']**2-psf['f606w']**2)/(2*np.sqrt(2*np.log(2))) # gaussian error sigma is different with FWHM
#     # print('get sigma of '+band+': '+str(sigma))
#     # kernel = Gaussian2DKernel(x_stddev=sigma,y_stddev=sigma,x_size=41,y_size=41)#int(10*sigma)
#     # print('get kernel of '+band)
#     hdu = fits.open(original_path+'goodsn_'+band+'/'+img_list[num])[0]
#     img = hdu.data
#     # astropy's convolution replace the NaN pixels with a kernel_weighted interpolation from their neighbors
#     # astropy_conv = convolve(img,kernel)#,allow_huge=True)
#     astropy_conv = convolve(img,psf_kernel)
#     # export the data which has been convolved
#     hdu.data = astropy_conv
#     # filename = 'convolution_'+ band + '.fits'
#     hdu.writeto(expath+'goodsn_'+band+'_photutils/'+img_list[num])
#     print(img_list[num]+' convolution is done')
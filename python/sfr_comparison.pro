hdu1=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/egs_Huangall_candels_radec_van_sfr_8sfr.fits',1)
hdu2=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/goodsn_Huangall_candels_radec_van_sfr_8sfr.fits',1)
hdu3=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/goodss_Huangall_candels_radec_van_sfr_8sfr.fits',1)
binsize=10
hdu1=hdu1[where(hdu1.z_used gt 0.8 and hdu1.z_used lt 1.3 and hdu1.separation_candels_16 lt 1 and hdu1.id ne -1)]
hdu2=hdu2[where(hdu2.z_used gt 0.8 and hdu2.z_used lt 1.3 and hdu2.separation_candels_16 lt 1 and hdu2.id ne -1)]
hdu3=hdu3[where(hdu3.z_used gt 0.8 and hdu3.z_used lt 1.3 and hdu3.separation_candels_16 lt 1 and hdu3.id ne -1)]
; use l8um star forming component to calculate sfr
; disk_sfr=[hdu1[where(hdu1.n_f160w gt 0 and hdu1.n_f160w lt 2)].sfr_8umsfr,hdu2[where(hdu2.n_f160w gt 0 and hdu2.n_f160w lt 2)].sfr_8umsfr,hdu3[where(hdu3.n_f160w gt 0 and hdu3.n_f160w lt 2)].sfr_8umsfr]
; bulge_sfr=[hdu1[where(hdu1.n_f160w gt 2 and hdu1.n_f160w lt 8)].sfr_8umsfr,hdu2[where(hdu2.n_f160w gt 2 and hdu2.n_f160w lt 8)].sfr_8umsfr,hdu3[where(hdu3.n_f160w gt 2 and hdu3.n_f160w lt 8)].sfr_8umsfr]
; use l8um to calculate sfr
; disk_sfr=[hdu1[where(hdu1.n_f160w gt 0 and hdu1.n_f160w lt 2)].sfr_8um,hdu2[where(hdu2.n_f160w gt 0 and hdu2.n_f160w lt 2)].sfr_8um,hdu3[where(hdu3.n_f160w gt 0 and hdu3.n_f160w lt 2)].sfr_8um]
; bulge_sfr=[hdu1[where(hdu1.n_f160w gt 2 and hdu1.n_f160w lt 8)].sfr_8um,hdu2[where(hdu2.n_f160w gt 2 and hdu2.n_f160w lt 8)].sfr_8um,hdu3[where(hdu3.n_f160w gt 2 and hdu3.n_f160w lt 8)].sfr_8um]
; use total IR luminosity to calculate sfr
disk_sfr=10.^[hdu1[where(hdu1.n_f160w gt 0 and hdu1.n_f160w lt 2)].ltir_full,hdu2[where(hdu2.n_f160w gt 0 and hdu2.n_f160w lt 2)].ltir_full,hdu3[where(hdu3.n_f160w gt 0 and hdu3.n_f160w lt 2)].ltir_full]*1.73*1e-10
bulge_sfr=10.^[hdu1[where(hdu1.n_f160w gt 2 and hdu1.n_f160w lt 8)].ltir_full,hdu2[where(hdu2.n_f160w gt 2 and hdu2.n_f160w lt 8)].ltir_full,hdu3[where(hdu3.n_f160w gt 2 and hdu3.n_f160w lt 8)].ltir_full]*1.73*1e-10
disk_sfr=disk_sfr[where(disk_sfr gt 0)]
bulge_sfr=bulge_sfr[where(bulge_sfr gt 0)]
disk_sfr=alog10(disk_sfr*0.63)
bulge_sfr=alog10(bulge_sfr*0.63)
; disk_sfr=disk_sfr[where(disk_sfr gt 0)];for TIR calculation
; bulge_sfr=bulge_sfr[where(bulge_sfr gt 0)]
max_sfr=max([disk_sfr,bulge_sfr])
min_sfr=min([disk_sfr,bulge_sfr])
; print,min_sfr,max_sfr
set_plot,'ps'
device,xsize=25,ysize=25,filename='/Users/lpr/Data/lirg_project/output/images/sfr_comparison.eps',/encapsulated
cghistoplot,disk_sfr,mininput=min_sfr,maxinput=max_sfr,binsize=(max_sfr-min_sfr)/binsize,color='blue',/frequency,xran=[0,2.7],yran=[0,0.3],yticks=4,ytickformat='(F4.2)',xtitle='!Nlog!S!R!I10'+'!NSFR!N('+'!NM!S!R!I'+sunsymbol()+'!N/!Nyr!N)',thick=3,charsize=2,charthick=3,title='!NSFR:L!S!!ESFR!I8um'
cghistoplot,bulge_sfr,mininput=min_sfr,maxinput=max_sfr,binsize=(max_sfr-min_sfr)/binsize,color='red',/frequency,/oplot,thick=3
cglegend,psyms=[0,0],length=0.05,thick=3,colors=['blue','red'],location=[1.7,0.25],titles=['disk galaxies','bulge galaxies'],/data
cgtext,0.01,0.25,'max SFR='+string(10.^max_sfr,format='(F6.2)'),charthick=3
cgtext,0.01,0.24,'min SFR='+string(10.^min_sfr,format='(F4.2)'),charthick=3
device,/close_file
; h_disk=histogram(disk_sfr,min=min_sfr,max=max_sfr,nbins=10);salpeter to chabrier correction
; h_bulge=histogram(bulge_sfr,min=min_sfr,max=max_sfr,nbins=10)
; cgplot,findgen(11)*(max_sfr-min_sfr)/10+min_sfr,h_disk/total(h_disk),color='blue',psym=10,yran=[0,0.5]
; cgoplot,findgen(11)*(max_sfr-min_sfr)/10+min_sfr,h_bulge/total(h_bulge),color='red',psym=10
end
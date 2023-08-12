hdu_gdn=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/goodsn_Huangall_candels_radec_van_zmasssfr8sfrhuangsfr.fits',1)
hdu_gds=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/goodss_Huangall_candels_radec_van_zmasssfr8sfrhuangsfr.fits',1)
hdu_egs=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/egs_Huangall_candels_radec_van_zmasssfr8sfrhuangsfr.fits',1)
hdu_gdn=hdu_gdn[where(hdu_gdn.z_used gt 0.8 and hdu_gdn.z_used lt 1.3 and hdu_gdn.id ne -1 and hdu_gdn.separation_candels_16 lt 1)]
hdu_gds=hdu_gds[where(hdu_gds.z_used gt 0.8 and hdu_gds.z_used lt 1.3 and hdu_gds.id ne -1 and hdu_gds.separation_candels_16 lt 1)]
hdu_egs=hdu_egs[where(hdu_egs.z_used gt 0.8 and hdu_egs.z_used lt 1.3 and hdu_egs.id ne -1 and hdu_egs.separation_candels_16 lt 1)]
; ensure the star forming 4.5um luminosity vs. mass relation do not contaminated by AGN
nonagn_gdn=hdu_gdn[where(hdu_gdn.l4p5ex/hdu_gdn.e4p5ex lt 3 and hdu_gdn.tmp_class[0] ne 1)]
nonagn_gds=hdu_gds[where(hdu_gds.l4p5ex/hdu_gds.e4p5ex lt 3 and hdu_gds.tmp_class[0] ne 1)]
nonagn_egs=hdu_egs[where(hdu_egs.l4p5ex/hdu_egs.e4p5ex lt 3 and hdu_egs.tmp_class[0] ne 1)]
agn_gdn=hdu_gdn[where(hdu_gdn.tmp_class[0] eq 1)]
agn_gds=hdu_gds[where(hdu_gds.tmp_class[0] eq 1)]
agn_egs=hdu_egs[where(hdu_egs.tmp_class[0] eq 1)]
; calculate correlation
pearson_coeff=correlate([alog10(10.^(nonagn_gdn.l4p5[0])-nonagn_gdn.l4p5ex),alog10(10.^(nonagn_gds.l4p5[0])-nonagn_gds.l4p5ex),alog10(10.^(nonagn_egs.l4p5[0])-nonagn_egs.l4p5ex)],[nonagn_gdn.lmass_candels,nonagn_gds.lmass_candels,nonagn_egs.lmass_candels])
linear_result=linfit([alog10((10.^(nonagn_gdn.l4p5[0])-nonagn_gdn.l4p5ex)),alog10((10.^(nonagn_gds.l4p5[0])-nonagn_gds.l4p5ex)),alog10((10.^(nonagn_egs.l4p5[0])-nonagn_egs.l4p5ex))],[nonagn_gdn.lmass_candels,nonagn_gds.lmass_candels,nonagn_egs.lmass_candels]);,sigma=sigma,MEASURE_ERRORS=[(sqrt((10.^nonagn_gdn.l4p5[0]*alog(10)*nonagn_gdn.e4p5[0])^2+(nonagn_gdn.e4p5ex)^2))/((10.^nonagn_gdn.l4p5[0]-nonagn_gdn.l4p5ex)*alog(10)),(sqrt((10.^nonagn_gds.l4p5[0]*alog(10)*nonagn_gds.e4p5[0])^2+(nonagn_gds.e4p5ex)^2))/((10.^nonagn_gds.l4p5[0]-nonagn_gds.l4p5ex)*alog(10)),(sqrt((10.^nonagn_egs.l4p5[0]*alog(10)*nonagn_egs.e4p5[0])^2+(nonagn_egs.e4p5ex)^2))/((10.^nonagn_egs.l4p5[0]-nonagn_egs.l4p5ex)*alog(10))])
; plot fitting results
set_plot, 'ps'
device,xsize=25,ysize=17,filename='/Users/lpr/Data/lirg_project/output/images/mass_correction_1.eps',/encapsulated
cgplot,[alog10(10.^(nonagn_gdn.l4p5[0])-nonagn_gdn.l4p5ex),alog10(10.^(nonagn_gds.l4p5[0])-nonagn_gds.l4p5ex),alog10(10.^(nonagn_egs.l4p5[0])-nonagn_egs.l4p5ex)],[nonagn_gdn.lmass_candels,nonagn_gds.lmass_candels,nonagn_egs.lmass_candels],color='red',title='mass MIR relation',xtitle='!Nlog!S!R!I10'+'!N(!NM!S!R!I*'+'!N/!NM'+sunsymbol()+')',ytitle='!Nlog!S!R!I10'+'!N(!NL!S!ESFR+stellar!R!I4.5'+'       '+'!N/!NL'+sunsymbol()+')',psym=4,xran=[8.5,11.2],yran=[9,11.3],charsize=2,thick=4;,err_yhigh=[(sqrt((10.^nonagn_gdn.l4p5[0]*alog(10)*nonagn_gdn.e4p5[0])^2+(nonagn_gdn.e4p5ex)^2))/((10.^nonagn_gdn.l4p5[0]-nonagn_gdn.l4p5ex)*alog(10)),(sqrt((10.^nonagn_gds.l4p5[0]*alog(10)*nonagn_gds.e4p5[0])^2+(nonagn_gds.e4p5ex)^2))/((10.^nonagn_gds.l4p5[0]-nonagn_gds.l4p5ex)*alog(10)),(sqrt((10.^nonagn_egs.l4p5[0]*alog(10)*nonagn_egs.e4p5[0])^2+(nonagn_egs.e4p5ex)^2))/((10.^nonagn_egs.l4p5[0]-nonagn_egs.l4p5ex)*alog(10))],err_ylow=[(sqrt((10.^nonagn_gdn.l4p5[0]*alog(10)*nonagn_gdn.e4p5[0])^2+(nonagn_gdn.e4p5ex)^2))/((10.^nonagn_gdn.l4p5[0]-nonagn_gdn.l4p5ex)*alog(10)),(sqrt((10.^nonagn_gds.l4p5[0]*alog(10)*nonagn_gds.e4p5[0])^2+(nonagn_gds.e4p5ex)^2))/((10.^nonagn_gds.l4p5[0]-nonagn_gds.l4p5ex)*alog(10)),(sqrt((10.^nonagn_egs.l4p5[0]*alog(10)*nonagn_egs.e4p5[0])^2+(nonagn_egs.e4p5ex)^2))/((10.^nonagn_egs.l4p5[0]-nonagn_egs.l4p5ex)*alog(10))]
cgoplot,[alog10(10.^(agn_gdn.l4p5[0])-agn_gdn.l4p5ex),alog10(10.^(agn_gds.l4p5[0])-agn_gds.l4p5ex),alog10(10.^(agn_egs.l4p5[0])-agn_egs.l4p5ex)],[agn_gdn.lmass_candels,agn_gds.lmass_candels,agn_egs.lmass_candels],color='purple',psym=4,thick=4;err_yhigh=[(sqrt((10.^agn_gdn.l4p5[0]*alog(10)*agn_gdn.e4p5[0])^2+(agn_gdn.e4p5ex)^2))/((10.^agn_gdn.l4p5[0]-agn_gdn.l4p5ex)*alog(10)),(sqrt((10.^agn_gds.l4p5[0]*alog(10)*agn_gds.e4p5[0])^2+(agn_gds.e4p5ex)^2))/((10.^agn_gds.l4p5[0]-agn_gds.l4p5ex)*alog(10)),(sqrt((10.^agn_egs.l4p5[0]*alog(10)*agn_egs.e4p5[0])^2+(agn_egs.e4p5ex)^2))/((10.^agn_egs.l4p5[0]-agn_egs.l4p5ex)*alog(10))],err_ylow=[(sqrt((10.^agn_gdn.l4p5[0]*alog(10)*agn_gdn.e4p5[0])^2+(agn_gdn.e4p5ex)^2))/((10.^agn_gdn.l4p5[0]-agn_gdn.l4p5ex)*alog(10)),(sqrt((10.^agn_gds.l4p5[0]*alog(10)*agn_gds.e4p5[0])^2+(agn_gds.e4p5ex)^2))/((10.^agn_gds.l4p5[0]-agn_gds.l4p5ex)*alog(10)),(sqrt((10.^agn_egs.l4p5[0]*alog(10)*agn_egs.e4p5[0])^2+(agn_egs.e4p5ex)^2))/((10.^agn_egs.l4p5[0]-agn_egs.l4p5ex)*alog(10))],
x=findgen(91)/20.+8
cgoplot,x,linear_result[0]+linear_result[1]*x,color='blue',thick=4
cgtext,8.6,11,'pearson correlation: '+string(pearson_coeff,format='(F4.2)'),color='blue'
cgtext,10,9.1,'linear fitting: y='+string(linear_result[0],format='(F4.2)')+'+'+string(linear_result[1],format='(F4.2)')+'x',color='blue'
cglegend,psyms=[4,4],length=0,thick=4,symthick=4,symsize=2,symcolors=['red','purple'],location=[10.7,9.4],titles=['nonAGNs','AGNs'],tcolors=['red','purple'],/data
cglegend,psyms=0,length=0.05,thick=4,symsize=2,colors='blue',location=[10.5,9.25],titles=['linear fitting'],tcolors='blue',/data
device,/close_file
; use fitting results to do mass correction

; agn_massoffset=(10.^[agn_gdn.lmass_candels,agn_gds.lmass_candels,agn_egs.lmass_candels]-10.^(([alog10(10.^(agn_gdn.l4p5[0])-agn_gdn.l4p5ex),alog10(10.^(agn_gds.l4p5[0])-agn_gds.l4p5ex),alog10(10.^(agn_egs.l4p5[0])-agn_egs.l4p5ex)]-linear_result[0])/linear_result[1]))/(10.^(([alog10(10.^(agn_gdn.l4p5[0])-agn_gdn.l4p5ex),alog10(10.^(agn_gds.l4p5[0])-agn_gds.l4p5ex),alog10(10.^(agn_egs.l4p5[0])-agn_egs.l4p5ex)]-linear_result[0])/linear_result[1]))
; print,max(10.^agn_gdn.lmass_candels,10.^((alog10(10.^(agn_gdn.l4p5[0])-agn_gdn.l4p5ex)-linear_result[0])/linear_result[1]))
; cghistoplot,agn_massoffset;,xran=[]
; cgplot,[agn_gdn.lmass_candels,agn_gds.lmass_candels,agn_egs.lmass_candels],([alog10(10.^(agn_gdn.l4p5[0])-agn_gdn.l4p5ex),alog10(10.^(agn_gds.l4p5[0])-agn_gds.l4p5ex),alog10(10.^(agn_egs.l4p5[0])-agn_egs.l4p5ex)]-linear_result[0])/linear_result[1],psym=4,xran=[9,12],yran=[9,12],xtitle='candels mass',ytitle='linear fit mass'
; x=findgen(91)/20.+8
; cgoplot,x,x
end
; ; error calculation
; (sqrt((10.^nonagn_gdn.l4p5[0]*alog(10)*nonagn_gdn.e4p5[0])^2+(nonagn_gdn.e4p5ex)^2))/((10.^nonagn_gdn.l4p5[0]-nonagn_gdn.l4p5ex)*alog(10))
; (sqrt((10.^nonagn_gds.l4p5[0]*alog(10)*nonagn_gds.e4p5[0])^2+(nonagn_gds.e4p5ex)^2))/((10.^nonagn_gds.l4p5[0]-nonagn_gds.l4p5ex)*alog(10))
; (sqrt((10.^nonagn_egs.l4p5[0]*alog(10)*nonagn_egs.e4p5[0])^2+(nonagn_egs.e4p5ex)^2))/((10.^nonagn_egs.l4p5[0]-nonagn_egs.l4p5ex)*alog(10))
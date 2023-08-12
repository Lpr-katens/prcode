hdu_gdn=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/goodsn_Huangall_candels_radec_van_sfr.fits',1)&$
hdu_gds=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/goodss_Huangall_candels_radec_van_sfr.fits',1)&$
hdu_egs=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/egs_Huangall_candels_radec_van_sfr.fits',1)&$
hdu_gdn=hdu_gdn[where(hdu_gdn.z_used gt 0.8 and hdu_gdn.z_used lt 1.3 and hdu_gdn.id ne -1 and hdu_gdn.separation_candels_16 lt 1)]&$
hdu_gds=hdu_gds[where(hdu_gds.z_used gt 0.8 and hdu_gds.z_used lt 1.3 and hdu_egs.id ne -1 and hdu_gds.separation_candels_16 lt 1)]&$
hdu_egs=hdu_egs[where(hdu_egs.z_used gt 0.8 and hdu_egs.z_used lt 1.3 and hdu_egs.id ne -1 and hdu_egs.separation_candels_16 lt 1)]

agn_45_gdn=hdu_gdn[where(hdu_gdn.l4p5ex gt 10.^9.5)]&$
agn_45_gds=hdu_gds[where(hdu_gds.l4p5ex gt 10.^9.5)]&$
agn_45_egs=hdu_egs[where(hdu_egs.l4p5ex gt 10.^9.5)]

; agn_45_gdn=hdu_gdn[where(hdu_gdn.l4p5ex gt 10.^9.5 and hdu_gdn.tmp_class[0] eq 1)]&$
; agn_45_gds=hdu_gds[where(hdu_gds.l4p5ex gt 10.^9.5 and hdu_gds.tmp_class[0] eq 1)]&$
; agn_45_egs=hdu_egs[where(hdu_egs.l4p5ex gt 10.^9.5 and hdu_egs.tmp_class[0] eq 1)]


l4p5ex_45=[agn_45_gdn.l4p5ex,agn_45_gds.l4p5ex,agn_45_egs.l4p5ex]&$
e4p5ex_45=[agn_45_gdn.e4p5ex,agn_45_gds.e4p5ex,agn_45_egs.e4p5ex]&$
n_45=[agn_45_gdn.n_f160w,agn_45_gds.n_f160w,agn_45_egs.n_f160w]&$
log_l4p5ex_45=alog10(l4p5ex_45)&$
log_e4p5ex_45=e4p5ex_45/(l4p5ex_45*alog(10))&$
result_45=linfit(n_45,log_l4p5ex_45,MEASURE_ERRORS=log_e4p5ex_45)&$
pearson_coeff_45=correlate(n_45,log_l4p5ex_45)

cgplot,n_45,log_l4p5ex_45,err_yhigh=log_e4p5ex_45,err_ylow=log_e4p5ex_45,color='red',title='AGN strength vs. n',xtitle='n',ytitle='!NL!S!EExc!R!I4.5',psym=7,yran=[8,12.5],xran=[-0.2,8.5],charsize=2

x=findgen(100)/100*8&$
cgoplot,x,result_45[0]+result_45[1]*x,color='black'

cgtext,1,12,'pearson correlatio: '+string(pearson_coeff_45,format='(F4.2)')
cglegend,psyms=[0,0],symcolors=['blu7','org7'],symsize=1,location=[0.25,0.75],titles=['with error','without error'],colors=['blu7','org7'],tcolors=['blu7','org7'],/box
saveimage,'/Users/lpr/Data/lirg_project/output/images/AGN_45_n.png'

agn_sed_gdn=hdu_gdn[where(hdu_gdn.tmp_class[0] eq 1)]&$
agn_sed_gds=hdu_gds[where(hdu_gds.tmp_class[0] eq 1)]&$
agn_sed_egs=hdu_egs[where(hdu_egs.tmp_class[0] eq 1)]

l4p5ex_sed=[agn_sed_gdn.l4p5ex,agn_sed_gds.l4p5ex,agn_sed_egs.l4p5ex]&$
n_sed=[agn_sed_gdn.n_f160w,agn_sed_gds.n_f160w,agn_sed_egs.n_f160w]&$
log_l4p5ex_sed=alog10(l4p5ex_sed)&$
inf_id=where(~finite(log_l4p5ex_sed))&$
remove,inf_id,log_l4p5ex_sed&$
remove,inf_id,n_sed&$
result_sed=linfit(n_sed,log_l4p5ex_sed)&$
pearson_coeff_sed=correlate(n_sed,log_l4p5ex_sed)

cgplot,n_sed,log_l4p5ex_sed,color='red',title='AGN strength vs. n',xtitle='n',ytitle='!NL!S!EExc!R!I4.5',psym=7,yran=[8,12.5],charsize=2

cgoplot,x,result_sed[0]+result_sed[1]*x,color='black'

cgtext,1,12,'pearson correlatio: '+string(pearson_coeff_sed,format='(F4.2)')

saveimage,'/Users/lpr/Data/lirg_project/output/images/AGN_sed_n.png'


; ================== with power law ====================
hdu_gdn=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/goodsn_Huangall_candels_radec_van_sfr.fits',1)&$
hdu_gds=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/goodss_Huangall_candels_radec_van_sfr.fits',1)&$
hdu_egs=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/egs_Huangall_candels_radec_van_sfr.fits',1)&$
hdu_gdn=hdu_gdn[where(hdu_gdn.z_used gt 0.8 and hdu_gdn.z_used lt 1.3 and hdu_gdn.id ne -1 and hdu_gdn.separation_candels_16 lt 1)]&$
hdu_gds=hdu_gds[where(hdu_gds.z_used gt 0.8 and hdu_gds.z_used lt 1.3 and hdu_egs.id ne -1 and hdu_gds.separation_candels_16 lt 1)]&$
hdu_egs=hdu_egs[where(hdu_egs.z_used gt 0.8 and hdu_egs.z_used lt 1.3 and hdu_egs.id ne -1 and hdu_egs.separation_candels_16 lt 1)]&$
agn_45_gdn=hdu_gdn[where(hdu_gdn.l4p5ex/hdu_gdn.e4p5ex gt 3 and hdu_gdn.tmp_class[0] eq 1 and hdu_gdn.n_f160w ne 8)]&$
agn_45_gds=hdu_gds[where(hdu_gds.l4p5ex/hdu_gds.e4p5ex gt 3 and hdu_gds.tmp_class[0] eq 1 and hdu_gds.n_f160w ne 8)]&$
agn_45_egs=hdu_egs[where(hdu_egs.l4p5ex/hdu_egs.e4p5ex gt 3 and hdu_egs.tmp_class[0] eq 1 and hdu_egs.n_f160w ne 8)]&$
l4p5ex_45=[agn_45_gdn.l4p5ex,agn_45_gds.l4p5ex,agn_45_egs.l4p5ex]&$
e4p5ex_45=[agn_45_gdn.e4p5ex,agn_45_gds.e4p5ex,agn_45_egs.e4p5ex]&$
n_45=[agn_45_gdn.n_f160w,agn_45_gds.n_f160w,agn_45_egs.n_f160w]&$
log_l4p5ex_45=alog10(l4p5ex_45)&$
log_e4p5ex_45=e4p5ex_45/(l4p5ex_45*alog(10))&$
result_45=linfit(n_45,log_l4p5ex_45,MEASURE_ERRORS=log_e4p5ex_45)&$
result_45_noerr=linfit(n_45,log_l4p5ex_45)&$
weights=1.0/e4p5ex_45&$
a=[10.^result_45[0],result_45[1]]&$
yfit=CURVEFIT(n_45,l4p5ex_45,weights,A,SIGMA,FUNCTION_NAME='powfunct')&$
cgplot,n_45,log_l4p5ex_45,err_yhigh=log_e4p5ex_45,err_ylow=log_e4p5ex_45,color='red',title='AGN strength vs. n',xtitle='n',ytitle='!NL!S!EExc!R!I4.5',psym=7,yran=[8.7,12.5],xran=[-0.2,8.5],charsize=2&$
x=findgen(100)/100*8&$
cgoplot,x,result_45[0]+result_45[1]*x,color='orange'&$
cgoplot,x,result_45_noerr[0]+result_45_noerr[1]*x,color='black'&$
cgoplot,x,alog10(a[0])+a[1]*x,color='blue'&$
cglegend,psyms=[0,0,0],symcolors=['org7','blk7','blu7'],symsize=1,location=[0.2,0.8],titles=['linear (with error)','linear (without error)','power-law'],colors=['org7','blk7','blu7'],tcolors=['org7','blk7','blu7'],/box
saveimage,'/Users/lpr/Data/lirg_project/output/images/AGN_n_pl.png'
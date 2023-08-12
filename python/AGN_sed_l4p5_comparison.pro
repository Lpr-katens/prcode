hdu_gdn=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/goodsn_Huangall_candels_radec_van_sfr.fits',1)&$
hdu_gds=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/goodss_Huangall_candels_radec_van_sfr.fits',1)&$
hdu_egs=mrdfits('/Users/lpr/Data/lirg_project/output/catalog/egs_Huangall_candels_radec_van_sfr.fits',1)&$
hdu_gdn=hdu_gdn[where(hdu_gdn.z_used gt 0.8 and hdu_gdn.z_used lt 1.3 and hdu_gdn.id ne -1 and hdu_gdn.separation_candels_16 lt 1)]&$
hdu_gds=hdu_gds[where(hdu_gds.z_used gt 0.8 and hdu_gds.z_used lt 1.3 and hdu_egs.id ne -1 and hdu_gds.separation_candels_16 lt 1)]&$
hdu_egs=hdu_egs[where(hdu_egs.z_used gt 0.8 and hdu_egs.z_used lt 1.3 and hdu_egs.id ne -1 and hdu_egs.separation_candels_16 lt 1)]

agn_sed_gdn=hdu_gdn[where(hdu_gdn.tmp_class[0] eq 1)]&$
agn_sed_gds=hdu_gds[where(hdu_gds.tmp_class[0] eq 1)]&$
agn_sed_egs=hdu_egs[where(hdu_egs.tmp_class[0] eq 1)]&$
agn_45_gdn=hdu_gdn[where(hdu_gdn.l4p5ex gt 10.^9.5)]&$
agn_45_gds=hdu_gds[where(hdu_gds.l4p5ex gt 10.^9.5)]&$
agn_45_egs=hdu_egs[where(hdu_egs.l4p5ex gt 10.^9.5)]&$
gdn_overlap_sedbased=agn_sed_gdn[where(agn_sed_gdn.l4p5ex gt 10.^9.5)]&$ 
gds_overlap_sedbased=agn_sed_gds[where(agn_sed_gds.l4p5ex gt 10.^9.5)]&$
egs_overlap_sedbased=agn_sed_egs[where(agn_sed_egs.l4p5ex gt 10.^9.5)]&$
gdn_overlap_45based=agn_45_gdn[where(agn_45_gdn.tmp_class[0] eq 1)]&$
gds_overlap_45based=agn_45_gds[where(agn_45_gds.tmp_class[0] eq 1)]&$
egs_overlap_45based=agn_45_egs[where(agn_45_egs.tmp_class[0] eq 1)]

help,gdn_overlap_sedbased,agn_sed_gdn&$
help,gds_overlap_sedbased,agn_sed_gds&$
help,egs_overlap_sedbased,agn_sed_egs
help,gdn_overlap_45based,agn_45_gdn&$
help,gds_overlap_45based,agn_45_gds&$
help,egs_overlap_45based,agn_45_egs

window,0
cghistoplot,agn_sed_gdn.n_f160w,color='red',mininput=0,maxinput=8,binsize=1,title='gdn-n',xtitle='n'&$ 
cghistoplot,agn_45_gdn.n_f160w,color='green',mininput=0,maxinput=8,binsize=1,/oplot&$
cghistoplot,gdn_overlap_sedbased.n_f160w,color='blue',mininput=0,maxinput=8,binsize=1,/oplot&$
cglegend,psyms=[4,4,4],symcolors=['red7','blu7','grn7'],symsize=1,location=[0.7,0.7],titles=['sed','overlap','l4p5ex']&$
saveimage,'/Users/lpr/Data/lirg_project/output/images/sed-45-n.png'

cghistoplot,agn_45_gdn.re_f160w,color='green',mininput=0,maxinput=2,nbins=6,/oplot&$
cghistoplot,gdn_overlap_sedbased.re_f160w,color='blue',mininput=0,maxinput=2,nbins=6,/oplot&$
cglegend,psyms=[4,4,4],symcolors=['red7','blu7','grn7'],symsize=1,location=[0.7,0.7],titles=['sed','overlap','l4p5ex']&$
saveimage,'/Users/lpr/Data/lirg_project/output/images/sed-45-re.png'
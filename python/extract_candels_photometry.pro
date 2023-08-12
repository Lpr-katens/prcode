fields = ['goodsn','goodss','egs']
match_catalog = '/Users/lpr/Data/lirg_project/output/catalog/'
match_catalog_suffix = '_Huangall_radec_candels.fits'
candels_catalog = '/Users/lpr/Data/lirg_project/intake/CANDELS/catalog/JFang_CANDELS_Data/'
inner = {goodsn:['_APER_3','_APER_10'],goodss:['_APER_4','_APER_11'],egs:['_APER_4','_APER_11']}
filter_list = ['F435W','F606W','F775W','F814W','F850LP','F098M','F105W','F125W','F140W','F160W']
for num1 = 0, 2 do begin
match_hdu = mrdfits(match_catalog+fields[num1]+match_catalog_suffix,1)
if fields[num1] eq 'egs' then begin
candels_hdu = mrdfits(candels_catalog+fields[num1]+'_all.fits',1)
endif else begin
candels_hdu = mrdfits(candels_catalog+strmid(fields[num1],0,1)+strmid(fields[num1],3,1)+strmid(fields[num1],5,1)+'_all.fits',1)
endelse
column_list = tag_names(candels_hdu)
; make sure which filter this field has
if fields[num1] eq 'egs' then begin
array = {ID_Huang:0ll,ID_CANDELS:0ll,F606W_APER_4:-999.d,F606W_APER_11:-999.d,F814W_APER_4:-999.d,F814W_APER_11:-999.d,F125W_APER_4:-999.d,F125W_APER_11:-999.d,F140W_APER_4:-999.d,F140W_APER_11:-999.d,F160W_APER_4:-999.d,F160W_APER_11:-999.d}
endif else begin
if fields[num1] eq 'goodss' then begin
array = {ID_Huang:0ll,ID_CANDELS:0ll,F435W_APER_4:-999.d,F435W_APER_11:-999.d,F606W_APER_4:-999.d,F606W_APER_11:-999.d,F775W_APER_4:-999.d,F775W_APER_11:-999.d,F814W_APER_4:-999.d,F814W_APER_11:-999.d,F850LP_APER_4:-999.d,F850LP_APER_11:-999.d,F098M_APER_4:-999.d,F098M_APER_11:-999.d,F105W_APER_4:-999.d,F105W_APER_11:-999.d,F125W_APER_4:-999.d,F125W_APER_11:-999.d,F160W_APER_4:-999.d,F160W_APER_11:-999.d}
endif
if fields[num1] eq 'goodsn' then begin
array = {ID_Huang:0ll,ID_CANDELS:0ll,F435W_APER_4:-999.d,F435W_APER_11:-999.d,F606W_APER_4:-999.d,F606W_APER_11:-999.d,F775W_APER_4:-999.d,F775W_APER_11:-999.d,F814W_APER_4:-999.d,F814W_APER_11:-999.d,F850LP_APER_4:-999.d,F850LP_APER_11:-999.d,F105W_APER_4:-999.d,F105W_APER_11:-999.d,F125W_APER_4:-999.d,F125W_APER_11:-999.d,F160W_APER_4:-999.d,F160W_APER_11:-999.d,ACS_F435W_FLUX:-999.d,ACS_F606W_FLUX:-999.d,ACS_F775W_FLUX:-999.d,ACS_F814W_FLUX:-999.d,ACS_F850LP_FLUX:-999.d,WFC3_F105W_FLUX:-999.d,WFC3_F125W_FLUX:-999.d,WFC3_F160W_FLUX:-999.d}
endif
endelse
print,'---------------------------------------- loading initial array ----------------------------------------' + fields[num1]
array = replicate(array,size(match_hdu,/n_elements))
id_huang_col = make_array(size(match_hdu,/n_elements),/L64,VALUE=0)
id_candels_col = make_array(size(match_hdu,/n_elements),/L64,VALUE=0)
if fields[num1] eq 'egs' then begin
f606w_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f606w_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f814w_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f814w_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f125w_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f125w_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f140w_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f140w_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f160w_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f160w_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
endif else begin
f435w_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f435w_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f606w_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f606w_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f775w_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f775w_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f814w_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f814w_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f850lp_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f850lp_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
if fields[num1] eq 'goodss' then begin
f098m_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f098m_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
endif
f105w_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f105w_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f125w_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f125w_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f160w_4 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
f160w_11 = make_array(size(match_hdu,/n_elements),/double,value=-999.)
if fields[num1] eq 'goodsn' then begin
tphot_f435w = make_array(size(match_hdu,/n_elements),/double,value=-999.)
tphot_f606w = make_array(size(match_hdu,/n_elements),/double,value=-999.)
tphot_f775w = make_array(size(match_hdu,/n_elements),/double,value=-999.)
tphot_f814w = make_array(size(match_hdu,/n_elements),/double,value=-999.)
tphot_f850lp = make_array(size(match_hdu,/n_elements),/double,value=-999.)
tphot_f105w = make_array(size(match_hdu,/n_elements),/double,value=-999.)
tphot_f125w = make_array(size(match_hdu,/n_elements),/double,value=-999.)
tphot_f160w = make_array(size(match_hdu,/n_elements),/double,value=-999.)
endif
endelse
print,'---------------------------------------- extrating ----------------------------------------' + fields[num1]
for num2 = 0, size(match_hdu,/n_elements)-1 do begin
idx = where(candels_hdu.id eq match_hdu[num2].id_candels)
id_huang_col[num2] = match_hdu[num2].id
id_candels_col[num2] = match_hdu[num2].id_candels
if fields[num1] eq 'egs' then begin
f606w_4[num2] = candels_hdu[idx].flux_aper_4_f606w
f606w_11[num2] = candels_hdu[idx].flux_aper_11_f606w
f814w_4[num2] = candels_hdu[idx].flux_aper_4_f814w
f814w_11[num2] = candels_hdu[idx].flux_aper_11_f814w
f125w_4[num2] = candels_hdu[idx].flux_aper_4_f125w
f125w_11[num2] = candels_hdu[idx].flux_aper_11_f125w
f140w_4[num2] = candels_hdu[idx].flux_aper_4_f140w
f140w_11[num2] = candels_hdu[idx].flux_aper_11_f140w
f160w_4[num2] = candels_hdu[idx].flux_aper_4_f160w
f160w_11[num2] = candels_hdu[idx].flux_aper_11_f160w
endif else begin
; goodsn aperture is not the same as egs and goodss
if fields[num1] eq 'goodsn' then begin
f435w_4[num2] = candels_hdu[idx].flux_aper_3_f435w
f435w_11[num2] = candels_hdu[idx].flux_aper_10_f435w
f606w_4[num2] = candels_hdu[idx].flux_aper_3_f606w
f606w_11[num2] = candels_hdu[idx].flux_aper_10_f606w
f775w_4[num2] = candels_hdu[idx].flux_aper_3_f775w
f775w_11[num2] = candels_hdu[idx].flux_aper_10_f775w
f814w_4[num2] = candels_hdu[idx].flux_aper_3_f814w
f814w_11[num2] = candels_hdu[idx].flux_aper_10_f814w
f850lp_4[num2] = candels_hdu[idx].flux_aper_3_f850lp
f850lp_11[num2] = candels_hdu[idx].flux_aper_10_f850lp
f105w_4[num2] = candels_hdu[idx].flux_aper_3_f105w
f105w_11[num2] = candels_hdu[idx].flux_aper_10_f105w
f125w_4[num2] = candels_hdu[idx].flux_aper_3_f125w
f125w_11[num2] = candels_hdu[idx].flux_aper_10_f125w
f160w_4[num2] = candels_hdu[idx].flux_aper_3_f160w
f160w_11[num2] = candels_hdu[idx].flux_aper_10_f160w
tphot_f435w[num2] = candels_hdu[idx].acs_f435w_flux
tphot_f606w[num2] = candels_hdu[idx].acs_f606w_flux
tphot_f775w[num2] = candels_hdu[idx].acs_f775w_flux
tphot_f814w[num2] = candels_hdu[idx].acs_f814w_flux
tphot_f850lp[num2] = candels_hdu[idx].acs_f850lp_flux
tphot_f105w[num2] = candels_hdu[idx].wfc3_f105w_flux
tphot_f125w[num2] = candels_hdu[idx].wfc3_f125w_flux
tphot_f160w[num2] = candels_hdu[idx].wfc3_f160w_flux
endif
if fields[num1] eq 'goodss' then begin
f435w_4[num2] = candels_hdu[idx].flux_aper_4_f435w
f435w_11[num2] = candels_hdu[idx].flux_aper_11_f435w
f606w_4[num2] = candels_hdu[idx].flux_aper_4_f606w
f606w_11[num2] = candels_hdu[idx].flux_aper_11_f606w
f775w_4[num2] = candels_hdu[idx].flux_aper_4_f775w
f775w_11[num2] = candels_hdu[idx].flux_aper_11_f775w
f814w_4[num2] = candels_hdu[idx].flux_aper_4_f814w
f814w_11[num2] = candels_hdu[idx].flux_aper_11_f814w
f850lp_4[num2] = candels_hdu[idx].flux_aper_4_f850lp
f850lp_11[num2] = candels_hdu[idx].flux_aper_11_f850lp
f098m_4[num2] = candels_hdu[idx].flux_aper_4_f098m
f098m_11[num2] = candels_hdu[idx].flux_aper_11_f098m
f105w_4[num2] = candels_hdu[idx].flux_aper_4_f105w
f105w_11[num2] = candels_hdu[idx].flux_aper_11_f105w
f125w_4[num2] = candels_hdu[idx].flux_aper_4_f125w
f125w_11[num2] = candels_hdu[idx].flux_aper_11_f125w
f160w_4[num2] = candels_hdu[idx].flux_aper_4_f160w
f160w_11[num2] = candels_hdu[idx].flux_aper_11_f160w
endif
endelse
endfor
array.id_huang = id_huang_col
array.id_candels = id_candels_col
if fields[num1] eq 'egs' then begin
array.F606W_aper_4 = f606w_4
array.F606W_aper_11 = f606w_11
array.F814W_aper_4 = f814w_4
array.F814W_aper_11 = f814w_11
array.F125W_aper_4 = f125w_4
array.F125W_aper_11 = f125w_11
array.F140W_aper_4 = f140w_4
array.F140W_aper_11 = f140w_11
array.F160W_aper_4 = f160w_4
array.F160W_aper_11 = f160w_11
endif else begin

array.F435W_aper_4 = f435w_4
array.F435W_aper_11 = f435w_11
array.F606W_aper_4 = f606w_4
array.F606W_aper_11 = f606w_11
array.F775W_aper_4 = f775w_4
array.F775W_aper_11 = f775w_11
array.F814W_aper_4 = f814w_4
array.F814W_aper_11 = f814w_11
array.F850LP_aper_4 = f850lp_4
array.F850LP_aper_11 = f850lp_11
if fields[num1] eq 'goodss' then begin
array.F098M_aper_4 = f098m_4
array.F098M_aper_11 = f098m_11
endif
array.F105W_aper_4 = f105w_4
array.F105W_aper_11 = f105w_11
array.F125W_aper_4 = f125w_4
array.F125W_aper_11 = f125w_11
array.F160W_aper_4 = f160w_4
array.F160W_aper_11 = f160w_11
if fields[num1] eq 'goodsn' then begin
array.acs_f435w_flux = tphot_f435w
array.acs_f606w_flux = tphot_f606w
array.acs_f775w_flux = tphot_f775w
array.acs_f814w_flux = tphot_f814w
array.acs_f850lp_flux = tphot_f850lp
array.wfc3_f105w_flux = tphot_f105w
array.wfc3_f125w_flux = tphot_f125w
array.wfc3_f160w_flux = tphot_f160w
endif
endelse
print,'---------------------------------------- writing array ----------------------------------------' + fields[num1]
mwrfits,array,match_catalog+fields[num1]+'_extracted_candels_photometry.fits';
endfor
end
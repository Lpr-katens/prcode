; hc = 6.62607*1e-27*3*1e10
; sed_path = '/Users/lpr/Data/fits/pridata/BrownSEDtemplates/Brown_sed'
; rest_filter_path = '/Users/lpr/Data/fits/pridata/filter/rest_frame/'
; obs_filter_path = '/Users/lpr/Data/fits/pridata/filter/obs_frame/'
; sed_file = file_search(sed_path,'*.dat',count=sed_num)
; res_filter_file = file_search(rest_filter_path,'*.dat',count=rest_filter_num)
; obs_filter_file = file_search(obs_filter_path,'*.dat',count=obs_filter_num)
; galaxy_file = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/color_difference/uvi_color_obs.fits'
; zmin = 0.8
; zmax = 1.3
; bin = 100
; export_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/color_difference/'
; name = export_path + 'UV_606125.ps'
; readcol,rest_filter_path+'Generic_Bessell.U.dat',trans_lambda_u,transmittance_u,format='f',/sil
; readcol,rest_filter_path+'Generic_Bessell.V.dat',trans_lambda_v,transmittance_v,format='f',/sil
; readcol,obs_filter_path+'HST_ACS_WFC.F606W.dat',trans_lambda_606,transmittance_606,format='f',/sil
; readcol,obs_filter_path+'HST_WFC3_IR.F125W.dat',trans_lambda_125,transmittance_125,format='f',/sil
; y_axis = dblarr(sed_num)
; x_axis = dblarr(sed_num)
; z = 0.8
; set_plot, 'ps'
; device, filename=name
; for num1 = 0 , sed_num-1 do begin
; 	readcol,sed_file[num1],sed_lambda,sed_flux,format='f',/sil
; 	interpolated_u = interpol(sed_flux,sed_lambda,trans_lambda_u)
; 	interpolated_v = interpol(sed_flux,sed_lambda,trans_lambda_v)
; 	convolve_u = int_tabulated(trans_lambda_u,interpolated_u*(hc*trans_lambda_u)^(-1)*transmittance_u)
; 	convolve_v = int_tabulated(trans_lambda_v,interpolated_v*(hc*trans_lambda_v)^(-1)*transmittance_v)
; 	convolve_denominator_u = int_tabulated(trans_lambda_u,transmittance_u*3631*1e-23*(hc*trans_lambda_u)^(-1)/(trans_lambda_u^2))
; 	convolve_denominator_v = int_tabulated(trans_lambda_v,transmittance_v*3631*1e-23*(hc*trans_lambda_v)^(-1)/(trans_lambda_v^2))
; 	mag_z_u = -2.5*alog10(convolve_u/convolve_denominator_u)
; 	mag_z_v = -2.5*alog10(convolve_v/convolve_denominator_v)
; 	color_uv = mag_z_u - mag_z_v
; 	y_axis[num1] = color_uv
; 	interpolated_606 = interpol(sed_flux,sed_lambda*(1+z),trans_lambda_606)
; 	interpolated_125 = interpol(sed_flux,sed_lambda*(1+z),trans_lambda_125)
; 	convolve_606 = int_tabulated(trans_lambda_606,interpolated_606*(hc*trans_lambda_606)^(-1)*transmittance_606)
; 	convolve_125 = int_tabulated(trans_lambda_125,interpolated_125*(hc*trans_lambda_125)^(-1)*transmittance_125)
; 	convolve_denominator_606 = int_tabulated(trans_lambda_606,transmittance_606*3631*1e-23*(hc*trans_lambda_606)^(-1)/(trans_lambda_606^2))
; 	convolve_denominator_125 = int_tabulated(trans_lambda_125,transmittance_125*3631*1e-23*(hc*trans_lambda_125)^(-1)/(trans_lambda_125^2))
; 	mag_z_606 = -2.5*alog10(convolve_606/convolve_denominator_606)
; 	mag_z_125 = -2.5*alog10(convolve_125/convolve_denominator_125)
; 	color_606125 = mag_z_606 - mag_z_125
; 	x_axis[num1] = color_606125
; endfor
; plot,x_axis,y_axis,psym=2,xtitle='F606W - F125W',ytitle='U - V',title='K-correction' ;,color=cgcolor(string(color[num1]))
; device, /close_file
; end
;Here to calculate the UVI diagram
; hc = 6.62607*1e-27*3*1e10
; sed_path = '/Users/lpr/Data/fits/pridata/BrownSEDtemplates/Brown_sed'
; rest_filter_path = '/Users/lpr/Data/fits/pridata/filter/rest_frame/'
; obs_filter_path = '/Users/lpr/Data/fits/pridata/filter/obs_frame/'
; export_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/color_difference/'
; sed_file = file_search(sed_path,'*.dat',count=sed_num)
; res_filter_file = file_search(rest_filter_path,'*.dat',count=rest_filter_num)
; obs_filter_file = file_search(obs_filter_path,'*.dat',count=obs_filter_num)
; galaxy_file = mrdfits('/Users/lpr/Data/fits/expdata/HST/goodsn_all/color_difference/uvi_color_obs.fits',1)
; ; aper_color_cata = mrdfits('/Users/lpr/Desktop/goodsn_restcolor.fits',1)
; galaxy_info = mrdfits('/Users/lpr/Data/fits/expdata/HST/goodsn_all/goodsn_Huangall_radecmatch_modifyz.fits',1)
; readcol,rest_filter_path+'Generic_Bessell.U.dat',bessell_lambda_u,transmittance_u_bessell,format='f',/sil
; readcol,rest_filter_path+'Generic_Bessell.V.dat',bessell_lambda_v,transmittance_v_bessell,format='f',/sil
; readcol,rest_filter_path+'Generic_Bessell.I.dat',bessell_lambda_i,transmittance_i_bessell,format='f',/sil
; array = {ID_Huang:0ll,UV_bulge:-999.d,UV_disk:-999.d,VI_bulge:0.d,VI_disk:0.d}
; array = replicate(array,size(galaxy_file,/n_elements))
; id_col = MAKE_ARRAY(size(galaxy_file,/n_elements),/L64,VALUE=0)
; uv_bulge_res = dblarr(size(galaxy_file,/n_elements))
; vi_bulge_res = dblarr(size(galaxy_file,/n_elements))
; uv_disk_res = dblarr(size(galaxy_file,/n_elements))
; vi_disk_res = dblarr(size(galaxy_file,/n_elements))
; for num2 = 0 , size(galaxy_file,/n_elements)-1 do begin
; 	idx = galaxy_file[num2].id_huang
; 	if where(galaxy_info.id eq idx) ne -1 then begin
; 		z = galaxy_info[where(galaxy_info.id eq idx)].z_used
; 		u_hst_filter = galaxy_file[num2].u_hst_filter
; 		v_hst_filter = galaxy_file[num2].v_hst_filter
; 		i_hst_filter = galaxy_file[num2].i_hst_filter
; 		; uv_bulge_obs = aper_color_cata[where(aper_color_cata.id_huang eq idx)].uv_bulge
; 		; vi_bulge_obs = aper_color_cata[where(aper_color_cata.id_huang eq idx)].vi_bulge
; 		; uv_disk_obs = aper_color_cata[where(aper_color_cata.id_huang eq idx)].uv_disk
; 		; vi_disk_obs = aper_color_cata[where(aper_color_cata.id_huang eq idx)].vi_disk
; 		uv_bulge_obs = galaxy_file[num2].uv_bulge
; 		vi_bulge_obs = galaxy_file[num2].vi_bulge
; 		uv_disk_obs = galaxy_file[num2].uv_disk
; 		vi_disk_obs = galaxy_file[num2].vi_disk
; 		; here to open hst filter, which this galaxy used
; 		for num3 = 0 , obs_filter_num-1 do begin
; 			if strpos(obs_filter_file[num3],strtrim(u_hst_filter,1)) ne -1 then begin
; 				readcol,obs_filter_file[num3],trans_lambda_u_obs,transmittance_u_obs,format='f'
; 			endif else begin
; 				if strpos(obs_filter_file[num3],strtrim(v_hst_filter,1)) ne -1 then begin
; 					readcol,obs_filter_file[num3],trans_lambda_v_obs,transmittance_v_obs,format='f'
; 				endif else begin
; 					if strpos(obs_filter_file[num3],strtrim(i_hst_filter,1)) ne -1 then begin
; 						readcol,obs_filter_file[num3],trans_lambda_i_obs,transmittance_i_obs,format='f'
; 					endif
; 				endelse
; 			endelse
; 		endfor
; 		; here comes the loop of every sed template
; 		uv_res = dblarr(sed_num)
; 		vi_res = dblarr(sed_num)
; 		uv_obs = dblarr(sed_num)
; 		vi_obs = dblarr(sed_num)
; 		for num1 = 0 , sed_num-1 do begin
; 			readcol,sed_file[num1],sed_lambda,sed_flux,format='f',/sil
; 			interpolated_u = interpol(sed_flux,sed_lambda,bessell_lambda_u)
; 			interpolated_v = interpol(sed_flux,sed_lambda,bessell_lambda_v)
; 			interpolated_i = interpol(sed_flux,sed_lambda,bessell_lambda_i)
; 			convolve_u = int_tabulated(bessell_lambda_u,interpolated_u*(hc*bessell_lambda_u)^(-1)*transmittance_u_bessell)
; 			convolve_v = int_tabulated(bessell_lambda_v,interpolated_v*(hc*bessell_lambda_v)^(-1)*transmittance_v_bessell)
; 			convolve_i = int_tabulated(bessell_lambda_i,interpolated_i*(hc*bessell_lambda_i)^(-1)*transmittance_i_bessell)
; 			convolve_denominator_u = int_tabulated(bessell_lambda_u,transmittance_u_bessell*3631*1e-23*(hc*bessell_lambda_u)^(-1)/(bessell_lambda_u^2))
; 			convolve_denominator_v = int_tabulated(bessell_lambda_v,transmittance_v_bessell*3631*1e-23*(hc*bessell_lambda_v)^(-1)/(bessell_lambda_v^2))
; 			convolve_denominator_i = int_tabulated(bessell_lambda_i,transmittance_u_bessell*3631*1e-23*(hc*bessell_lambda_i)^(-1)/(bessell_lambda_i^2))
; 			mag_u = -2.5*alog10(convolve_u/convolve_denominator_u)
; 			mag_v = -2.5*alog10(convolve_v/convolve_denominator_v)
; 			mag_i = -2.5*alog10(convolve_i/convolve_denominator_i)
; 			uv_res[num1] = mag_u - mag_v
; 			vi_res[num1] = mag_v - mag_i
; 			interpolated_u_obs = interpol(sed_flux,sed_lambda*(1+z),trans_lambda_u_obs)
; 			interpolated_v_obs = interpol(sed_flux,sed_lambda*(1+z),trans_lambda_v_obs)
; 			interpolated_i_obs = interpol(sed_flux,sed_lambda*(1+z),trans_lambda_i_obs)
; 			convolve_u_obs = int_tabulated(trans_lambda_u_obs,interpolated_u_obs*(hc*trans_lambda_u_obs)^(-1)*transmittance_u_obs)
; 			convolve_v_obs = int_tabulated(trans_lambda_v_obs,interpolated_v_obs*(hc*trans_lambda_v_obs)^(-1)*transmittance_v_obs)
; 			convolve_i_obs = int_tabulated(trans_lambda_i_obs,interpolated_i_obs*(hc*trans_lambda_i_obs)^(-1)*transmittance_i_obs)
; 			convolve_denominator_u_obs = int_tabulated(trans_lambda_u_obs,transmittance_u_obs*3631*1e-23*(hc*trans_lambda_u_obs)^(-1)/(trans_lambda_u_obs^2))
; 			convolve_denominator_v_obs = int_tabulated(trans_lambda_v_obs,transmittance_v_obs*3631*1e-23*(hc*trans_lambda_v_obs)^(-1)/(trans_lambda_v_obs^2))
; 			convolve_denominator_i_obs = int_tabulated(trans_lambda_i_obs,transmittance_i_obs*3631*1e-23*(hc*trans_lambda_i_obs)^(-1)/(trans_lambda_i_obs^2))
; 			mag_u_obs = -2.5*alog10(convolve_u_obs/convolve_denominator_u_obs)
; 			mag_v_obs = -2.5*alog10(convolve_v_obs/convolve_denominator_v_obs)
; 			mag_i_obs = -2.5*alog10(convolve_i_obs/convolve_denominator_i_obs)
; 			uv_obs[num1] = mag_u_obs - mag_v_obs
; 			vi_obs[num1] = mag_v_obs - mag_i_obs
; 		endfor
; 		id_col[num2] = idx
; 		uv_bulge_res[num2] = interpol(uv_res,uv_obs,uv_bulge_obs)
; 		vi_bulge_res[num2] = interpol(vi_res,vi_obs,vi_bulge_obs)
; 		uv_disk_res[num2] = interpol(uv_res,uv_obs,uv_disk_obs)
; 		vi_disk_res[num2] = interpol(vi_res,vi_obs,vi_disk_obs)
; 	endif
; endfor
; array.id_huang = id_col
; array.UV_bulge = uv_bulge_res
; array.UV_disk = uv_disk_res
; array.VI_bulge = vi_bulge_res
; array.VI_disk = vi_disk_res
; mwrfits,array,export_path+'uvi_color_res.fits';_CANDELS_aper
; end

; ; sed_path = '/Users/lpr/Data/fits/pridata/BrownSEDtemplates/Brown_sed'
; ; rest_filter_path = '/Users/lpr/Data/fits/pridata/filter/rest_frame/'
; ; obs_filter_path = '/Users/lpr/Data/fits/pridata/filter/obs_frame/'
; export_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/uvi/'
; ; sed_file = file_search(sed_path,'*.dat',count=sed_num)
; ; res_filter_file = file_search(rest_filter_path,'*.dat',count=rest_filter_num)
; ; obs_filter_file = file_search(obs_filter_path,'*.dat',count=obs_filter_num)
; galaxy_file = mrdfits('/Users/lpr/Data/fits/expdata/HST/goodsn_all/uvi/uvi_obs_CANDELS.fits',1,/sil); galaxy_file = mrdfits('/Users/lpr/Data/fits/expdata/HST/goodsn_all/uvi/uvi_obs.fits',1,/sil)
; ; aper_color_cata = mrdfits('/Users/lpr/Desktop/goodsn_restcolor.fits',1)
; galaxy_info = mrdfits('/Users/lpr/Data/fits/expdata/HST/goodsn_all/goodsn_Huangall_radecmatch_modifyz.fits',1,/sil)
; res_obs_rela_left1 = file_search('/Users/lpr/Data/fits/expdata/twocolor_kcc/left1','*.fits',count=rela_left1_num)
; res_obs_rela_right1 = file_search('/Users/lpr/Data/fits/expdata/twocolor_kcc/right1','*.fits',count=rela_right1_num)
; rela_coef_left1 = {UV:-999.d,VI:-999.d};dblarr(2,2)
; rela_coef_right1 = {UV:-999.d,VI:-999.d};dblarr(3,2)
; rela_coef_left1 = replicate(rela_coef_left1,2)
; rela_coef_right1 = replicate(rela_coef_right1,3)
; for num4 = 0 , rela_left1_num - 1 do begin
; if strmid(res_obs_rela_left1[num4],48,2) eq 'UV' then begin
; 	temp = mrdfits(res_obs_rela_left1[num4],1,/sil)
; 	ind = where(temp.obs_color gt 0)
; 	result = linfit(temp[ind].obs_color,temp[ind].res_color)
; 	rela_coef_left1.UV = double(result)
; endif else begin
; 	temp = mrdfits(res_obs_rela_left1[num4],1,/sil)
; 	ind = where(temp.obs_color gt 0)
; 	result = linfit(temp[ind].obs_color,temp[ind].res_color)
; 	rela_coef_left1.VI = double(result)
; endelse
; endfor
; for num4 = 0 , rela_right1_num - 1 do begin
; if strmid(res_obs_rela_right1[num4],49,2) eq 'UV' then begin
; 	temp = mrdfits(res_obs_rela_right1[num4],1,/sil)
; 	ind = where(temp.obs_color gt 0)
; 	result = poly_fit(temp[ind].obs_color,temp[ind].res_color,2)
; 	rela_coef_right1.UV = double(result[0:2])
; endif else begin
; 	temp = mrdfits(res_obs_rela_right1[num4],1,/sil)
; 	ind = where(temp.obs_color gt 0)
; 	result = poly_fit(temp[ind].obs_color,temp[ind].res_color,2)
; 	rela_coef_right1.VI = double(result[0:2])
; endelse
; endfor
; array = {ID_Huang:0ll,UV_bulge:-999.d,UV_disk:-999.d,VI_bulge:0.d,VI_disk:0.d,z_used:-999.d}
; array = replicate(array,size(galaxy_file,/n_elements))
; id_col = MAKE_ARRAY(size(galaxy_file,/n_elements),/L64,VALUE=0)
; uv_bulge_res = dblarr(size(galaxy_file,/n_elements))
; vi_bulge_res = dblarr(size(galaxy_file,/n_elements))
; uv_disk_res = dblarr(size(galaxy_file,/n_elements))
; vi_disk_res = dblarr(size(galaxy_file,/n_elements))
; z_used = dblarr(size(galaxy_file,/n_elements))
; for num2 = 0 , size(galaxy_file,/n_elements)-1 do begin
; 	idx = galaxy_file[num2].id_huang
; 	if where(galaxy_info.id eq idx) ne -1 then begin
; 		z = galaxy_info[where(galaxy_info.id eq idx)].z_used
; 		uv_bulge_obs = galaxy_file[num2].uv_bulge_CANDELS_obs
; 		vi_bulge_obs = galaxy_file[num2].vi_bulge_CANDELS_obs
; 		uv_disk_obs = galaxy_file[num2].uv_disk_CANDELS_obs
; 		vi_disk_obs = galaxy_file[num2].vi_disk_CANDELS_obs
; 		if z lt 1 then begin
; 			id_col[num2] = idx
; 			; print,uv_bulge_obs,uv_disk_obs
; 			uv_bulge_res[num2] = (rela_coef_left1.UV)[0] + (rela_coef_left1.UV)[1] * uv_bulge_obs
; 			vi_bulge_res[num2] = (rela_coef_left1.VI)[0] + (rela_coef_left1.VI)[1] * vi_bulge_obs
; 			uv_disk_res[num2] = (rela_coef_left1.UV)[0] + (rela_coef_left1.UV)[1] * uv_disk_obs
; 			vi_disk_res[num2] = (rela_coef_left1.VI)[0] + (rela_coef_left1.VI)[1] * vi_disk_obs
; 			z_used[num2] = z
; 		endif else begin
; 			id_col[num2] = idx
; 			; print,uv_bulge_obs,uv_disk_obs
; 			uv_bulge_res[num2] = (rela_coef_right1.UV)[0] + (rela_coef_right1.UV)[1] * uv_bulge_obs + (rela_coef_right1.UV)[2] * uv_bulge_obs^2
; 			vi_bulge_res[num2] = (rela_coef_right1.VI)[0] + (rela_coef_right1.VI)[1] * vi_bulge_obs + (rela_coef_right1.VI)[2] * vi_bulge_obs^2
; 			uv_disk_res[num2] = (rela_coef_right1.UV)[0] + (rela_coef_right1.UV)[1] * uv_disk_obs + (rela_coef_right1.UV)[2] * uv_disk_obs^2
; 			vi_disk_res[num2] = (rela_coef_right1.VI)[0] + (rela_coef_right1.VI)[1] * vi_disk_obs + (rela_coef_right1.VI)[2] * vi_disk_obs^2
; 			z_used[num2] = z
; 		endelse
; 	endif
; endfor
; array.id_huang = id_col
; array.UV_bulge = uv_bulge_res
; array.UV_disk = uv_disk_res
; array.VI_bulge = vi_bulge_res
; array.VI_disk = vi_disk_res
; array.z_used = z_used
; mwrfits,array,export_path+'uvi_res_fit_CANDELS.fits';_aper
; end

print,'         *              *'
print,'        ***            ***'
print,'       *****          *****'
print,'      **********************'
print,'    **************************'
print,'   **********  code  **********'
print,'   *********   begin  *********'
print,'    ********   now    ********'
print,'      **********************             ***'
print,'        ******************           **********'
print,'          **************        ******************'
print,'            **********      **************************'
print,'            *********************************************'
print,'          ************************************************'
print,'        ***************************************************'
print,'      ******************************************************'
print,'      *******      *******            *******     *******'
print,'       *****        *****              *****       *****'
print,'        ***          ***                ***         ***'
print,'         *            *                  *           *'
fields_list = ['goodsn','goodss','egs']
twocol_relation_path = '/Users/lpr/Data/lirg_project/output/twocolor_kcc/f_nu/'
export_path = '/Users/lpr/Data/lirg_project/output/catalog/'
observe_suffix = '_uvi_obs_CANDELS_onlygdnfrac.fits'
; sed_path = '/Users/lpr/Data/Brown/Brown_sed/'
; rest_filter_path = '/Users/lpr/Data/fits/pridata/filter/rest_frame/'
; obs_filter_path = '/Users/lpr/Data/fits/pridata/filter/obs_frame/'
; sed_file = file_search(sed_path,'*.dat',count=sed_num)
; res_filter_file = file_search(rest_filter_path,'*.dat',count=rest_filter_num)
; obs_filter_file = file_search(obs_filter_path,'*.dat',count=obs_filter_num)
for field =0,size(fields_list,/n_elements)-1 do begin
	galaxy_file = mrdfits(export_path+fields_list[field]+observe_suffix,1,/silent)
	relation = file_search(twocol_relation_path+fields_list[field],'*.fits',count=relation_counts)
	print,n_elements(galaxy_file)
	; res_obs_rela_left1 = file_search('/Users/lpr/Data/lirg_project/output/twocolor_kcc/f_nu/left1','*.fits',count=rela_left1_num)
	; res_obs_rela_right1 = file_search('/Users/lpr/Data/lirg_project/output/twocolor_kcc/f_nu/right1','*.fits',count=rela_right1_num)
	; for num4 = 0 , rela_left1_num - 1 do begin
	; if strmid(res_obs_rela_left1[num4],48,2) eq 'UV' then begin
	; 	temp = mrdfits(res_obs_rela_left1[num4],1,/silent)
	; 	ind = where(temp.obs_color gt 0)
	; 	result = linfit(temp[ind].obs_color,temp[ind].res_color)
	; 	rela_coef_left1.UV = double(result)
	; endif else begin
	; 	temp = mrdfits(res_obs_rela_left1[num4],1,/silent)
	; 	ind = where(temp.obs_color gt 0)
	; 	result = linfit(temp[ind].obs_color,temp[ind].res_color)
	; 	rela_coef_left1.VI = double(result)
	; endelse
	; endfor
	; for num4 = 0 , rela_right1_num - 1 do begin
	; if strmid(res_obs_rela_right1[num4],49,2) eq 'UV' then begin
	; 	temp = mrdfits(res_obs_rela_right1[num4],1,/silent)
	; 	ind = where(temp.obs_color gt 0)
	; 	result = poly_fit(temp[ind].obs_color,temp[ind].res_color,2)
	; 	rela_coef_right1.UV = double(result[0:2])
	; endif else begin
	; 	temp = mrdfits(res_obs_rela_right1[num4],1,/silent)
	; 	ind = where(temp.obs_color gt 0)
	; 	result = poly_fit(temp[ind].obs_color,temp[ind].res_color,2)
	; 	rela_coef_right1.VI = double(result[0:2])
	; endelse
	; endfor
	array = {ID_Huang:0ll,UV_bulge:-999.d,UV_disk:-999.d,VI_bulge:0.d,VI_disk:0.d,z_used:-999.d}
	array = replicate(array,size(galaxy_file,/n_elements))
	id_col = make_array(size(galaxy_file,/n_elements),/L64,VALUE=0)
	uv_bulge_res = make_array(size(galaxy_file,/n_elements),/double,value=-999.)
	vi_bulge_res = make_array(size(galaxy_file,/n_elements),/double,value=-999.)
	uv_disk_res = make_array(size(galaxy_file,/n_elements),/double,value=-999.)
	vi_disk_res = make_array(size(galaxy_file,/n_elements),/double,value=-999.)
	z_used = make_array(size(galaxy_file,/n_elements),/double,value=-999.)
	for num2 = 0 , size(galaxy_file,/n_elements)-1 do begin
		idx = galaxy_file[num2].id_huang
		z = galaxy_file[num2].z_used
		uv_bulge_obs = galaxy_file[num2].uv_bulge_CANDELS_obs
		vi_bulge_obs = galaxy_file[num2].vi_bulge_CANDELS_obs
		uv_disk_obs = galaxy_file[num2].uv_disk_CANDELS_obs
		vi_disk_obs = galaxy_file[num2].vi_disk_CANDELS_obs
		if z gt 0.8 and z lt 1.3 and FINITE(uv_bulge_obs,/NAN) ne 1 and FINITE(vi_bulge_obs,/NAN) ne 1 and FINITE(uv_disk_obs,/NAN) ne 1 and FINITE(vi_disk_obs,/NAN) ne 1 then begin
			rela_coef_left1 = {UV:-999.d,VI:-999.d};dblarr(2,2)
			rela_coef_right1 = {UV:-999.d,VI:-999.d};dblarr(3,2)
			rela_coef_left1 = replicate(rela_coef_left1,2)
			rela_coef_right1 = replicate(rela_coef_right1,3)
			; uv_res = dblarr(sed_num)
			; vi_res = dblarr(sed_num)
			; uv_obs = dblarr(sed_num)
			; vi_obs = dblarr(sed_num)
			for num4 = 0 , relation_counts - 1 do begin
				name = (strsplit(relation[num4],'/',/extract))[n_elements(strsplit(relation[num4],'/',/extract))-1]
				if (strsplit(name,'_',/extract))[2] eq 'UV' and (strsplit(name,'_',/extract))[1] eq strtrim(idx,1) then begin
					temp = mrdfits(relation[num4],1,/silent)
					ind = where(temp.obs_color gt 0)
					if z lt 1 then begin
						result = linfit(temp[ind].obs_color,temp[ind].res_color)
						rela_coef_left1.UV = double(result)
					endif
					if z gt 1 then begin
						result = poly_fit(temp[ind].obs_color,temp[ind].res_color,2)
						rela_coef_right1.UV = double(result[0:2])
					endif
				endif

				if (strsplit(name,'_',/extract))[2] eq 'VI' and (strsplit(name,'_',/extract))[1] eq strtrim(idx,1) then begin
					temp = mrdfits(relation[num4],1,/silent)
					ind = where(temp.obs_color gt 0)
					if z lt 1 then begin
						result = linfit(temp[ind].obs_color,temp[ind].res_color)
						rela_coef_left1.VI = double(result)
					endif
					if z gt 1 then begin
						result = poly_fit(temp[ind].obs_color,temp[ind].res_color,2)
						rela_coef_right1.VI = double(result[0:2])
					endif
				endif
			endfor
			
				if z lt 1 then begin
					if total(rela_coef_left1.UV eq float(-999.0)) eq 0 and total(rela_coef_left1.VI eq float(-999.0)) eq 0 then begin
						id_col[num2] = idx
						uv_bulge_res[num2] = (rela_coef_left1.UV)[0] + (rela_coef_left1.UV)[1] * uv_bulge_obs
						vi_bulge_res[num2] = (rela_coef_left1.VI)[0] + (rela_coef_left1.VI)[1] * vi_bulge_obs
						uv_disk_res[num2] = (rela_coef_left1.UV)[0] + (rela_coef_left1.UV)[1] * uv_disk_obs
						vi_disk_res[num2] = (rela_coef_left1.VI)[0] + (rela_coef_left1.VI)[1] * vi_disk_obs
						z_used[num2] = z
					endif
				endif else begin
					if total(rela_coef_right1.UV eq float(-999.0)) eq 0 and total(rela_coef_right1.VI eq float(-999.0)) eq 0 then begin
						id_col[num2] = idx
						uv_bulge_res[num2] = (rela_coef_right1.UV)[0] + (rela_coef_right1.UV)[1] * uv_bulge_obs + (rela_coef_right1.UV)[2] * uv_bulge_obs^2
						vi_bulge_res[num2] = (rela_coef_right1.VI)[0] + (rela_coef_right1.VI)[1] * vi_bulge_obs + (rela_coef_right1.VI)[2] * vi_bulge_obs^2
						uv_disk_res[num2] = (rela_coef_right1.UV)[0] + (rela_coef_right1.UV)[1] * uv_disk_obs + (rela_coef_right1.UV)[2] * uv_disk_obs^2
						vi_disk_res[num2] = (rela_coef_right1.VI)[0] + (rela_coef_right1.VI)[1] * vi_disk_obs + (rela_coef_right1.VI)[2] * vi_bulge_obs^2
						z_used[num2] = z
					endif
				endelse
			
		endif
	endfor
	array.id_huang = id_col
	array.UV_bulge = uv_bulge_res
	array.UV_disk = uv_disk_res
	array.VI_bulge = vi_bulge_res
	array.VI_disk = vi_disk_res
	array.z_used = z_used
	mwrfits,array,export_path+fields_list[field]+'_uvi_res_CANDELS_c1_onlygdnfrac.fits';
end
end
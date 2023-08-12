hc = 6.62607*1e-27*3*1e10
sed_path = '/Users/lpr/Data/fits/pridata/BrownSEDtemplates/Brown_sed'
rest_filter_path = '/Users/lpr/Data/fits/pridata/filter/rest_frame/'
obs_filter_path = '/Users/lpr/Data/fits/pridata/filter/obs_frame/'
export_path = '/Users/lpr/Data/fits/expdata/HST/goodsn_all/color_difference/'
sed_file = file_search(sed_path,'*.dat',count=sed_num)
res_filter_file = file_search(rest_filter_path,'*.dat',count=rest_filter_num)
obs_filter_file = file_search(obs_filter_path,'*.dat',count=obs_filter_num)
galaxy_file = mrdfits('/Users/lpr/Data/fits/expdata/HST/goodsn_all/uvi/uvi_obs_CANDELS_useauto.fits',1)
aper_color_cata = mrdfits('/Users/lpr/Desktop/goodsn_CANDELScolor.fits',1)
galaxy_info = mrdfits('/Users/lpr/Data/fits/expdata/HST/goodsn_all/goodsn_Huangall_radecmatch_modifyz.fits',1)
readcol,rest_filter_path+'Generic_Bessell.U.dat',bessell_lambda_u,transmittance_u_bessell,format='f'
readcol,rest_filter_path+'Generic_Bessell.V.dat',bessell_lambda_v,transmittance_v_bessell,format='f'
readcol,rest_filter_path+'Generic_Bessell.I.dat',bessell_lambda_i,transmittance_i_bessell,format='f'
array = {ID_Huang:0ll,UV_bulge:-999.d,UV_disk:-999.d,VI_bulge:0.d,VI_disk:0.d}
; array = replicate(array,size(galaxy_file,/n_elements))
; id_col = MAKE_ARRAY(size(galaxy_file,/n_elements),/L64,VALUE=0)
; uv_bulge_res = dblarr(size(galaxy_file,/n_elements))
; vi_bulge_res = dblarr(size(galaxy_file,/n_elements))
; uv_disk_res = dblarr(size(galaxy_file,/n_elements))
; vi_disk_res = dblarr(size(galaxy_file,/n_elements))
array = replicate(array,size(aper_color_cata,/n_elements))
id_col = MAKE_ARRAY(size(aper_color_cata,/n_elements),/L64,VALUE=0)
uv_bulge_res = dblarr(size(aper_color_cata,/n_elements))
vi_bulge_res = dblarr(size(aper_color_cata,/n_elements))
uv_disk_res = dblarr(size(aper_color_cata,/n_elements))
vi_disk_res = dblarr(size(aper_color_cata,/n_elements))
for num2 = 0 , size(aper_color_cata,/n_elements)-1 do begin; for num2 = 0 , size(galaxy_file,/n_elements)-1 do begin
	idx = aper_color_cata[num2].id_huang; 	idx = galaxy_file[num2].id_huang
	if where(galaxy_info.id eq idx) ne -1 then begin
		z = galaxy_info[where(galaxy_info.id eq idx)].z_used
		; u_hst_filter = galaxy_file[num2].u_hst_filter
		; v_hst_filter = galaxy_file[num2].v_hst_filter
		; i_hst_filter = galaxy_file[num2].i_hst_filter
		u_hst_filter = galaxy_file[where(galaxy_file.id_huang eq idx)].u_hst_filter
		v_hst_filter = galaxy_file[where(galaxy_file.id_huang eq idx)].v_hst_filter
		i_hst_filter = galaxy_file[where(galaxy_file.id_huang eq idx)].i_hst_filter
		uv_bulge_obs = aper_color_cata[where(aper_color_cata.id_huang eq idx)].uv_bulge_CANDELS_obs
		vi_bulge_obs = aper_color_cata[where(aper_color_cata.id_huang eq idx)].vi_bulge_CANDELS_obs
		uv_disk_obs = aper_color_cata[where(aper_color_cata.id_huang eq idx)].uv_disk_CANDELS_obs
		vi_disk_obs = aper_color_cata[where(aper_color_cata.id_huang eq idx)].vi_disk_CANDELS_obs
		; uv_bulge_obs = galaxy_file[num2].uv_bulge
		; vi_bulge_obs = galaxy_file[num2].vi_bulge
		; uv_disk_obs = galaxy_file[num2].uv_disk
		; vi_disk_obs = galaxy_file[num2].vi_disk
		; here to open hst filter, which this galaxy used
		for num3 = 0 , obs_filter_num-1 do begin
			if strpos(obs_filter_file[num3],strtrim(u_hst_filter,1)) ne -1 then begin
				readcol,obs_filter_file[num3],trans_lambda_u_obs,transmittance_u_obs,format='f'
			endif else begin
				if strpos(obs_filter_file[num3],strtrim(v_hst_filter,1)) ne -1 then begin
					readcol,obs_filter_file[num3],trans_lambda_v_obs,transmittance_v_obs,format='f'
				endif else begin
					if strpos(obs_filter_file[num3],strtrim(i_hst_filter,1)) ne -1 then begin
						readcol,obs_filter_file[num3],trans_lambda_i_obs,transmittance_i_obs,format='f'
					endif
				endelse
			endelse
		endfor
		; here comes the loop of every sed template
		uv_res = dblarr(sed_num)
		vi_res = dblarr(sed_num)
		uv_obs = dblarr(sed_num)
		vi_obs = dblarr(sed_num)
		for num1 = 0 , sed_num-1 do begin
			readcol,sed_file[num1],sed_lambda,sed_flux,format='f'
			interpolated_u = interpol(sed_flux,sed_lambda,bessell_lambda_u)
			interpolated_v = interpol(sed_flux,sed_lambda,bessell_lambda_v)
			interpolated_i = interpol(sed_flux,sed_lambda,bessell_lambda_i)
			convolve_u = int_tabulated(bessell_lambda_u,interpolated_u*(hc*bessell_lambda_u)^(-1)*transmittance_u_bessell)
			convolve_v = int_tabulated(bessell_lambda_v,interpolated_v*(hc*bessell_lambda_v)^(-1)*transmittance_v_bessell)
			convolve_i = int_tabulated(bessell_lambda_i,interpolated_i*(hc*bessell_lambda_i)^(-1)*transmittance_i_bessell)
			convolve_denominator_u = int_tabulated(bessell_lambda_u,transmittance_u_bessell*3631*1e-23*(hc*bessell_lambda_u)^(-1)/(bessell_lambda_u^2))
			convolve_denominator_v = int_tabulated(bessell_lambda_v,transmittance_v_bessell*3631*1e-23*(hc*bessell_lambda_v)^(-1)/(bessell_lambda_v^2))
			convolve_denominator_i = int_tabulated(bessell_lambda_i,transmittance_u_bessell*3631*1e-23*(hc*bessell_lambda_i)^(-1)/(bessell_lambda_i^2))
			mag_u = -2.5*alog10(convolve_u/convolve_denominator_u)
			mag_v = -2.5*alog10(convolve_v/convolve_denominator_v)
			mag_i = -2.5*alog10(convolve_i/convolve_denominator_i)
			uv_res[num1] = mag_u - mag_v
			vi_res[num1] = mag_v - mag_i
			interpolated_u_obs = interpol(sed_flux,sed_lambda*(1+z),trans_lambda_u_obs)
			interpolated_v_obs = interpol(sed_flux,sed_lambda*(1+z),trans_lambda_v_obs)
			interpolated_i_obs = interpol(sed_flux,sed_lambda*(1+z),trans_lambda_i_obs)
			convolve_u_obs = int_tabulated(trans_lambda_u_obs,interpolated_u_obs*(hc*trans_lambda_u_obs)^(-1)*transmittance_u_obs)
			convolve_v_obs = int_tabulated(trans_lambda_v_obs,interpolated_v_obs*(hc*trans_lambda_v_obs)^(-1)*transmittance_v_obs)
			convolve_i_obs = int_tabulated(trans_lambda_i_obs,interpolated_i_obs*(hc*trans_lambda_i_obs)^(-1)*transmittance_i_obs)
			convolve_denominator_u_obs = int_tabulated(trans_lambda_u_obs,transmittance_u_obs*3631*1e-23*(hc*trans_lambda_u_obs)^(-1)/(trans_lambda_u_obs^2))
			convolve_denominator_v_obs = int_tabulated(trans_lambda_v_obs,transmittance_v_obs*3631*1e-23*(hc*trans_lambda_v_obs)^(-1)/(trans_lambda_v_obs^2))
			convolve_denominator_i_obs = int_tabulated(trans_lambda_i_obs,transmittance_i_obs*3631*1e-23*(hc*trans_lambda_i_obs)^(-1)/(trans_lambda_i_obs^2))
			mag_u_obs = -2.5*alog10(convolve_u_obs/convolve_denominator_u_obs)
			mag_v_obs = -2.5*alog10(convolve_v_obs/convolve_denominator_v_obs)
			mag_i_obs = -2.5*alog10(convolve_i_obs/convolve_denominator_i_obs)
			uv_obs[num1] = mag_u_obs - mag_v_obs
			vi_obs[num1] = mag_v_obs - mag_i_obs
		endfor
		id_col[num2] = idx
		uv_bulge_res[num2] = interpol(uv_res,uv_obs,uv_bulge_obs)
		vi_bulge_res[num2] = interpol(vi_res,vi_obs,vi_bulge_obs)
		uv_disk_res[num2] = interpol(uv_res,uv_obs,uv_disk_obs)
		vi_disk_res[num2] = interpol(vi_res,vi_obs,vi_disk_obs)
	endif
endfor
array.id_huang = id_col
array.UV_bulge = uv_bulge_res
array.UV_disk = uv_disk_res
array.VI_bulge = vi_bulge_res
array.VI_disk = vi_disk_res
mwrfits,array,export_path+'uvi_color_res_CANDELS_aper.fits';
end
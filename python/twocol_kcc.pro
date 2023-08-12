; Here to calculate the rest-frame color versus observe-frame color,
; then save the results to fits table and images.
; Due to the redshift range and I don't want to do colcol_kcc for every redshift.
; So I divide our sample to two redshift bins: [0.8-1], [1-1.3].
; In [0.8-1], U - V used F606W - F125W, in [1-1.3], U - V used F606W - F125W;
; In [0.8-1], V _ I used F850LP - F160W, in [1-1.3], V _ I used F125W - F160W;

z_range_list = {left1:[0.8,1],right1:[1,1.3]}
; All bands use brown templates.
templates = '/Users/lpr/Data/Brown/Brown_sed/'
; ---------------- low_z ----------------
; bins = 10
; ; U - V
; res_short = '/Users/lpr/Data/fits/pridata/filter/rest_frame/Generic_Bessell.U.dat'
; res_long = '/Users/lpr/Data/fits/pridata/filter/rest_frame/Generic_Bessell.V.dat'
; obs_short = '/Users/lpr/Data/fits/pridata/filter/obs_frame/HST_ACS_WFC.F606W.dat'
; obs_long = '/Users/lpr/Data/fits/pridata/filter/obs_frame/HST_WFC3_IR.F125W.dat'
; output = ['/Users/lpr/Data/fits/expdata/twocolor_kcc/f_nu/UV_f606f125.eps','/Users/lpr/Data/fits/expdata/twocolor_kcc/f_nu/UV_f606f125.fits','F606w-F125w','U-V'] ;_left1_left1
; colcol_kcc,z_range=z_range_list.low_z,res_short=res_short,res_long=res_long, $ 
; 	obs_short=obs_short,obs_long=obs_long,templates=templates,bins=bins,output=output ;,xrange=[-1,4]
; print,'--------- U-V low-z done --------'
; ; ; V - I
; res_short = '/Users/lpr/Data/fits/pridata/filter/rest_frame/Generic_Bessell.V.dat'
; res_long = '/Users/lpr/Data/fits/pridata/filter/rest_frame/Generic_Bessell.I.dat'
; obs_short = '/Users/lpr/Data/fits/pridata/filter/obs_frame/HST_ACS_WFC.F850LP.dat'
; obs_long = '/Users/lpr/Data/fits/pridata/filter/obs_frame/HST_WFC3_IR.F160W.dat'
; output = ['/Users/lpr/Data/fits/expdata/twocolor_kcc/f_nu/VI_f850f160.eps','/Users/lpr/Data/fits/expdata/twocolor_kcc/f_nu/VI_f850f160.fits','F850l-F160w','V-I'];_left1_left1
; colcol_kcc,z_range=z_range_list.low_z,res_short=res_short,res_long=res_long, $ 
; 	obs_short=obs_short,obs_long=obs_long,templates=templates,bins=bins,output=output ;,xrange=[-1,2]
; res_short = '/Users/lpr/Data/fits/pridata/filter/rest_frame/Generic_Bessell.V.dat'
; res_long = '/Users/lpr/Data/fits/pridata/filter/rest_frame/Generic_Bessell.I.dat'
; obs_short = '/Users/lpr/Data/fits/pridata/filter/obs_frame/HST_WFC3_IR.F125W.dat'
; obs_long = '/Users/lpr/Data/fits/pridata/filter/obs_frame/HST_WFC3_IR.F160W.dat'
; output = ['/Users/lpr/Data/fits/expdata/twocolor_kcc/f_nu/VI_f125f160_left1.eps','/Users/lpr/Data/fits/expdata/twocolor_kcc/f_nu/VI_f125f160_left1.fits','F125w-F160w','V-I']
; colcol_kcc,z_range=z_range_list.low_z,res_short=res_short,res_long=res_long, $ 
; 	obs_short=obs_short,obs_long=obs_long,templates=templates,bins=bins,output=output ;,xrange=[-1,1]
; print,'--------- V-I low-z done --------'
; ---------------- high_z ----------------
; bins = 6
; ; U - V
; res_short = '/Users/lpr/Data/fits/pridata/filter/rest_frame/Generic_Bessell.U.dat'
; res_long = '/Users/lpr/Data/fits/pridata/filter/rest_frame/Generic_Bessell.V.dat'
; obs_short = '/Users/lpr/Data/fits/pridata/filter/obs_frame/HST_ACS_WFC.F606W.dat'
; obs_long = '/Users/lpr/Data/fits/pridata/filter/obs_frame/HST_WFC3_IR.F105W.dat'
; output = ['/Users/lpr/Data/fits/expdata/twocolor_kcc/f_nu/UV_f606f105_right1.eps','/Users/lpr/Data/fits/expdata/twocolor_kcc/f_nu/UV_f606f105_right1.fits','F606w-F105w','U-V']
; colcol_kcc,z_range=z_range_list.high_z,res_short=res_short,res_long=res_long, $ 
; 	obs_short=obs_short,obs_long=obs_long,templates=templates,bins=bins,output=output ;,xrange=[-1,4.5]
; print,'--------- U-V high-z done --------'
; ; V - I
; res_short = '/Users/lpr/Data/fits/pridata/filter/rest_frame/Generic_Bessell.V.dat'
; res_long = '/Users/lpr/Data/fits/pridata/filter/rest_frame/Generic_Bessell.I.dat'
; obs_short = '/Users/lpr/Data/fits/pridata/filter/obs_frame/HST_WFC3_IR.F125W.dat'
; obs_long = '/Users/lpr/Data/fits/pridata/filter/obs_frame/HST_WFC3_IR.F160W.dat'
; output = ['/Users/lpr/Data/fits/expdata/twocolor_kcc/f_nu/VI_f125f160_right1.eps','/Users/lpr/Data/fits/expdata/twocolor_kcc/f_nu/VI_f125f160_right1.fits','F125w-F160w','V-I']
; colcol_kcc,z_range=z_range_list.high_z,res_short=res_short,res_long=res_long, $
; 	obs_short=obs_short,obs_long=obs_long,templates=templates,bins=bins,output=output ;xrange=[-0.4,2],
; print,'--------- V-I high-z done --------'
; =====================================================================================
; Fit one relation for every galaxy, which is every redshift.
; =====================================================================================
fields_list = ['goodsn','goodss','egs']
catalog_path = '/Users/lpr/Data/lirg_project/output/catalog/'
catalog_fix = '_Huangall_candels_radec_van_modifyz.fits'
export_path = '/Users/lpr/Data/lirg_project/output/twocolor_kcc/f_nu/'
filter_path = '/Users/lpr/Data/filter/'
uv_filters_res = {left1:['Generic_Bessell.U.dat','Generic_Bessell.V.dat'],right1:['Generic_Bessell.U.dat','Generic_Bessell.V.dat']}
vi_filters_res = {left1:['Generic_Bessell.V.dat','Generic_Bessell.I.dat'],right1:['Generic_Bessell.V.dat','Generic_Bessell.I.dat']}
uv_filters_obs = {left1:['HST_ACS_WFC.F606W.dat','HST_WFC3_IR.F125W.dat'],right1:['HST_ACS_WFC.F606W.dat','HST_WFC3_IR.F125W.dat']}
vi_filters_obs = {left1:['HST_ACS_WFC.F814W.dat','HST_WFC3_IR.F160W.dat'],right1:['HST_WFC3_IR.F125W.dat','HST_WFC3_IR.F160W.dat']}
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
for fields = 2,2 do begin; for fields = 0 , n_elements(fields_list) - 1 do begin
	; file_mkdir,export_path+fields_list[fields]
	catalog = mrdfits(catalog_path+fields_list[fields]+catalog_fix,1)
	for num1 = 0 , n_elements(catalog) - 1 do begin
		galaid = catalog[num1].id
		z = catalog[num1].z_used
		temp = file_search(export_path+fields_list[fields],'*.fits')
		if z gt 1. and where(temp eq export_path+fields_list[fields]+'/'+fields_list[fields]+'_'+strtrim(galaid,1)+'_UV_'+strmid(uv_filters_obs.right1[0],strpos(uv_filters_obs.right1[0],'.dat')-5,5)+strmid(uv_filters_obs.right1[1],strpos(uv_filters_obs.right1[1],'.dat')-5,5)+'.fits') ne -1 then begin
			file_delete,export_path+fields_list[fields]+'/'+fields_list[fields]+'_'+strtrim(galaid,1)+'_UV_'+strmid(uv_filters_obs.right1[0],strpos(uv_filters_obs.right1[0],'.dat')-5,5)+strmid(uv_filters_obs.right1[1],strpos(uv_filters_obs.right1[1],'.dat')-5,5)+'.fits'
			; First do UV k-correction.
			if z gt z_range_list.left1[0] and z lt z_range_list.left1[1] then begin
				print,'================================= ENTER LEFT REGION ================================='
				; ========================================== U - V color ==========================================
				xtitle = strmid(uv_filters_obs.left1[0],strpos(uv_filters_obs.left1[0],'.dat')-5,5)+' - '+strmid(uv_filters_obs.left1[1],strpos(uv_filters_obs.left1[1],'.dat')-5,5)
				ytitle = 'U - V'
				name = export_path+fields_list[fields]+'/'+fields_list[fields]+'_'+strtrim(galaid,1)+'_UV_'+strmid(uv_filters_obs.left1[0],strpos(uv_filters_obs.left1[0],'.dat')-5,5)+strmid(uv_filters_obs.left1[1],strpos(uv_filters_obs.left1[1],'.dat')-5,5)
				colcol_kcc,z_range=z,res_short=filter_path+'rest_frame/'+uv_filters_res.left1[0],res_long=filter_path+'rest_frame/'+uv_filters_res.left1[1],obs_short=filter_path+'obs_frame/'+uv_filters_obs.left1[0],obs_long=filter_path+'obs_frame/'+uv_filters_obs.left1[1],templates=templates,bins=1,output=[name+'.eps',name+'.fits',xtitle,ytitle]
				; ========================================== V - I color ==========================================
				xtitle = strmid(vi_filters_obs.left1[0],strpos(vi_filters_obs.left1[0],'.dat')-5,5)+' - '+strmid(vi_filters_obs.left1[1],strpos(vi_filters_obs.left1[1],'.dat')-5,5)
				ytitle = 'V - I'
				name = export_path+fields_list[fields]+'/'+fields_list[fields]+'_'+strtrim(galaid,1)+'_VI_'+strmid(vi_filters_obs.left1[0],strpos(vi_filters_obs.left1[0],'.dat')-5,5)+strmid(vi_filters_obs.left1[1],strpos(vi_filters_obs.left1[1],'.dat')-5,5)
				colcol_kcc,z_range=z,res_short=filter_path+'rest_frame/'+vi_filters_res.left1[0],res_long=filter_path+'rest_frame/'+vi_filters_res.left1[1],obs_short=filter_path+'obs_frame/'+vi_filters_obs.left1[0],obs_long=filter_path+'obs_frame/'+vi_filters_obs.left1[1],templates=templates,bins=1,output=[name+'.eps',name+'.fits',xtitle,ytitle]
			endif
			if z gt z_range_list.right1[0] and z lt z_range_list.right1[1] then begin
				print,'================================= ENTER RIGHT REGION ================================='
				; ========================================== U - V color ==========================================
				xtitle = strmid(uv_filters_obs.right1[0],strpos(uv_filters_obs.right1[0],'.dat')-5,5)+' - '+strmid(uv_filters_obs.right1[1],strpos(uv_filters_obs.right1[1],'.dat')-5,5)
				ytitle = 'U - V'
				name = export_path+fields_list[fields]+'/'+fields_list[fields]+'_'+strtrim(galaid,1)+'_UV_'+strmid(uv_filters_obs.right1[0],strpos(uv_filters_obs.right1[0],'.dat')-5,5)+strmid(uv_filters_obs.right1[1],strpos(uv_filters_obs.right1[1],'.dat')-5,5)
				colcol_kcc,z_range=z,res_short=filter_path+'rest_frame/'+uv_filters_res.right1[0],res_long=filter_path+'rest_frame/'+uv_filters_res.right1[1],obs_short=filter_path+'obs_frame/'+uv_filters_obs.right1[0],obs_long=filter_path+'obs_frame/'+uv_filters_obs.right1[1],templates=templates,bins=1,output=[name+'.eps',name+'.fits',xtitle,ytitle]
				; ; ========================================== V - I color ==========================================
				; xtitle = strmid(vi_filters_obs.right1[0],strpos(vi_filters_obs.right1[0],'.dat')-5,5)+' - '+strmid(vi_filters_obs.right1[1],strpos(vi_filters_obs.right1[1],'.dat')-5,5)
				; ytitle = 'V - I'
				; name = export_path+fields_list[fields]+'/'+fields_list[fields]+'_'+strtrim(galaid,1)+'_VI_'+strmid(vi_filters_obs.right1[0],strpos(vi_filters_obs.right1[0],'.dat')-5,5)+strmid(vi_filters_obs.right1[1],strpos(vi_filters_obs.right1[1],'.dat')-5,5)
				; colcol_kcc,z_range=z,res_short=filter_path+'rest_frame/'+vi_filters_res.right1[0],res_long=filter_path+'rest_frame/'+vi_filters_res.right1[1],obs_short=filter_path+'obs_frame/'+vi_filters_obs.right1[0],obs_long=filter_path+'obs_frame/'+vi_filters_obs.right1[1],templates=templates,bins=1,output=[name+'.eps',name+'.fits',xtitle,ytitle]
			endif
		endif
	endfor
	print,'======================================================= '+fields_list[fields]+' done ======================================================='
endfor
end
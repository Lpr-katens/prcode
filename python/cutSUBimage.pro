; ; --------------------------------------------------------------------------
; ; This code is used to extract subimages from CANDELS whole field images.
; ; Because IDL will read a fits image into an array by function readfits.
; ; It will treat float as an inter index to subtract image information. And,
; ; function AD2XY or ADXY do not give round result.
; ; For example: 2.3 -> 2, 2.8 -> 2. So, when I creat sub-images with hextract and
; ; ad2xy, I give round(ad2xy_result) to hextract to ensure it will cut the
; ; correct sub-images with right galaxy centers.
; ; --------------------------------------------------------------------------
fields_list = ['goodsn','goodss','egs'];
images_original_path = '/Users/lpr/Data/lirg_project/intake/CANDELS/'
catalog_original_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
export_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
; image_halfsize_inpixels = 50 ; the whole image with size 101x101
image_halfsize_inpixels = fix(5/0.06)
; a=[[1,1,1],[1,1,1],[1,1,1]]
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
for num1 = 1 , n_elements(fields_list) do begin
	file_mkdir,export_path + fields_list[num1-1] + '_10'
	catalog = mrdfits(catalog_original_path + fields_list[num1-1] + '_Huangall_radec_candels.fits',1,/silent)
	catalog = catalog[where(catalog.separation lt 1)]
	drz_ima_list = file_search(images_original_path + fields_list[num1-1],'*_drz.fits',count = drz_ima_counts)
	wht_ima_list = file_search(images_original_path + fields_list[num1-1],'*wht.fits',count = wht_ima_counts)
	if fields_list[num1-1] ne 'egs' then begin
		for num2 = 1 , drz_ima_counts do begin
			band = strmid(drz_ima_list[num2-1],strpos(drz_ima_list[num2-1],'f',70),5)
			file_mkdir,export_path + fields_list[num1-1] + '_10/' + fields_list[num1-1] + '_' + band
			whole_field_ima = readfits(drz_ima_list[num2-1],hdr,/silent)
			for num3 = 1 , size(catalog,/n_elements) do begin
				galaid = catalog[num3-1].ID
				; if galaid eq 32411 and band eq 'f160w' then begin
					; print,band
				subimages_list = file_search(export_path+fields_list[num1-1]+'/'+fields_list[num1-1]+'_'+band,'*.fits',count=check_counts)
				check = 1
				if band ne 'f160w' then check = 0
				; for number_check = 1 , check_counts do begin
				; 	if (strsplit(subimages_list[number_check-1],'/',/extract))[n_elements(strsplit(subimages_list[number_check-1],'/',/extract))-1] eq fields_list[num1-1]+'_'+band+'_'+strtrim(galaid,1)+'.fits' then check = 0
				; endfor
				if check eq 1 then begin
					ra = catalog[num3-1].ra_candels
					dec = catalog[num3-1].dec_candels
					extast,hdr,astr
					AD2XY,ra,dec,astr,x,y
					if round(x) lt (size(whole_field_ima,/dimensions))[0]-image_halfsize_inpixels and round(x) gt image_halfsize_inpixels and round(y) lt (size(whole_field_ima,/dimensions))[1]-image_halfsize_inpixels and round(y) gt image_halfsize_inpixels then begin
						print,'full image'
						print,galaid
						hextract,whole_field_ima,hdr,newima,newhdr,round(x)-image_halfsize_inpixels,round(x)+image_halfsize_inpixels,round(y)-image_halfsize_inpixels,round(y)+image_halfsize_inpixels,/silent
						; writefits,'/Users/lpr/Data/lirg_project/output/van_re_half_light_radius/goodsn_f160w_32411.fits',newima,newhdr,/silent
						writefits,export_path+fields_list[num1-1]+'_10/'+fields_list[num1-1]+'_'+band+'/'+fields_list[num1-1]+'_'+band+'_'+strtrim(galaid,1)+'.fits',newima,newhdr,/silent
						print,'================== drizzle '+fields_list[num1-1]+band+'_'+strtrim(galaid,1)+' is done=================='
					endif else begin
						if round(x) lt (size(whole_field_ima,/dimensions))[0] and round(x) gt 0 and round(y) lt (size(whole_field_ima,/dimensions))[1] and round(y) gt 0 then begin
							print,'!!! not full image !!!'
							hextract,whole_field_ima,hdr,newima,newhdr,max([round(x)-image_halfsize_inpixels,0]),min([round(x)+image_halfsize_inpixels,(size(whole_field_ima,/dimensions))[0]])-1,max([round(y)-image_halfsize_inpixels,0]),min([round(y)+image_halfsize_inpixels,(size(whole_field_ima,/dimensions))[1]])-1,/silent
							; writefits,'/Users/lpr/Data/lirg_project/output/van_re_half_light_radius/goodsn_f160w_32411.fits',newima,newhdr,/silent
							writefits,export_path+fields_list[num1-1]+'_10/'+fields_list[num1-1]+'_'+band+'/'+fields_list[num1-1]+'_'+band+'_'+strtrim(galaid,1)+'.fits',newima,newhdr,/silent
							print,'================== drizzle '+fields_list[num1-1]+band+'_'+strtrim(galaid,1)+' is done=================='
						endif
					endelse
				endif
				; endif
			endfor
		endfor
		for num4 = 1 , wht_ima_counts do begin
			band = strmid(wht_ima_list[num4-1],strpos(wht_ima_list[num4-1],'f',70),5)
			whole_field_ima = readfits(wht_ima_list[num4-1],hdr,/silent)
			for num5 = 1 , size(catalog,/n_elements) do begin
				galaid = catalog[num5-1].ID
				subimages_list = file_search(export_path+fields_list[num1-1]+'/'+fields_list[num1-1]+'_'+band,'*_wht.fits',count=check_counts)
				check = 1
				if band ne 'f160w' then check = 0
				; for number_check = 1 , check_counts do begin
				; 	if (strsplit(subimages_list[number_check-1],'/',/extract))[n_elements(strsplit(subimages_list[number_check-1],'/',/extract))-1] eq fields_list[num1-1]+'_'+band+'_'+strtrim(galaid,1)+'_wht.fits' then check = 0
				; endfor
				if check eq 1 then begin
					ra = catalog[num5-1].ra_candels
					dec = catalog[num5-1].dec_candels
					extast,hdr,astr
					AD2XY,ra,dec,astr,x,y
					if round(x) lt (size(whole_field_ima,/dimensions))[0]-image_halfsize_inpixels and round(x) gt image_halfsize_inpixels and round(y) lt (size(whole_field_ima,/dimensions))[1]-image_halfsize_inpixels and round(y) gt image_halfsize_inpixels then begin
						print,'full image'
						hextract,whole_field_ima,hdr,newima,newhdr,round(x)-image_halfsize_inpixels,round(x)+image_halfsize_inpixels,round(y)-image_halfsize_inpixels,round(y)+image_halfsize_inpixels,/silent
						writefits,export_path+fields_list[num1-1]+'_10/'+fields_list[num1-1]+'_'+band+'/'+fields_list[num1-1]+'_'+band+'_'+strtrim(galaid,1)+'_wht.fits',a,hdr,/silent
						print,'================== weight '+fields_list[num1-1]+band+'_'+strtrim(galaid,1)+' is done=================='
					endif else begin
						if round(x) lt (size(whole_field_ima,/dimensions))[0] and round(x) gt 0 and round(y) lt (size(whole_field_ima,/dimensions))[1] and round(y) gt 0 then begin
							print,'!!! not full image !!!'
							hextract,whole_field_ima,hdr,newima,newhdr,max([round(x)-image_halfsize_inpixels,0]),min([round(x)+image_halfsize_inpixels,(size(whole_field_ima,/dimensions))[0]])-1,max([round(y)-image_halfsize_inpixels,0]),min([round(y)+image_halfsize_inpixels,(size(whole_field_ima,/dimensions))[1]])-1,/silent
							writefits,export_path+fields_list[num1-1]+'_10/'+fields_list[num1-1]+'_'+band+'/'+fields_list[num1-1]+'_'+band+'_'+strtrim(galaid,1)+'_wht.fits',newima,newhdr,/silent
							print,'================== weight '+fields_list[num1-1]+band+'_'+strtrim(galaid,1)+' is done=================='
						endif
					endelse
				endif
			endfor
		endfor
	endif else begin
		for num2 = 1 , drz_ima_counts do begin
			band = strmid(drz_ima_list[num2-1],strpos(drz_ima_list[num2-1],'f',62),5)
			; band = 'f160w'
			file_mkdir,export_path + fields_list[num1-1] + '_10/' + fields_list[num1-1] + '_' + band
			whole_field_ima = readfits(drz_ima_list[num2-1],hdr,/silent)
			; whole_field_ima = readfits('/Users/lpr/Data/lirg_project/intake/CANDELS/egs/aegis_3dhst.v4.0.F160W_orig_sci.fits',hdr,/silent)
			for num3 = 1 , size(catalog,/n_elements) do begin
				galaid = catalog[num3-1].ID
				; if galaid eq 13025385 and band eq 'f160w' then begin
				subimages_list = file_search(export_path+fields_list[num1-1]+'/'+fields_list[num1-1]+'_'+band,'*.fits',count=check_counts)
				check = 1
				if band ne 'f160w' then check = 0
				; for number_check = 1 , check_counts do begin
				; 	if (strsplit(subimages_list[number_check-1],'/',/extract))[n_elements(strsplit(subimages_list[number_check-1],'/',/extract))-1] eq fields_list[num1-1]+'_'+band+'_'+strtrim(galaid,1)+'.fits' then check = 0
				; endfor
				if check eq 1 then begin
					ra = catalog[num3-1].ra_candels
					dec = catalog[num3-1].dec_candels
					extast,hdr,astr
					AD2XY,ra,dec,astr,x,y
					print,ra,dec,x,y
					if round(x) lt (size(whole_field_ima,/dimensions))[0]-image_halfsize_inpixels and round(x) gt image_halfsize_inpixels and round(y) lt (size(whole_field_ima,/dimensions))[1]-image_halfsize_inpixels and round(y) gt image_halfsize_inpixels then begin
						print,'full image'
						print,galaid
						hextract,whole_field_ima,hdr,newima,newhdr,round(x)-image_halfsize_inpixels,round(x)+image_halfsize_inpixels,round(y)-image_halfsize_inpixels,round(y)+image_halfsize_inpixels,/silent
						writefits,export_path+fields_list[num1-1]+'_10/'+fields_list[num1-1]+'_'+band+'/'+fields_list[num1-1]+'_'+band+'_'+strtrim(galaid,1)+'.fits',newima,newhdr,/silent
						print,'================== drizzle '+fields_list[num1-1]+band+'_'+strtrim(galaid,1)+' is done=================='
					endif else begin
						if round(x) lt (size(whole_field_ima,/dimensions))[0] and round(x) gt 0 and round(y) lt (size(whole_field_ima,/dimensions))[1] and round(y) gt 0 then begin
							print,'!!! not full image !!!'
							hextract,whole_field_ima,hdr,newima,newhdr,max([round(x)-image_halfsize_inpixels,0]),min([round(x)+image_halfsize_inpixels,(size(whole_field_ima,/dimensions))[0]])-1,max([round(y)-image_halfsize_inpixels,0]),min([round(y)+image_halfsize_inpixels,(size(whole_field_ima,/dimensions))[1]])-1,/silent
							writefits,export_path+fields_list[num1-1]+'_10/'+fields_list[num1-1]+'_'+band+'/'+fields_list[num1-1]+'_'+band+'_'+strtrim(galaid,1)+'.fits',newima,newhdr,/silent
							print,'================== drizzle '+fields_list[num1-1]+band+'_'+strtrim(galaid,1)+' is done=================='
						endif
					endelse
				endif
				; endif
			endfor
		endfor
		for num4 = 1 , wht_ima_counts do begin
			band = strmid(wht_ima_list[num4-1],strpos(wht_ima_list[num4-1],'f',62),5)
			whole_field_ima = readfits(wht_ima_list[num4-1],hdr,/silent)
			for num5 = 1 , size(catalog,/n_elements) do begin
				galaid = catalog[num5-1].ID
				subimages_list = file_search(export_path+fields_list[num1-1]+'/'+fields_list[num1-1]+'_'+band,'*_wht.fits',count=check_counts)
				check = 1
				if band ne 'f160w' then check = 0
				; for number_check = 1 , check_counts do begin
				; 	if (strsplit(subimages_list[number_check-1],'/',/extract))[n_elements(strsplit(subimages_list[number_check-1],'/',/extract))-1] eq fields_list[num1-1]+'_'+band+'_'+strtrim(galaid,1)+'_wht.fits' then check = 0
				; endfor
				if check eq 1 then begin
					ra = catalog[num5-1].ra_candels
					dec = catalog[num5-1].dec_candels
					extast,hdr,astr
					AD2XY,ra,dec,astr,x,y
					if round(x) lt (size(whole_field_ima,/dimensions))[0]-image_halfsize_inpixels and round(x) gt image_halfsize_inpixels and round(y) lt (size(whole_field_ima,/dimensions))[1]-image_halfsize_inpixels and round(y) gt image_halfsize_inpixels then begin
						print,'full image'
						hextract,whole_field_ima,hdr,newima,newhdr,round(x)-image_halfsize_inpixels,round(x)+image_halfsize_inpixels,round(y)-image_halfsize_inpixels,round(y)+image_halfsize_inpixels,/silent
						writefits,export_path+fields_list[num1-1]+'_10/'+fields_list[num1-1]+'_'+band+'/'+fields_list[num1-1]+'_'+band+'_'+strtrim(galaid,1)+'_wht.fits',a,hdr,/silent
						print,'================== weight '+fields_list[num1-1]+band+'_'+strtrim(galaid,1)+' is done=================='
					endif else begin
						if round(x) lt (size(whole_field_ima,/dimensions))[0] and round(x) gt 0 and round(y) lt (size(whole_field_ima,/dimensions))[1] and round(y) gt 0 then begin
							print,'!!! not full image !!!'
							hextract,whole_field_ima,hdr,newima,newhdr,max([round(x)-image_halfsize_inpixels,0]),min([round(x)+image_halfsize_inpixels,(size(whole_field_ima,/dimensions))[0]])-1,max([round(y)-image_halfsize_inpixels,0]),min([round(y)+image_halfsize_inpixels,(size(whole_field_ima,/dimensions))[1]])-1,/silent
							writefits,export_path+fields_list[num1-1]+'_10/'+fields_list[num1-1]+'_'+band+'/'+fields_list[num1-1]+'_'+band+'_'+strtrim(galaid,1)+'_wht.fits',newima,newhdr,/silent
							print,'================== weight '+fields_list[num1-1]+band+'_'+strtrim(galaid,1)+' is done=================='
						endif
					endelse
				endif
			endfor
		endfor
	endelse
endfor
end
; ; extast
; ; ad2xy,ra,dec,ast,x,y
; ; image_halfsize_inpixels = 50
; ; hextract,img,hdr,newimg,newhdr,round(x)-image_halfsize_inpixels,round(y)-image_halfsize_inpixels
; ; python photometry
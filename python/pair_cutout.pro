; --------------------------------------------------------------------------
; 找到了pair，把有pair的星系cut出来，半径是7个角秒。要包含所有潜在可能的pair
; --------------------------------------------------------------------------
fields_list = ['goodsn','goodss','egs'];
images_original_path = '/Users/lpr/Data/lirg_project/intake/CANDELS/'
catalog_original_path = '/Users/lpr/Data/lirg_project/output/catalog_radec/'
export_path = '/Users/lpr/Data/lirg_project/output/pair/'
image_halfsize_inpixels =  6/0.06 ;我们寻找pair的时候是以5.35"的半径进行搜寻的，为了包含进pair的整个星系，取7"的半径进行图片裁剪

for num1 = 1 , n_elements(fields_list) do begin
	file_mkdir,export_path + fields_list[num1-1] + '_1000kms'
	pair_ctg = mrdfits(catalog_original_path + fields_list[num1-1] + '_pair_1000kms.fits',1,/silent)
	infor_ctg = mrdfits(catalog_original_path + fields_list[num1-1] + '_Huangall_candels_van_params_zmasshuangsfr.fits',1,/silent)
	drz_ima_list = file_search(images_original_path + fields_list[num1-1],'*drz.fits',count = drz_ima_counts)
	if fields_list[num1-1] ne 'egs' then begin
		for num2 = 1 , drz_ima_counts do begin
			band = strmid(drz_ima_list[num2-1],strpos(drz_ima_list[num2-1],'f',70),5)
			if band eq 'f160w' then begin
				file_mkdir,export_path + fields_list[num1-1] + '_1000kms/' + fields_list[num1-1] + '_' + band
				whole_field_ima = readfits(drz_ima_list[num2-1],hdr,/silent)
				for num3 = 1 , size(pair_ctg,/n_elements) do begin
					obj_id = pair_ctg[num3-1].ID ;目标源的id
					if pair_ctg[num3-1].count_spec+pair_ctg[num3-1].count_phot gt 0 then begin
						print,obj_id
						subimages_list = file_search(export_path+fields_list[num1-1]+'_1000kms/'+fields_list[num1-1]+'_'+band,'*.fits',count=check_counts)
						check = 1 ;检测是否已经有cutoutimage了
						for number_check = 1 , check_counts do begin
							if (strsplit(subimages_list[number_check-1],'_1000kms/',/extract))[n_elements(strsplit(subimages_list[number_check-1],'/',/extract))-1] eq fields_list[num1-1]+'_'+band+'_'+strtrim(obj_id,1)+'.fits' then check = 0 ;如果已经存在了，那么check设为0
						endfor
						;如果是check=1，也就是cutout不存在，那么就从大图上裁剪小图
						if check eq 1 and pair_ctg[num3-1].count_spec + pair_ctg[num3-1].count_phot gt 0 then begin
							ra = infor_ctg[where(infor_ctg.id eq obj_id)].ra_candels
							dec = infor_ctg[where(infor_ctg.id eq obj_id)].dec_candels
							extast,hdr,astr
							AD2XY,ra,dec,astr,x,y
							if round(x) lt (size(whole_field_ima,/dimensions))[0]-image_halfsize_inpixels and round(x) gt image_halfsize_inpixels and round(y) lt (size(whole_field_ima,/dimensions))[1]-image_halfsize_inpixels and round(y) gt image_halfsize_inpixels then begin
								print,'full image'
								print,obj_id
								hextract,whole_field_ima,hdr,newima,newhdr,round(x)-image_halfsize_inpixels,round(x)+image_halfsize_inpixels,round(y)-image_halfsize_inpixels,round(y)+image_halfsize_inpixels,/silent
								writefits,export_path+fields_list[num1-1]+'_1000kms/'+fields_list[num1-1]+'_'+band+'/'+fields_list[num1-1]+'_'+band+'_'+strtrim(obj_id,1)+'.fits',newima,newhdr,/silent
								print,'================== drizzle '+fields_list[num1-1]+band+'_'+strtrim(obj_id,1)+' is done=================='
							endif else begin
								if round(x) lt (size(whole_field_ima,/dimensions))[0] and round(x) gt 0 and round(y) lt (size(whole_field_ima,/dimensions))[1] and round(y) gt 0 then begin
									print,'!!! not full image !!!'
									hextract,whole_field_ima,hdr,newima,newhdr,max([round(x)-image_halfsize_inpixels,0]),min([round(x)+image_halfsize_inpixels,(size(whole_field_ima,/dimensions))[0]])-1,max([round(y)-image_halfsize_inpixels,0]),min([round(y)+image_halfsize_inpixels,(size(whole_field_ima,/dimensions))[1]])-1,/silent
									writefits,export_path+fields_list[num1-1]+'_1000kms/'+fields_list[num1-1]+'_'+band+'/'+fields_list[num1-1]+'_'+band+'_'+strtrim(obj_id,1)+'.fits',newima,newhdr,/silent
									print,'================== drizzle '+fields_list[num1-1]+band+'_'+strtrim(obj_id,1)+' is done=================='
								endif
							endelse
						endif
					endif
				endfor
			endif
		endfor
	endif else begin
		for num2 = 1 , drz_ima_counts do begin
			; band = strmid(drz_ima_list[num2-1],strpos(drz_ima_list[num2-1],'f',62),5)
			band = 'f160w'
			if band eq 'f160w' then begin
				file_mkdir,export_path + fields_list[num1-1] + '_1000kms/' + fields_list[num1-1] + '_' + band
				whole_field_ima = readfits(drz_ima_list[num2-1],hdr,/silent)
				for num3 = 1 , size(pair_ctg,/n_elements) do begin
					obj_id = pair_ctg[num3-1].ID
					if pair_ctg[num3-1].count_spec+pair_ctg[num3-1].count_phot gt 0 then begin
						print,obj_id
						print,'         *              *'
						print,'        %$@            !@#'
						print,'       &%*&%          ^$#!@'
						print,'      <:{!+@"|&^*&%$#!@%"}#)'
						print,'     #!@%$**************&%$@~'
						print,'    &++&@%***  code  ***%&%*&%'
						print,'    @%*#!&**   begin  **"?)^&%'
						print,'     %*#!&**   now    **?%*@#'
						print,'      ++&@**************@~#^'
						print,'          #!*&&^$#!@%*()'
						print,'            $%(_+(%@?^'
						print,'         )!(#*&++&@%%!*_#+'
						print,'      !@#$)(*&^$#!@%*()*&%$@'
						print,'    (*&^$#!@%*()*&%$@?#^%*(*&*'
						print,'  #$)(*&^$#<{}"?)^&%$@?#^%*(*&*)'
						print,' !@#$)(*&^$#!@%*&%*&%$@~#^%*(*&*)'
						print,'  )(*&^$#!@%*$#*&^$#!@%*()*&%$@~'
						print,'   $)(*&^$#!@%*!@*&%$@|#^%&%$@~'
						print,'    (*&^$#!@%*@|_>"}#)$(*()*&%          **   **'
						print,'     (*&^$#!@%*#!*&&^$#!@%*()       **          **'
						print,'      )*&%$@|%*@%*&&^$#!@%*(      **      **     **'
						print,'       @|#^%&%$>:}?%*@#!@%*     **       **      **'
						print,'        @%*@#!@%*!@%*$#*&@*     **         **   **'
						print,'        !@%*()*&:>%#!~>:"& ** **'
						print,'        <:{!+@"|&^*&%$#!@%   *'
						print,'       *$#!+@"*&  &%$@~<{}"'
						subimages_list = file_search(export_path+fields_list[num1-1]+'_1000kms/'+fields_list[num1-1]+'_'+band,'*.fits',count=check_counts)
						check = 1
						for number_check = 1 , check_counts do begin
							if (strsplit(subimages_list[number_check-1],'_1000kms/',/extract))[n_elements(strsplit(subimages_list[number_check-1],'/',/extract))-1] eq fields_list[num1-1]+'_'+band+'_'+strtrim(obj_id,1)+'.fits' then check = 0
						endfor
						if check eq 1 and pair_ctg[num3-1].count_spec + pair_ctg[num3-1].count_phot gt 0 then begin
							ra = infor_ctg[where(infor_ctg.id eq obj_id)].ra_candels
							dec = infor_ctg[where(infor_ctg.id eq obj_id)].dec_candels
							extast,hdr,astr
							AD2XY,ra,dec,astr,x,y
							if round(x) lt (size(whole_field_ima,/dimensions))[0]-image_halfsize_inpixels and round(x) gt image_halfsize_inpixels and round(y) lt (size(whole_field_ima,/dimensions))[1]-image_halfsize_inpixels and round(y) gt image_halfsize_inpixels then begin
								print,'full image'
								print,obj_id
								hextract,whole_field_ima,hdr,newima,newhdr,round(x)-image_halfsize_inpixels,round(x)+image_halfsize_inpixels,round(y)-image_halfsize_inpixels,round(y)+image_halfsize_inpixels,/silent
								writefits,export_path+fields_list[num1-1]+'_1000kms/'+fields_list[num1-1]+'_'+band+'/'+fields_list[num1-1]+'_'+band+'_'+strtrim(obj_id,1)+'.fits',newima,newhdr,/silent
								print,'================== drizzle '+fields_list[num1-1]+band+'_'+strtrim(obj_id,1)+' is done=================='
							endif else begin
								if round(x) lt (size(whole_field_ima,/dimensions))[0] and round(x) gt 0 and round(y) lt (size(whole_field_ima,/dimensions))[1] and round(y) gt 0 then begin
									print,'!!! not full image !!!'
									hextract,whole_field_ima,hdr,newima,newhdr,max([round(x)-image_halfsize_inpixels,0]),min([round(x)+image_halfsize_inpixels,(size(whole_field_ima,/dimensions))[0]])-1,max([round(y)-image_halfsize_inpixels,0]),min([round(y)+image_halfsize_inpixels,(size(whole_field_ima,/dimensions))[1]])-1,/silent
									writefits,export_path+fields_list[num1-1]+'_1000kms/'+fields_list[num1-1]+'_'+band+'/'+fields_list[num1-1]+'_'+band+'_'+strtrim(obj_id,1)+'.fits',newima,newhdr,/silent
									print,'================== drizzle '+fields_list[num1-1]+band+'_'+strtrim(obj_id,1)+' is done=================='
								endif
							endelse
						endif
					endif
				endfor
			endif
		endfor
	endelse
endfor
end
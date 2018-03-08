 pro firasin,dsname,header,dsfield,product_type,subscr,faceno,badpixval,$
       data,frequency,res,frequnits,title,units,weight,weights,wunits,$
       wtitle,instrume

;
; This is not a stand alone procedure. It is called from DCONVERT/DATAIN.
; This reads in the FIRAS Project Data Set binary extension FITS
; files.

IF (!cgis_os EQ 'vms') THEN eod=']'
IF (!cgis_os EQ 'unix') THEN eod='/'
IF (!cgis_os EQ 'windows') THEN eod='\'

   If strtrim(instrume) eq 'FIRAS' then begin
            product   = strtrim(sxpar(header, 'PRODUCT'))
	    res       = sxpar(header, 'PIXRESOL')
	    nu_zero   = sxpar(header, 'NU_ZERO')
            delta_nu  = sxpar(header, 'DELTA_NU')
	    num_freq  = sxpar(header, 'NUM_FREQ')
	    if num_freq ne 0 then begin
		   frequency = fltarr(num_freq)
		   for i = 0, num_freq-1 do $
		       frequency(i) = nu_zero + (delta_nu * i)
		   frequnits = 'Frequency (GHz)'
	    endif else freq_in = 'N'

	    title = strmid(dsname,0,strpos(dsname,'.'))

	       ; read in the header for the fits BINARY TABLE extension
;	       fxbopen,unit,dsname,1,bin_header
            all_data = readfits(dsname,bin_header,/ext)

	    fields = fxpar(bin_header,'TTYPE*',count=fieldcnt)
	    unitlist  = fxpar(bin_header,'TUNIT*',count=unitcnt)

            if strtrim(dsfield,2) eq '?' then begin
               print,'This file contains the following fields:'
	       for i = 0,fieldcnt-1 do print, fields(i)
               dsfield = ''
               read,dsfield,prompt='Enter Field Name as shown: '
               dsfield = strupcase(dsfield)
               while strlen(dsfield) lt 8 do dsfield = dsfield + ' '
	    endif

; trap for coadded interferograms
            if strpos(dsfield,'COADDED_') ge 0 then begin
                  frequency=indgen(512)
                  frequnits='Position'
            endif

            y_cvec = 0
	    if strpos(dsfield,'+WEIGHT') ge 0 then begin
                   print,'NOTE: these weights are a good estimate of the '+$
                      'relative weights.' 
                   print,'However, for accurate absolute '+$
                      'errors/weights please read the FIRAS Explanatory '+$
                      'Supplement.' 
                   print,'You can then compute the exact weight '+$
                      'and Associate this weight with the Primary Data field.'
		   weight = '3'
		   wtitle = 'weight_of_'+strmid(dsname,0,strpos(dsname,'.'))
		   wname = dsname
                   wunits='sr^2/MJy^2'
;
; read in weights - these are different for different products
                   case 1 of
                       product eq 'CHI-SQ': begin
                            print,'The data product '+ product +' has no weight.'
                            weight = 'N'
                            end
                       product eq 'PIXEL-ZM': begin
                            print,'The data product '+ product +' has no weight.'
                            weight = 'N'
                            end
                       product eq 'COADD-ZM': begin
                            print,'The data product '+ product +' has no weight.'
                            weight = 'N'
                            end
                       product eq 'DIRBE-GK': begin
                            print,'The data product '+ product +' has no weight.'
                            weight = 'N'
                            end
                       product eq 'TEMP-MAP': begin
                            print,'The weight for data product TEMP-MAP is '+$
                                  'computed by 1./(TEMP_SIG)^2'
                            tweights = tbget(bin_header,all_data,'TEMP_SIG')
                            good = where(tweights ne 0)
                            tweights(good) = 1./tweights(good)^2
                            dsfield = 'TEMP    '
                            end
                       product eq 'SKY-IFG' or product eq 'CAL-IFG': begin
                            print,'The weight for data product '+product+$
                                  ' is computed by ADJ_NUM_IFGS'
                            tweights = tbget(bin_header,all_data,'ADJ_NUM_')
                            dsfield = 'COADDED_'
                            end
                       product eq 'CAL-SPEC' or product eq 'CAL-DIFF' or $
                          product eq 'SKY-DIFF' or product eq 'SKY-SPEC': begin
                            print,'The weight for data product '+product+$
                                  ' is computed by ADJ_NUM_IFGS/SIGMAS^2'
                            adj_num = tbget(bin_header,all_data,'ADJ_NUM_') 
                            sigmas = tbget(bin_header,all_data,'SIGMAS')
                            good = where(sigmas ne 0)
                            tweights = adj_num * 0.
                            tweights(good) = adj_num(good)/sigmas(good)^2
                            dsfield = 'REAL_SPE'
                            end
                       product eq 'SKY-DEST' or product eq 'DUST-MAP' or $
                           product eq 'CONT-MAP': begin
                            dsfield = 'SPECTRUM'
                            y_cvec = 1
                            if product eq 'CONT-MAP' then begin
                            print,'The weight for data product '+product+$
                                  ' is computed by NUM_IFGS/C-VECTOR^2'
                               tweights = tbget(bin_header,all_data,'NUM_IFGS')
                            endif else begin
                            print,'The weight for data product '+product+$
                                  ' is computed by WEIGHT/C-VECTOR^2'
                               tweights = tbget(bin_header,all_data,'WEIGHT')
                            endelse 
;
                           j=0
                           while (j ne -1) do begin
                              j=strpos(wname,eod,j)
                              if j ne -1 then begin cnt=j
                                j=j+1 
                              end
                           endwhile
                           if n_elements(cnt) eq 0 then begin
                             j=0
                             alt_eod=':'
                             while (j ne -1) do begin
                               j=strpos(wname,alt_eod,j)
                               if j ne -1 then begin cnt=j
                                 j=j+1 
                               end
                             endwhile
                           endif
                            
                           dir=strmid(wname,0,cnt+1)
                           file_spec=dir+'fi*.fit*'
                           status = lstdsnam( 'FITS', file_spec, flist, $
                               extension, fcount )
                           fcount = n_elements( flist )
                           file_list = strarr( fcount + 1)
                           file_list(0)='Select C-VECTOR for Weights- '+dir 
                           file_list(1:fcount)=flist(0:fcount-1)
mod_menu:
                     	   status = umenu( file_list, title = 0, init=1, xpos=0)
	                   wname = file_list( status )
                     	   wheader = headfits(wname)
                           mod_type=sxpar(wheader,'PRODUCT')
                           if strtrim(mod_type) ne 'C-VEC' and strtrim(mod_type,2) $
                             ne '0' then begin
                             !quiet=0
                             message,/informational,'file chosen is not a FIRAS '+$
                               'C-VECTOR solution'
                             !quiet=old_quiet
                             goto, mod_menu
                           endif
                           w_data = readfits(wname,w_header,/ext)
                           c_vec = tbget(w_header,w_data,'C_VECTOR')
	        	   c_vec = c_vec(0:num_freq-1)
                           w_data = 0
                           end
                       else:
                   endcase
            endif

	    choice = where(dsfield Eq fields)
	    units =  unitlist(choice(0))

;	       IDL arrays start at 0, fits columns start at 1 so we have
;	       to add 1 to the variable choice to point to the correct one
	    readfield = choice(0) + 1
	    fieldpointer = choice(0)

	    print,'reading in data ',systime(0)
;	    fxbread,unit,data_list,readfield
            data_list = tbget(bin_header,all_data,readfield)

		print, 'reading in pixel data ',systime(0)
;	       if (choice(0) ne 0) then fxbread,unit,pixel,'PIXEL'
            if (choice(0) ne 0) then pixel = $
                         tbget(bin_header,all_data,'PIXEL') 

;	    fxbclose,unit
            all_data = 0
            dum = temporary(all_data)

	    ; subset the FIRAS data by frequency if asked
	    if subscr ne '' then begin
	          pos = strpos(subscr,':')

	          if pos gt 0 then begin
	    	     num1 = fix(strmid(subscr,0,pos)) - 1 
	    	     len = strlen(subscr) - (pos+1)
	    	     num2 = fix(strmid(subscr,pos+1,len)) -1
		     data_list = data_list(num1:num2,*)
		     frequency = frequency(num1:num2)
		     if y_cvec then c_vec = c_vec(num1:num2)
		     num_freq = num2 - num1 + 1

	          ;else we have a scalar value
	          endif else begin
	    	     num = fix(subscr) - 1 
		     data_list = data_list(num,*)
		     if y_cvec then c_vec = c_vec(num)
		     frequency = frequency(num)
		     num_freq = 1
	          endelse
	    endif 

	    if choice(0) eq 0 then pixel = data_list

	    print, 'beginning pixelization ',systime(0)
;
	    if faceno ne 7 then begin
		  if weight eq '3' then begin
		     pixface, pixel, data_list, moredata=tweights, $
                              faceno=faceno, res=res
		     if !err eq -1 then goto, abort
		     pix2xy, pixel, res=res, data=data_list, raster=data, $
			     bad_pixval=badpixval, /face
		     pix2xy, pixel, res=res, data=tweights, raster=rweights,$
                              /face
		  endif else begin
		     pixface, pixel, data_list, faceno=faceno, res=res
		     if !err eq -1 then goto, abort
		     pix2xy, pixel, res=res, data=data_list, raster=data, $
			     bad_pixval=badpixval, /face
		  endelse
	    endif else begin
		  pix2xy, pixel, res=res, data=data_list, raster=data, $
		          bad_pixval=badpixval
		  if weight eq '3' then $
		     pix2xy, pixel, res=res, data=tweights, raster=rweights
	    endelse

	    if y_cvec then begin
		  dim = size(data)
	          weights = fltarr(dim(1),dim(2),num_freq)
		  good = where(c_vec gt 0.)
		  if good(0) ne -1 then $
		     for i = 0, dim(1)-1 do $
		       for j = 0,dim(2)-1 do $
			   weights(i,j,good) = rweights(i,j)/c_vec(good)^2
	    endif else if weight eq '3' then weights = rweights

	    if fields(fieldpointer) eq 'SIGNAL  ' or $
                  fields(fieldpointer) eq 'SERROR  ' or $
                  fields(fieldpointer) eq 'SIGMAS  ' or $
                  fields(fieldpointer) eq 'REAL_SPE' or $
                  fields(fieldpointer) eq 'IMAG_SPE' or $
                  fields(fieldpointer) eq 'REAL_VAR' or $
                  fields(fieldpointer) eq 'IMAG_VAR' or $
                  fields(fieldpointer) eq 'REAL_IMA' or $
                  fields(fieldpointer) eq 'SPECTRUM' or $
	          fields(fieldpointer) eq 'IM_SPEC ' then $
	          data = data(*,*,0:num_freq-1)

	    data_list = 0
	    pixel = 0
	    if weight ne 'N' then tweights = 0
	    if weight ne 'N' then rweights = 0
   Endif                       ; end of handling firas data

   If instrume eq 'HLIN-MAP' or instrume eq 'LLIN-MAP' then begin
	    res       = sxpar(header, 'PIXRESOL')

	    frequency = sxpar(header, 'LINFRQ*')
	    frequnits = 'Frequency (GHz)'

	    title = strmid(dsname,0,strpos(dsname,'.'))

	    ; read in the header for the fits BINARY TABLE extension
;           fxbopen,unit,dsname,1,bin_header
            all_data = readfits(dsname,bin_header,/ext)

	    fields = fxpar(bin_header,'TTYPE*',count=fieldcnt)
	    unitlist  = fxpar(bin_header,'TUNIT*',count=unitcnt)

            if strtrim(dsfield,2) eq '?' then begin
                print,'This file contains the following fields:'
	        for i = 0,fieldcnt-1 do print, fields(i)
                dsfield = ''
                read,dsfield,prompt='Enter Field Name as shown: '
                dsfield = strupcase(dsfield)
        	while strlen(dsfield) lt 8 do dsfield = dsfield + ' '
	    endif

	    if dsfield Eq 'LINE_FLU+WEIGHT' then begin
                   print,'NOTE: these weights are a good estimate of the '+$
                      'relative weights.' 
                   print,'However, for accurate absolute '+$
                      'errors/weights please read the FIRAS Explanatory '+$
                      'Supplement.' 
                   print,'You can then compute the exact weight '+$
                      'and Associate this weight with the Primary Data field.'
		   dsfield = 'LINE_FLU'
		   weight = '3'
		   wtitle = 'weight_of_'+strmid(dsname,0,strpos(dsname,'.'))
		   wname = dsname
            endif

	    choice = where(STRUPCASE(dsfield) Eq STRUPCASE(fields))
	    units =  unitlist(choice(0))

;	    IDL arrays start at 0, fits columns start at 1 so we have
;	    to add 1 to the variable choice to point to the correct one
	    readfield = choice(0) + 1
            fieldpointer=choice(0)

		print,'reading in data ',systime(0)
;	       fxbread,unit,data_list,readfield
            data_list = tbget(bin_header,all_data,readfield)

		print, 'reading in pixel data ',systime(0)
;	    if (choice(0) ne 0) then fxbread,unit,pixel,'PIXEL'
            if (choice(0) ne 0) then pixel = $
                         tbget(bin_header,all_data,'PIXEL') 

    	    ; bring in the weighting array
	    if weight eq '3' then begin
	           print, 'reading in the weights array ',systime(0)
;		   fxbread,unit,tweights,'LINE_FL2'
                   tweights = tbget(bin_header,all_data,'LINE_FL2')
                   wunits='sr^2m^4/watt^2'
	    endif

;	    fxbclose,unit
            all_data = 0
            dum = temporary(all_data)

	    ; subset the FIRAS line map data by wavelength if that 
	    ; has been requested
	    if subscr ne '' then begin
	          pos = strpos(subscr,':')

	          if pos gt 0 then begin
		     ; if we're dealing with ranges we took the numbers
		     ; in indexed from 1-8 as idl is 0 based we have to
		     ; subtract 1 from each of the index numbers so that
		     ; we use the proper portions of the array
	    	     num1 = fix(strmid(subscr,0,pos)) -1
	    	     len = strlen(subscr) - (pos+1)
	    	     num2 = fix(strmid(subscr,pos+1,len)) -1
		     data_list = data_list(num1:num2,*)
		     if weight eq '3' then tweights = tweights(num1:num2,*)
                     szf = SIZE(frequency)
                     if num2 gt szf(1) then frequency=findgen(num2+1)
		     frequency = frequency(num1:num2)

	          ;else we have a scalar value
	          endif else begin
	    	     num = fix(subscr) - 1
		     data_list = data_list(num,*)
		     if weight eq '3' then tweights = tweights(num,*)
		     frequency = frequency(num)
	          endelse
	    endif
 
            if weight eq '3' then begin
                   good = where(tweights gt 0)
                   tweights(good)=1./tweights(good)^2
            endif
	    if choice(0) eq 0 then pixel = data_list
 	    print, 'beginning pixelization ',systime(0)
 
	    if faceno ne 7 then begin
		  if weight eq '3' then begin
		     pixface, pixel, data_list, moredata=tweights, $
                              faceno=faceno, res=res
		     if !err eq -1 then goto, abort
		     pix2xy, pixel, res=res, data=data_list, raster=data, $
			     bad_pixval=badpixval, /face
		     pix2xy, pixel, res=res, data=tweights, raster=weights,$
                              /face
		  endif else begin
		     pixface, pixel, data_list, faceno=faceno, res=res
		     if !err eq -1 then goto, abort
		     pix2xy, pixel, res=res, data=data_list, raster=data, $
			     bad_pixval=badpixval, /face
		  endelse
	    endif else begin
		  pix2xy, pixel, res=res, data=data_list, raster=data, $
		          bad_pixval=badpixval
		  if weight eq '3' then $
		     pix2xy, pixel, res=res, data=tweights, raster=weights
	    endelse

	    data_list = 0
	    pixel = 0
	    if weight ne 'N' then tweights = 0
   Endif               ; handling of FIRAS line map data
abort:
end

 pro dirbein,dsname,header,dsfield,product_type,subscr,faceno,badpixval,$
       data,frequency,res,frequnits,title,units,weight,dbsig,dbqual

;
; This is not a stand alone procedure. It is called from DCONVERT/DATAIN.
; This reads in the DIRBE Project Data Set binary extension FITS
; files.
; 1997 Feb 1 Changd FXBOPEN/FXBREAD to READFITS/TBGET for speed. J. Newmark
; 1997 Aug 18 - VERSION 2.0 - J. Newmark - major changes -
;              1) Error if read Pass 2B product
;              2) Handle warm products, FSM, ZSMA

 timediff  = 4.103136e08
 maxsecs   = 1.58112e07
 weight    = 'N'   
 version   = sxpar(header, 'VERSION')
 version   = strupcase(strtrim(version))
 if version eq 'PASS 2B' then begin
     error ='DIRBE Pass 2B data is now superseded and is no longer '+$
            'supported by this software. Please contact NSSDC at '+$
            'http://www.gsfc.nasa.gov/cobe/cobe_home.html to obtain '+$
            'the latest data products.'
     message,error
 endif
 res       = sxpar(header, 'PIXRESOL')
              
 frequency = sxpar(header, 'WAVE*')
 frequency = float(strmid(frequency,0,4))
 frequnits = 'Wavelength (microns)'
              
 title = strmid(dsname,0,strpos(dsname,'.'))
              
; read in the header for the fits BINARY TABLE extension
; fxbopen,unit,dsname,1,bin_header
 all_data = readfits(dsname,bin_header,/ext) 
              
 fields = fxpar(bin_header,'TTYPE*',count=fieldcnt)
 unitlist  = fxpar(bin_header,'TUNIT*',count=unitcnt)

 if strtrim(dsfield,2) eq '?' then begin
      print,'This file contains the following fields:'
      for i = 0,fieldcnt-1 do print, fields(i)
      dsfield = ''
      read,dsfield,prompt='Enter Field Name as shown: '
      while strlen(dsfield) lt 8 do dsfield = dsfield + ' '
 endif
              
 choice = where(STRUPCASE(dsfield) Eq STRUPCASE(fields))
 units =  unitlist(choice(0))
              
; IDL arrays start at 0, fits columns start at 1 so we have
; to add 1 to the variable choice to point to the correct one
 readfield = choice(0) + 1
 fieldpointer=choice(0)
              
 print,'reading in DATA ',systime(0)
; fxbread,unit,data_list,readfield
 data_list = tbget(bin_header,all_data,readfield)
              
 print, 'reading in PIXEL_NO data ',systime(0)
; if (STRUPCASE(dsfield) ne 'PIXEL_NO') then fxbread,unit,pixel,'Pixel_no'
 if (STRUPCASE(dsfield) ne 'PIXEL_NO') then $
      pixel=tbget(bin_header,all_data,'Pixel_no')
              
;             
; IF read in photqual field select either SIGMA or QUAL_FLAG
;
 IF (strupcase(fields(fieldpointer)) EQ 'PHOTQUAL') THEN BEGIN
        if keyword_set(dbsig) then dsel=1
        if keyword_set(dbqual) then dsel=2
        if n_elements(dsel) eq 0 then $
               dsel=one_line_menu(['Select','Absolute SIGMA','QUAL_FLAG'],$
                     init=1)
        if (dsel eq 2) then begin
                 data_list = data_list - data_list/256*256
        endif ELSE BEGIN
;
; stuff to extract proper bytes
                data_list=LONG(data_list)
                sz=size(data_list)
                dbsig=fltarr(sz(1),sz(2))
                blow=where(data_list LT 0,cnt)
                if cnt gt 0 then data_list(blow)=data_list(blow)+65536

                print,'reading in PHOTOMET data for sigma...',systime(0)
   	        choice = where('PHOTOMET' Eq STRUPCASE(fields))
	        readfield = choice(0) + 1
;                fxbread,unit,photomet,readfield
                photomet = tbget(bin_header,all_data,readfield)
                print,'unpacking standard deviations',systime(0)
                For ib=0,9 do begin
                  itemp = data_list(ib,*)/256
                  n = 3
                  IF (strupcase(strtrim(product_type)) eq 'GPM' and $
                     (ib eq 4 or ib eq 5)) THEN n=4
                  IF (ib GT 7) THEN n = 1
                  dim=SIZE(itemp)
                  for k=0l,dim(2)-1 do $ 
                       CASE itemp(k) OF
                         0: dbsig(ib,k) = 10.0^(-n) 
                         255: dbsig(ib,k) = 10.0^(-n+4) 
                         else: dbsig(ib,k) = $
                                  10.0^((itemp(k) - 0.5)/63.5 - n)
                       ENDCASE
                  IF (ib LE 7) THEN dbsig(ib,*) = $
                           dbsig(ib,*)*ABS(photomet(ib,*))
                endfor
                data_list=dbsig
                photomet=0
                dbsig=0
        ENDELSE
  ENDIF

  plen=strlen(strtrim(product_type))

  IF (strupcase(strmid(product_type,0,7)) eq 'WEEKMAP') THEN BEGIN
       warm = strpos(strupcase(product_type),'WARM')
       ssdname=strupcase(dsfield)

       IF (ssdname eq 'STDDEV  ') THEN BEGIN
                   sz=size(data_list)
                   dbsig=fltarr(sz(1),sz(2))
                   print,'unpacking standard deviations'
                   if warm eq -1 then top_ib = 9 else top_ib = 3
                   For ib=0,top_ib do begin
                     itemp = data_list(ib,*)
                     dim=SIZE(itemp)
                     CASE ib OF
                         4: n=3
                         5: n=3
                         6: n=2
                         7: n=2
                         8: n=0
                         9: n=1
                         else: n=4
                     ENDCASE
                     for k=0l,dim(2)-1 do $ 
                       CASE itemp(k) OF
                         0: dbsig(ib,k) = 10.0^(-n) 
                         255: dbsig(ib,k) = 10.0^(-n+4) 
                         else: dbsig(ib,k) = $
                                  10.0^((itemp(k) - 0.5)/63.5 - n)
                       ENDCASE
                   endfor
                   data_list=dbsig
                   dbsig=0
       ENDIF
       IF (ssdname eq 'STOKESSD') THEN BEGIN
                   sz=size(data_list)
                   dbsig=fltarr(sz(1),sz(2))
                   n=4
                   print,'unpacking Stokes standard deviations'
                   For ib=0,5 do begin
                     itemp = data_list(ib,*)
                     dim=SIZE(itemp)
                     for k=0l,dim(2)-1 do $ 
                       CASE itemp(k) OF
                         0: dbsig(ib,k) = 10.0^(-n) 
                         255: dbsig(ib,k) = 10.0^(-n+4) 
                         else: dbsig(ib,k) = $
                                  10.0^((itemp(k) - 0.5)/63.5 - n)
                       ENDCASE
                   endfor
                   data_list=dbsig
                   dbsig=0
       ENDIF
  ENDIF
  IF (strupcase(strmid(product_type,plen-3,3)) eq 'E90') THEN BEGIN
       IF (strupcase(dsfield) eq 'STDDEV  ') THEN BEGIN
                   print,'unpacking standard deviations'
                   yabs=0
                   n=3
                   bandname=strupcase(strmid(product_type,0,3))
                   if  bandname eq 'B05' or bandname eq 'B06' then n=4
                   if  bandname eq 'B09' or bandname eq 'B10' then begin
                        n=1
                        yabs=1
                   endif
                   sz=size(data_list)
                   dbsig=fltarr(sz(1))
                   dim=SIZE(data_list)
                   for k=0l,dim(1)-1 do $ 
                       CASE data_list(k) OF
                         0: dbsig(k) = 10.0^(-n) 
                         255: dbsig(k) = 10.0^(-n+4) 
                         else: dbsig(k) = $
                                  10.0^((data_list(k) - 0.5)/63.5 - n)
                       ENDCASE
                   IF (yabs ne 1) THEN BEGIN
                     print,'reading in PHOTOMET data for sigma...',systime(0)
         	     choice = where('PHOTOMET' Eq STRUPCASE(fields))
                     readfield = choice(0) + 1
;                     fxbread,unit,photomet,readfield
                     photomet = tbget(bin_header,all_data,readfield)
                     dbsig = dbsig*ABS(photomet)
                   ENDIF
                   data_list=dbsig
                   dbsig=0
       ENDIF
       IF (STRUPCASE(dsfield) ne 'TIME    ') THEN BEGIN
                  print,'reading in TIME field...',systime(0)
                  choice = where('TIME    ' Eq STRUPCASE(fields))
	          readfield = choice(0) + 1
;                  fxbread,unit,time,readfield
                  time=tbget(bin_header,all_data,readfield)
       ENDIF ELSE time=data_list
       start_sec=min(time)+timediff
       stop_sec=max(time)+timediff
       start=timeconv(start_sec,infmt='sec',outfmt='v')
       stop=timeconv(stop_sec,infmt='sec',outfmt='v')
times: 
       print,'This file contains data from ',start,' to ',stop
       print,'Please select a time period (maximum of 6 months)'
       print,'Time format must be: DD-MMM-YYYY, e.g 01-jan-1990'
       stime=''
       etime=''
       read,'Enter start time: ',stime
       read,'Enter stop time or "max" to get 6 months from start' $
                     +' time: ',etime
       stime_sec=timeconv(stime,infmt='v',outfmt='sec')-timediff
       if (strupcase(etime) ne 'MAX') then $
                  etime_sec=timeconv(etime,infmt='v',outfmt='sec')-timediff $
       else etime_sec=stime_sec+maxsecs
       if (etime_sec-stime_sec gt maxsecs) then begin
                   print,'Time range was greater than 6 months.'
                   print,'Please select another time range.'
                   goto, times
       endif
       valid=where(time ge stime_sec and time le etime_sec)
       IF (valid(0) eq -1) THEN BEGIN
                   print,'The times selected are not within the proper range.'
                   print,'Please select another time range.'
                   goto, times
       ENDIF
       data_list=data_list(valid)
       if (STRUPCASE(dsfield) ne 'PIXEL_NO') then pixel=pixel(valid)
  ENDIF
;  fxbclose,unit
  all_data = 0
  dum=temporary(all_data)

; subset the DIRBE data by wavelength if that has been requested
  if subscr ne '' then begin
       pos = strpos(subscr,':')

       if pos gt 0 then begin
		     ; if we're dealing with ranges we took the numbers
		     ; in indexed from 1-10 as idl is 0 based we have to
		     ; subtract 1 from each of the index numbers so that
		     ; we use the proper portions of the array
	    	     num1 = fix(strmid(subscr,0,pos)) -1
	    	     len = strlen(subscr) - (pos+1)
	    	     num2 = fix(strmid(subscr,pos+1,len)) -1
		     data_list = data_list(num1:num2,*)
		     frequency = frequency(num1:num2)

       ;else we have a scalar value
       endif else begin
	    	     num = fix(subscr) - 1
		     data_list = data_list(num,*)
		     frequency = frequency(num)
       endelse
  endif

  if (STRUPCASE(dsfield) eq 'PIXEL_NO') then pixel = data_list

  print, 'beginning pixelization ',systime(0)
  if faceno ne 7 then begin
		     pixface, pixel, data_list, faceno=faceno, res=res
		     if !err eq -1 then goto, abort
		     pix2xy, pixel, res=res, data=data_list, raster=data, $
			     bad_pixval=badpixval, /face
  endif else begin
		  pix2xy, pixel, res=res, data=data_list, raster=data, $
		          bad_pixval=badpixval
  endelse

abort:
  data_list = 0
  pixel = 0
  end
;DISCLAIMER:
;
;This software was written at the Cosmology Data Analysis Center in
;support of the Cosmic Background Explorer (COBE) Project under NASA
;contract number NAS5-30750.
;
;This software may be used, copied, modified or redistributed so long
;as it is not sold and this disclaimer is distributed along with the
;software.  If you modify the software please indicate your
;modifications in a prominent place in the source code.  
;
;All routines are provided "as is" without any express or implied
;warranties whatsoever.  All routines are distributed without guarantee
;of support.  If errors are found in this code it is requested that you
;contact us by sending email to the address below to report the errors
;but we make no claims regarding timely fixes.  This software has been 
;used for analysis of COBE data but has not been validated and has not 
;been used to create validated data sets of any type.
;
;Please send bug reports to CGIS@ZWICKY.GSFC.NASA.GOV.



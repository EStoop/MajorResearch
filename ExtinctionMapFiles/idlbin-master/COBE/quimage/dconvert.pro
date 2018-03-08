;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    DCONVERT is a routine used to read one data format and write another
;
;DESCRIPTION:
;    DConvert will read FITS, IDL Save Sets, COBE SkyMaps, and UImage Objects
;    and output FITS, IDL Save Sets, and UImage Objects.
;
;    DCONVERT will read the data set according to input_format and write
;    the data set according to output_format
;
;    If the input_format is a file (CSM, CISS, FITS) then dsname must
;    contain a data set name if CSM then the fieldname must be filled in
;
;CALLING SEQUENCE:
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    param		I/O   	type		description
;    input_format       I      	string		CSM, FITS, SaveSet, etc
;    input_type         I	string		FIRAS, DIRBE, DMR, etc
;    dsname             I	string		full file specification
;    faceno             I	string		face number, 0-5, 7 if all
;    dsfield            I	string		full field name
;    subscr             I	string		subscript(s) (eg, '5:8')
;    weightf		I	string		name of file to write weights
;						array to/read it from
;    output_format      I	string		FITS, SaveSet, etc
;    outname            I	string		full file specification
;
;    All outputs from this routine are in the form of files or UImage
;    Objects.
;
;WARNINGS:
;
; Make sure that everything is in order, this is not a very forgiving
; program.
;
;EXAMPLE: 
;
; ____
;
;#
;COMMON BLOCKS: uimage_data, uimage_data2, info_3D
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;       Calls DIRBEIN to read in DIRBE FITS files
;
;MODIFICATION HISTORY
;
;	25-Nov-1992	SPR 10169	Kryszak	checks UPutData status
;	16-Dec-1992	SPR 10178	Kryszak put in pause when done
;	31-Dec-1992	SPR 10369	Kryszak added 3d meta-data when
;					retrieving UImage object
;	31-Dec-1992	SPR 10309	Kryszak added freq array to FITS
;					it's a kludge, but 
;	08-Jan-1993			Kryszak, removed some debug code
;	14-Jan-1993	SPR 10461	Kryszak removed extra XDisplay
;	03-Feb-1993	SPR 10493	Bazell changed order of Header
;					Keywords to conform to FITS std
;       11-Mar-1993     SPR 10672       Remove the PAUSE statement
;	26-Mar-1993			Ward, add the capability to read
;					FITS binary tables, clean up some
;					errors and wording of prompts
;       15-Apr-1993     SPR 10829       Ward,Newmark add capability to read
;                                       in proper weights along with signal.
;       28-Jun-1993     SPR 11101       Newmark added call to append_number
;                                       for all titles.
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;  SPR 11226 18 Aug 1993  Add large reprojections. J. Newmark
;  SPR 11311 14 Sep 1993  Write out weights in CISS. J. Newmark
;  SPR 11347 30 Sep 1993  Modify for DATAIN and DATAOUT J. Newmark
;  SPR 11375 13 Oct 1993  Modify for DIRBE PhotQual. J. Newmark
;  SPR 11491 10 Dec 1993  Add weight units field. J. Newmark
;  SPR 11759 17-May-1994   Ingest FIRAS PDS's. J. Newmark
;  SPR 11905 07-Sep-1994   Ingest DIRBE PDS's. J. Newmark
;  SPR 11981 08-Nov-1994  Read in skymap_info files.  J. Newmark 
; SPR 12140 16-Mar-1995   Ingest FIRAS line maps. J. Newmark
;
;  25-Aug-1997 J. Newmark - major upgrade - changed calls from
;                fxbopen/fxbread to readfits/tbget.
;  17-Nov-1997 J. Newmark - created subroutine FIRASIN for FIRAS input
;
;.TITLE
;DConvert
;-

Pro dconvert, input_format, input_type, dsname, fits_exten, faceno, $
              dsfield, subscr, weightf, output_format, outname,	$
	      data=data, weights=weights, frequency=frequency, $
	      instrume=instrume, badpixval=badpixval, res=res, $
	      units=units, wunits=wunits, rms=rms, coord=coord, $
              num_freq=num_freq, delta_nu=delta_nu, nu_zero=nu_zero, $
              dbsig=dbsig,dbqual=dbqual

; include the common block of all data
;
COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
COMMON uimage_data2,proj_map,proj2_map
COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
  freq3D_2,wt3d_0,wt3d_1,wt3d_2

; initialize some stuff
IF (!cgis_os EQ 'vms') THEN eod=']'
IF (!cgis_os EQ 'unix') THEN eod='/'
IF (!cgis_os EQ 'windows') THEN eod='\'
weight = 'N'
old_quiet=!quiet

; Save the fieldname for journaling
jdsfield = dsfield

;First, get the data into internal format from the existing format
; if User files selected then do that and return...
if input_type eq 'U' then Begin
    if input_format eq 'CSM' then Begin
	outmap = buildmap()		; call flexible i/o 

	; find out if we actually created anything!
	elem = n_elements(outmap)

	if (elem gt 0) then begin
	    uputdata, outmap, name=uname, error=error
	    If (!D.NAME EQ 'X' or !D.NAME EQ 'WIN') and (error EQ 0) Then Begin
			xdisplay,uname
	    endif else if (error eq 1) then begin
		message, /continue, 'No UImage Object Created.'
		message, /continue, 'Object is...'
		help, outmap
	    endif
	EndIf Else Begin
	    message, /continue, 'No map returned, no Object created.'
	Endelse
	goto, report
    EndIf                      ; if input_format = 'CSM'
EndIf			       ; if input_type = 'U'


;configured formats

CASE input_format OF
    'CSM': BEGIN
  	   Print, 'Reading COBE format Sky Map, '+dsname
           csmhdr, dsname, cname, description, cversion, cunits, channels, $
                   idxlev, bandfreq, cbandwidth, bandunits, bandnull

           csmdata, dsname, dsfield, subscr, faceno, idxlev(0), data
        
	;UNITS STRING SELECTION -------------------------------------
       
	 If (dsfield ne 'SIGNAL') or (dsfield ne 'SERROR') Then Begin
           ;units string needs to be retrieved

           If subscr ne '' Then Begin
              ; remember that bandno is 1 based (like FORTRAN)
              ; but IDL is 0 based,  should be changed in future to
              ; be consistent or list frequency
              bandno = Fix(subscr)
              EndIf $
           Else Begin
              bandno = 1
              EndElse

	   if cunits(0) eq -1 then begin
		if bandunits(0) le 0 then begin
		   message, /continue, 'Units defined for each record.'
		   message, /continue, 'See record definition for units field.'
		   message, /continue, 'Using field name for units.'
		   units = dsfield
		   endif $
		else begin
	       	   ;otherwise assume that units the same for each channel
		   ; this is a bad assumption, but necessary for UImage
		   units = UnitStr( bandunits(0) )
		   endelse
		endif $
	   else begin  ; units defined by scalar UNITS field
		units = UnitStr( cunits(0) )
		endelse

	   EndIf $
	Else Begin
	   message,/continue,'Units not available for this field.'
	   message,/continue,'Using field name.' 
	   units = dsfield
	EndElse

	; UNITS STRING SELECTION ----------------------------------------

	; make assumptions about ADB CSM (cobe sky maps)
        orient = 'R' ; always right handed T, skylooking cube
        coordinate_system = 'ECLIPTIC' 
        projection = 'SKY-CUBE'  
	badpixelval = 0.	

        instrume = input_type   ; for now, just first letter

        ; subscripted vars converted to scalars
        version = cversion(0)
	t = STRPOS( subscr, ':' )
	If t gt 0 Then Begin
           starti = FIX(STRMID(subscr,0,t))
	   endi   = FIX(STRMID(subscr,t+1,STRLEN(subscr)))
	   frequency=bandfreq(starti-1:endi-1) 
	   bandwidth=cbandwidth(starti-1:endi-1)
	   Endif $
	Else Begin ; no subscripts
	   dim3 = SIZE(data)
           If dim3(0) GT 3 Then Begin
              MESSAGE, 'Too many dimensions in data array.',/continue 
              goto, report 
              EndIf
           If dim3(0) EQ 3 Then Begin 
	      ; for 3D data, with whole array, but
	      ; remember that the header array always 
	      ; has 512 elements, not the same as the actual data
	      frequency = bandfreq(0:dim3(3)-1)
	      bandwidth = cbandwidth(0:dim3(3)-1)
	      Endif $
	   Else Begin
	      ; not 3D data, no subscripts
              frequency = bandfreq(bandno-1)
              bandwidth = cbandwidth(bandno-1)
	      EndElse
	   EndElse


        if subscr eq '' then begin
           title = cname(0)+' -- '+dsfield
           endif $
        else begin
           title = cname(0)+' -- '+dsfield+'('+subscr+')'
        endelse

        END
    ;
    'CISS': BEGIN
        print, 'Restoring save set: ',dsname
	uputdata, dsname, name=uname, error=error
  	If (!D.NAME EQ 'X' OR !D.NAME EQ 'WIN') and (error EQ 0) Then Begin
	   xdisplay,uname
	   endif $
	else if (error eq 1) then begin
	   message, /continue, 'No UImage object created from saveset.'
	   endif
  	
        goto, report 
	END
    ;
    'FITS': BEGIN
        Print, 'Reading FITS format file, '+dsname
	print, 'beginning fits processing ',systime(0)

	if (fits_exten Eq 'BINTABLE') then begin
	   ; get the primary FITS header
	   header = headfits(dsname)

           instrume  = sxpar(header, 'INSTRUME') 
           product_type = sxpar(header, 'PRODUCT')
           if string(product_type) eq 'HLIN-MAP' then instrume ='HLIN-MAP'
           if string(product_type) eq 'LLIN-MAP' then instrume ='LLIN-MAP'
           if string(product_type) eq 'PIXINFO6' or $
                string(product_type) eq 'PIXINFO9' then instrume='PXINFO'
	   weight = 'N'

	   ; setup some default values
	   if strtrim(instrume) eq 'PXINFO' then begin
	       res       = sxpar(header, 'PIXINDEX')
	       coord     = sxpar(header, 'COORDSYS')
	       if coord eq 'E' then coordinate_system = 'Ecliptic'
	       if coord eq 'C' then coordinate_system = 'Equatorial'
	       if coord eq 'G' then coordinate_system = 'Galactic'
	       if coord eq 'S' then coordinate_system = 'Special'
	       title = strmid(dsname,0,strpos(dsname,'.'))
	       ; read in the header for the fits BINARY TABLE extension
;	       fxbopen,unit,dsname,1,bin_header
               all_data = readfits(dsname,bin_header,/ext)
	       fields = fxpar(bin_header,'TTYPE*',count=fieldcnt)
	       unitlist  = fxpar(bin_header,'TUNIT*',count=unitcnt)
	       choice = where(dsfield Eq fields)
	       units =  unitlist(choice(0))

;	       IDL arrays start at 0, fits columns start at 1 so we have
;	       to add 1 to the variable choice to point to the correct one
	       readfield = choice(0) + 1

;	       fxbread,unit,data_list,readfield
;	       if (choice(0) ne 0) then fxbread,unit,pixel,'QSPIXEL'
;	       fxbclose,unit
               data_list = tbget(bin_header,all_data,readfield)
               if (choice(0) ne 0) then pixel = $
                         tbget(bin_header,all_data,'QSPIXEL')
               all_data = 0
               dum = temporary(all_data) 

	       if choice(0) eq 0 then pixel = data_list

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

;	       delete the data arrays that are no longer needed
	       data_list = 0
	       pixel = 0
           endif

	   if strtrim(instrume) eq 'DMR' then begin
	       badpixval = sxpar(header, 'BADPIX' )
	       res       = sxpar(header, 'PIXRESOL')
	       rms       = sxpar(header,  'RMS_FLGT')
               if rms eq 0 then rms = sxpar(header,  'RMS')
	       coord     = sxpar(header, 'COORDSYS')
	       if coord eq 'E' then coordinate_system = 'Ecliptic'
	       if coord eq 'C' then coordinate_system = 'Equatorial'
	       if coord eq 'G' then coordinate_system = 'Galactic'
	       if coord eq 'S' then coordinate_system = 'Special'

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
	          while strlen(dsfield) lt 8 do dsfield = dsfield + ' '
	       endif

	       if dsfield Eq 'SIGNAL+WEIGHTS' then begin
		   dsfield = 'SIGNAL  '
		   weight = 'Y'
	           wfield= 'N_OBS   '
		   wtitle = 'weight_of_'+strmid(dsname,0,strpos(dsname,'.'))
               endif

	       choice = where(strupcase(dsfield) Eq strupcase(fields))
	       units =  unitlist(choice(0))

               IF (weight EQ 'Y') THEN BEGIN
                 wchoice=where(wfield EQ fields)
                 wunits = unitlist(wchoice(0))
               ENDIF

;	       IDL arrays start at 0, fits columns start at 1 so we have
;	       to add 1 to the variable choice to point to the correct one
	       readfield = choice(0) + 1

;	       fxbread,unit,data_list,readfield
;	       if (choice(0) ne 0) then fxbread,unit,pixel,'PIXEL'
               data_list = tbget(bin_header,all_data,readfield)
               if (choice(0) ne 0) then pixel = $
                         tbget(bin_header,all_data,'PIXEL') 

	       ; bring in the weighting array
	       if weight eq 'Y' then begin
;		  fxbread,unit,tweights,'N_OBS'
                  tweights = tbget(bin_header,all_data,'N_OBS')
		  tweights = tweights/rms^2
	       endif

;	       fxbclose,unit
               all_data = 0
               dum = temporary(all_data)

	       if choice(0) eq 0 then pixel = data_list

	       print, 'beginning pixelization ',systime(0)
	       if faceno ne 7 then begin
		  if weight eq 'Y' then begin
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
		  if weight eq 'Y' then $
		     pix2xy, pixel, res=res, data=tweights, raster=weights
	       endelse

;	       delete the data arrays that are no longer needed
	       data_list = 0
	       pixel = 0
	       if weight ne 'N' then tweights = 0
	   Endif                ; end DMR instrument processing


	   If strtrim(instrume) eq 'FIRAS' or strtrim(instrume) eq 'HLIN-MAP' $
             or strtrim(instrume) eq 'LLIN-MAP' then $
             firasin,dsname,header,dsfield,product_type,subscr,faceno,$
               badpixval,data,frequency,res,frequnits,title,units,weight,$
               weights,wunits,wtitle,instrume

	   If strtrim(instrume) eq 'DIRBE' then $
             dirbein,dsname,header,dsfield,product_type,subscr,faceno,$
               badpixval,data,frequency,res,frequnits,title,units,weight,$
               dbsig,dbqual

	EndIf Else begin
	   ; get the primary FITS header
	   data = readfits( dsname, header, /NOSCALE )	

	   ; these are proper FITS keywords
	   units = sxpar( header, 'BUNITS' ) 
           instrume = sxpar ( header, 'INSTRUME') 
	   badpixval = sxpar( header, 'BADPIX' )

	   ; defaults (no UImage header records)
	   version = 'none'
	   faceno = 7
	   frequency = 0.
	   bandwidth = 0.
       	   title = dsname

	   freqrecs = where(strpos(header, 'FREQ') ge 0)
	   if freqrecs(0) ne -1 then frequency = sxpar(header, 'FREQ*')

	   histrecs = where (strpos( header, 'HISTORY') ge 0)

	   if histrecs(0) ne -1 then begin
		header =    header ( histrecs )
     	        ; then get the values
     		; version
		spos = where (strpos( header, 'VERSION') ge 0) 
		If spos(0) ge 0 Then Begin
			aversion =   header ( spos ) 
			version = STRMID( aversion(0), $
			  STRPOS(aversion(0), '=')+1, STRLEN(aversion(0)) )
			version = STRTRIM( version, 2 )
           		EndIf

     		; coord
		spos = where (strpos( header, 'COORD') ge 0) 
		If spos(0) ge 0 Then Begin
			acoord =   header ( spos ) 
			coord = STRMID( acoord(0), $
			  STRPOS(acoord(0), '=')+1, STRLEN(acoord(0)) )
			coord = STRTRIM( coord, 2 )
           		EndIf

		; faceno
       		spos = where (strpos( header, 'FACENO') ge 0) 
       		If spos(0) ge 0 Then Begin
       			afaceno =   header ( spos ) 
       			status = strnumber( STRMID( afaceno(0), $
       			  STRPOS(afaceno(0), '=')+1, STRLEN(afaceno(0)) ), $
				faceno )
			EndIf

		; res
       		spos = where (strpos( header, 'RES') ge 0) 
       		If spos(0) ge 0 Then Begin
       			afaceno =   header ( spos ) 
       			status = strnumber( STRMID( afaceno(0), $
       			  STRPOS(afaceno(0), '=')+1, STRLEN(afaceno(0)) ), $
				res )
			EndIf

		; rms
       		spos = where (strpos( header, 'RMS') ge 0) 
       		If spos(0) ge 0 Then Begin
       			afaceno =   header ( spos ) 
       			status = strnumber( STRMID( afaceno(0), $
       			  STRPOS(afaceno(0), '=')+1, STRLEN(afaceno(0)) ), $
				rms)
			EndIf

		; nu_zero
       		spos = where (strpos( header, 'NU_ZERO') ge 0) 
       		If spos(0) ge 0 Then Begin
       			afaceno =   header ( spos ) 
       			status = strnumber( STRMID( afaceno(0), $
       			  STRPOS(afaceno(0), '=')+1, STRLEN(afaceno(0)) ), $
				nu_zero)
			EndIf

		; delta_nu
       		spos = where (strpos( header, 'DELTA_NU') ge 0) 
       		If spos(0) ge 0 Then Begin
       			afaceno =   header ( spos ) 
       			status = strnumber( STRMID( afaceno(0), $
       			  STRPOS(afaceno(0), '=')+1, STRLEN(afaceno(0)) ), $
				delta_nu )
			EndIf

 		; num_freq
       		spos = where (strpos( header, 'NUM_FREQ') ge 0) 
       		If spos(0) ge 0 Then Begin
       			afaceno =   header ( spos ) 
       			status = strnumber( STRMID( afaceno(0), $
       			  STRPOS(afaceno(0), '=')+1, STRLEN(afaceno(0)) ), $
				num_freq )
			EndIf

		; bandwidth
		spos = where (strpos( header, 'BANDWIDTH') ge 0) 
		If spos(0) ge 0 Then Begin
			bwidth = header ( spos ) 
			status = strnumber( STRMID( bwidth(0), $
			   STRPOS(bwidth(0), '=')+1, $
			   STRLEN(bwidth(0)) ), bandwidth )
	        EndIf

		; title
       		spos = where (strpos( header, 'TITLE' ) ge 0) 
       		atitle =     header ( spos )
       		If spos(0) ge 0 Then Begin
       			title = STRMID( atitle(0), $ 
       				STRPOS( atitle(0), '=' )+1, STRLEN( atitle(0)) )
       			title = STRTRIM( title, 2 )
       		EndIf
	
	   EndIf ; proper history records

     	   If faceno eq 7 Then Begin  ;it's a full map
	      datasize = SIZE( data )

	      If 2*datasize(1) eq 3*datasize(2)  Then Begin
		data = sixunpack( data, badval=badpixval ) 
		EndIf $
              Else Begin 
		If (4*datasize(2) ne 3*datasize(1)) and $
    		   ((datasize(1) ne 512) or (datasize(2) ne 256)) Then Begin
		   Message, /continue, 'Invalid data array size.'
		   Return
		EndIf
              EndElse
	   EndIf

	   ; get the weights if they were asked for
	   if weightf eq 'yes' then begin
	      header = headfits(dsname)

	      ; get the path for the weights file (assumed to be in the same
	      ; directory as the datafile
	      dfile = fname(dsname,path=path)
	      if path eq '' then $
                 filename = sxpar(header,'WEIGHTS') $
		 else filename = path+sxpar(header,'WEIGHTS')

	      wunits = sxpar(header,'BUNITS')
	      weights = readfits(filename, wheader, /NOSCALE)

       	      spos = where(strpos(wheader, 'TITLE') ge 0) 
       	      atitle = wheader (spos)
       	      If spos(0) ge 0 Then Begin
       		 wtitle = STRMID( atitle(0), $ 
       			 STRPOS( atitle(0), '=' )+1, STRLEN( atitle(0)) )
       		 wtitle = STRTRIM( wtitle, 2 )
       	      EndIf

	      datasize = SIZE(weights)
	      if datasize(0) eq 3 then weight = '3' $
		 else weight = 'Y'

     	      If faceno eq 7 Then Begin  ;it's a full map
	         If 2*datasize(1) eq 3*datasize(2)  Then Begin
		    weights = sixunpack( weights, badval=badpixval ) 
		 EndIf Else Begin 
		    If (4*datasize(2) ne 3*datasize(1)) and $
    		       ((datasize(1) ne 512) or (datasize(2) ne 256)) Then Begin
		       Message, /continue, 'Invalid weights array size.'
		       Return
		    EndIf
                 EndElse
	      Endif
	   EndIf

       EndElse  
    END

    ;
    'UIMAGE': BEGIN

        i = Execute( 'title='+dsname+'.title' )
	Print, 'Using object called '+ title
        ;
	; These items are common to any UImage format data
	; including MAP, FACE, ZOOMED
        i = Execute( 'badpixval='+dsname+'.badpixval' )
        i = Execute( 'units='+dsname+'.units' )
        i = Execute( 'faceno ='+dsname+'.faceno' )
        i = Execute( 'instrume='+dsname+'.instrume' )
        i = Execute( 'version='+dsname+'.version' )
        i = Execute( 'orient ='+dsname+'.orient' )
        i = Execute( 'projection='+dsname+'.projection' )
        i = Execute( 'coordinate_system='+dsname+'.coordinate_system' )

        bandno = 0
        frequency= 0
        bandwidth= 0
        parent = 0
        scale_max = 0
        scale_min = 0

	; if weights are supposed to be written out then we grab them here
	; to be used later on in the program (the output section)
	if ((weightf ne '') OR (output_format EQ 'CISS')) then begin
	   weights=0
	   namew = ''
	   first_letter=STRMID(dsname,0,1)
	   j=EXECUTE('linkweight= ' +dsname+ '.linkweight')
	   ;
	   ; Handle Case where selected object is a 3-dim one.

	   IF first_letter EQ 'O' THEN BEGIN
	      str_index3d=STRMID(dsname,9,1)
	      index3d=FIX(str_index3d)
	      j=EXECUTE('havewt= wt3d_' + str_index3d)
	      sz=SIZE(havewt)
	      IF sz(0) EQ 3 THEN BEGIN
	         PRINT,'This object has a 3-dim weight array associated with it.'
	         weights=havewt
		 wunits = 'none'
		 wtitle = 'weight_of_'+dsname
	      ENDIF
	   ENDIF

	   IF first_letter NE 'G' THEN BEGIN 
	      IF linkweight NE -1 THEN BEGIN
	         namew=get_name(linkweight)
	         j = EXECUTE('wtitle = ' + namew + '.title')
                 i = Execute( 'wunits='+namew+'.units' )
	         PRINT, 'Weights will be taken from ' + wtitle
	         j = EXECUTE('weights = ' + namew + '.data')
	      ENDIF
	   ENDIF   ;ELSE BEGIN
	;
	;      Handle case where selected object is a 1-dim graph
	;      j=EXECUTE('num = ' +dsname+ '.num')
	;      IF (linkweight NE -1) THEN BEGIN
	;          namew=get_name(linkweight)
	;          j = EXECUTE('wtitle = ' + namew + '.title')
	;          PRINT, 'Weights will be taken from ' + wtitle
	;          j = EXECUTE('weights = ' + namew + '.data(0:num-1,1)')
	;      ENDIF
	;   ENDELSE
	endif

	Case 1 Of

    ( STRPOS( dsname, 'MAP' ) ne -1 ) or (STRPOS( dsname, 'PROJ' ) ne -1 ) $
           or (STRPOS( dsname, 'FACE' ) ne -1 ) : Begin

 	    i = Execute( 'data='+dsname+'.data' )
            i = Execute( 'bandno ='+dsname+'.bandno' )
            i = Execute( 'frequency='+dsname+'.frequency' )
            i = Execute( 'bandwidth='+dsname+'.bandwidth' )
            i = Execute( 'scale_max ='+dsname+'.scale_max' )
            i = Execute( 'scale_min ='+dsname+'.scale_min'  )

	   End

    	(STRPOS( dsname, 'ZOOMED' ) ne -1) : begin

	   i = Execute( 'winorig='+dsname+'.win_orig' )
	   
	   i = Execute( 'start_x=STRING('+dsname+'.start_x)' )
	   i = Execute( 'stop_x=STRING('+dsname+'.stop_x)' )
	   i = Execute( 'start_y=STRING('+dsname+'.start_y)' )
	   i = Execute( 'stop_y=STRING('+dsname+'.stop_y)' )

 	   i = Execute( 'naxis=SIZE( ' + dsname + '.data )' ) 
	   if start_x lt 0 then start_x = 0
	   if start_y lt 0 then start_y = 0
	   if stop_x gt naxis(1) then stop_x = naxis(1)
	   if stop_y gt naxis(2) then stop_y = naxis(2)

	   if start_x gt naxis(1) then start_x = naxis(1)
	   if start_y gt naxis(2) then start_y = naxis(2)
	   if stop_x lt 0 then stop_x = 0
	   if stop_y lt 0 then stop_y = 0

	   xdata = start_x + ':' + stop_x
	   ydata = start_y  + ':' + stop_y

	   dsname= get_name(winorig)
  
 	    i = Execute( 'data='+dsname+'.data('+ xdata + ',' + ydata + ')' )
            i = Execute( 'bandno ='+dsname+'.bandno' )
            i = Execute( 'frequency='+dsname+'.frequency' )
            i = Execute( 'bandwidth='+dsname+'.bandwidth' )
            i = Execute( 'scale_max ='+dsname+'.scale_max' )
            i = Execute( 'scale_min ='+dsname+'.scale_min'  )

	   End

	(STRPOS( dsname, 'OBJECT3D' ) ne -1 ) : begin
	    n = STRPOS( dsname, '(' ) + 1
	    n = FIX( STRMID( dsname, n, 1 ) )
	    rdsname = 'object3d(' + strtrim(string(n),2) + ')'

            i = Execute( 'data=data3D_'+STRTRIM(STRING(n),2) )
            i = Execute( 'frequency=freq3D_'+STRTRIM(STRING(n),2) )
            i = Execute( 'badpixval='+rdsname+'.badpixval' )
            i = Execute( 'units='+rdsname+'.units' )
            i = Execute( 'wunits='+rdsname+'.frequnits' )
            i = Execute( 'faceno ='+rdsname+'.faceno' )
            i = Execute( 'instrume='+rdsname+'.instrume' )
            i = Execute( 'version='+rdsname+'.version' )
            i = Execute( 'orient ='+rdsname+'.orient' )
            i = Execute( 'projection='+rdsname+'.projection' )
            i = Execute( 'coordinate_system='+rdsname+'.coordinate_system' )
	   End

	Else: Begin
		Message, /continue, 'Invalid data object selected.'
		return
		End

	EndCase ; of data object type (MAP,FACE,ZOOMED, 3D)
    End ; case of UImage

    'UIDL': Begin
            End
       
    ; default
    ELSE: BEGIN
        !quiet=0
        Message, /informational, 'Invalid data set type selected.'
        !err = 1
        !error = 1
        !quiet=old_quiet
        goto, report 
        END
    ENDCASE

CASE instrume OF
  'B' : instrume='DIRBE'
  'D' : instrume='DMR'
  'F' : instrume='FIRAS'
  ELSE:
ENDCASE

CASE output_format OF
    'CSM': BEGIN
        !quiet=0
        message, 'COBE Sky Map Format is not supported for output.', $
            /informational
        message, 'No file written.', /informational
        !quiet=old_quiet
        END

    'CISS': BEGIN
        print, 'Saving data in COBE IDL Save Set, ', outname
    
	If ((faceno lt 0) or (faceno ge 6)) then begin  ;full map 
	   datasize = SIZE( data )			;or some such

	   If 3*datasize(1) eq 4*datasize(2) Then Begin
		data = sixpack( data )
                IF (N_ELEMENTS(weights) GT 1) THEN weights=sixpack(weights)
		EndIf
	   EndIf
 
	if n_elements(title) eq 0 then title = 'COBE data for instrument ' $
           +instrume
	pos = strpos(title,'[')
	if pos gt 0 then title = strmid(title,0,pos-1)
        if n_elements(wunits) eq 0 then wunits='UNKNOWN'
        frequnits=wunits
;
; Check for weights
; 	
        IF (N_ELEMENTS(weights) GT 1) THEN $
	save, $     
                data, $
                weights, $
                badpixval, $
                bandno, $
                frequency, $
                bandwidth, $
                units, $
                frequnits, $
                faceno, $
                title, $
                instrume, $
                version, $
                orient, $
                projection, $
                coordinate_system, $
                scale_max, $
                scale_min, $
		filename=outname $
	ELSE save, $     
                data, $
                badpixval, $
                bandno, $
                frequency, $
                bandwidth, $
                units, $
                frequnits, $
                faceno, $
                title, $
                instrume, $
                version, $
                orient, $
                projection, $
                coordinate_system, $
                scale_max, $
                scale_min, $
		filename=outname 

        END
    ;
    'FITS': BEGIN
	Print, 'Writing data in FITS file, '+outname
	if (n_elements(faceno) EQ 0) OR (faceno eq -1) then faceno = 7

	if n_elements(version) eq 0 then version = 'none'
	if n_elements(frequency) eq 0 then frequency = 0.
	if n_elements(bandwidth) eq 0 then bandwidth = 0.

	naxis = size(data)

	; add FITS header parameters one-at-a-time
	sxaddpar, header, 'SIMPLE', 'T', ' Written by IDL:  ' + !STIME
	sxaddpar, header, 'BITPIX', -32
     	sxaddpar, header, 'NAXIS',  naxis(0)
	sxaddpar, header, 'NAXIS1', naxis(1)
	sxaddpar, header, 'NAXIS2', naxis(2)

	if naxis(0) GT 2 then begin
	   sxaddpar, header, 'NAXIS3', STRING(naxis(3))
           if (n_elements(frequency) ne naxis(3)) then $
              frequency=indgen(naxis(3))
	   for i = 0, naxis(3)-1 do begin
	       sxaddpar, header, strcompress('FREQ'+string(i+1),/remove_all),$
                         frequency(i)
	   endfor	
	endif else begin
	   sxaddpar, header, 'FREQ1',frequency
	endelse

	; put in a link to the weights file, if one exists
	if weightf ne '' then begin
	   wfile = fname(weightf)
	   sxaddpar, header, 'WEIGHTS', wfile
	endif
	sxaddpar, header, 'BADPIX',  badpixval
	sxaddpar, header, 'BUNITS', string(units)
	sxaddpar, header, 'TELESCOP', 'COBE'
	sxaddpar, header, 'INSTRUME', instrume
	sxaddpar, header, 'HISTORY', 'VERSION='+version     
        sxaddpar, header, 'HISTORY', 'FACENO='+string(faceno)
        sxaddpar, header, 'HISTORY', 'BANDWIDTH='+string(bandwidth)
        if n_elements(res) ne 0 then $
           sxaddpar, header, 'HISTORY', 'RES='+string(res)
        if n_elements(rms) ne 0 then $
           sxaddpar, header, 'HISTORY', 'RMS='+string(rms)
        if n_elements(coord) ne 0 then $
           sxaddpar, header, 'HISTORY', 'COORD='+coord
        if n_elements(nu_zero) ne 0 then $
           sxaddpar, header, 'HISTORY', 'NU_ZERO='+string(nu_zero)
        if n_elements(delta_nu) ne 0 then $
           sxaddpar, header, 'HISTORY', 'DELTA_NU='+string(delta_nu)
        if n_elements(num_freq) ne 0 then $
           sxaddpar, header, 'HISTORY', 'NUM_FREQ='+string(num_freq)
        if n_elements(wunits) ne 0 then $
           sxaddpar, header, 'HISTORY', 'FREQUNITS='+wunits
	if n_elements(title) eq 0 then title = 'COBE data for instrument ' $
           +instrume
	pos = strpos(title,'[')
	if pos gt 0 then title = strmid(title,0,pos-1)
	sxaddpar, header, 'HISTORY', 'TITLE='+title
	sxaddpar, header, 'HISTORY', 'Written by UImage DATA IO'
	sxaddpar, header, 'DATE', systime()

	writefits, outname, data, header

	; handle writing out weights if present and requested
	if weightf ne '' then begin
	   Print, 'Writing weights in FITS file, '+weightf
	   naxis = size(weights)

	   ; add FITS header parameters one-at-a-time
	   sxaddpar, wheader, 'SIMPLE', 'T', ' Written by IDL:  ' + !STIME
	   sxaddpar, wheader, 'BITPIX', -32
     	   sxaddpar, wheader, 'NAXIS',  naxis(0)
	   sxaddpar, wheader, 'NAXIS1', naxis(1)
	   sxaddpar, wheader, 'NAXIS2', naxis(2)

	   if naxis(0) GT 2 then begin
	      sxaddpar, wheader, 'NAXIS3', STRING(naxis(3))
              if (n_elements(frequency) ne naxis(3)) then $
                   frequency=indgen(naxis(3))
	      for i = 0, naxis(3)-1 do begin
	          sxaddpar, wheader, strcompress('FREQ'+string(i+1),$
                            /remove_all), frequency(i)
	      endfor
	   endif else begin
	      sxaddpar, wheader, 'FREQ1',frequency
	   endelse

	   dfile = fname(outname,path=dpath)
	   sxaddpar, wheader, 'DATAFILE', dfile
	   sxaddpar, wheader, 'BADPIX',  badpixval
	   sxaddpar, wheader, 'BUNITS', string(wunits)
	   sxaddpar, wheader, 'TELESCOP', 'COBE'
	   sxaddpar, wheader, 'INSTRUME', instrume
	   sxaddpar, wheader, 'HISTORY', 'VERSION='+version     
           sxaddpar, wheader, 'HISTORY', 'FACENO='+string(faceno)
           sxaddpar, wheader, 'HISTORY', 'BANDWIDTH='+string(bandwidth)

           if n_elements(res) ne 0 then $
              sxaddpar, header, 'HISTORY', 'RES='+string(res)
           if n_elements(rms) ne 0 then $
              sxaddpar, header, 'HISTORY', 'RMS='+string(rms)
           if n_elements(coord) ne 0 then $
              sxaddpar, header, 'HISTORY', 'COORD='+coord
           if n_elements(nu_zero) ne 0 then $
              sxaddpar, header, 'HISTORY', 'NU_ZERO='+string(nu_zero)
           if n_elements(delta_nu) ne 0 then $
              sxaddpar, header, 'HISTORY', 'DELTA_NU='+string(delta_nu)
           if n_elements(num_freq) ne 0 then $
              sxaddpar, header, 'HISTORY', 'NUM_FREQ='+string(num_freq)
           if n_elements(wunits) ne 0 then $
              sxaddpar, header, 'HISTORY', 'FREQUNITS='+wunits

	   ; be sure that wtitle is setup and ok.
	   if n_elements(wtitle) eq 0 then wtitle = 'weights for '+dfile
	   pos = strpos(wtitle,'[')
	   if pos gt 0 then wtitle = strmid(wtitle,0,pos-1)

	   sxaddpar, wheader, 'HISTORY', 'TITLE='+wtitle
	   sxaddpar, wheader, 'HISTORY', 'Written by UImage DATA IO'
	   sxaddpar, wheader, 'DATE', systime()

	   writefits, weightf, weights, wheader

	endif
        END
    ;
    'UIMAGE': BEGIN
	Print, 'Placing data into UImage object. ',systime(0)
        !err = 1
        !error = 1
        dim3 = SIZE(data)
        IF dim3(0) GT 3 THEN BEGIN
            MESSAGE, 'Too many dimensions in data array.',/continue 
            goto, report 
        ENDIF

	if faceno eq 7 then faceno = -1

        IF dim3(0) EQ 3 THEN Begin 
	   wlink = -1
	   if weight eq 'Y' then begin
	       wtitle=append_number(wtitle)
               uname = setup_image( $
                 data=weights, $
		 hidden=1, $
                 badpixval=0, $
                 bandno = bandno, $
                 frequency=freqency, $
                 bandwidth=bandwidth, $
                 units=wunits, $
                 faceno = faceno, $
                 title=wtitle, $
                 instrume=instrume, $
                 version=version, $
                 orient = orient, $
                 projection=projection, $
                 coordinate_system=coordinate_system $
                 )
		 j = execute('wlink = '+uname+'.window')
	   endif

	   if weight eq '3' then begin
             title=append_number(title)
             uname = setup_object3d( $
               data=data, $
               badpixval=badpixval, $
               frequency=frequency, $
	       frequnits=frequnits, $
               units=units, $
               wgtunits=wunits, $
               faceno = faceno, $
               title=title, $
	       linkweight=wlink, $
	       weight3d=weights, $	       
               instrume=instrume, $
               version=version, $
               orient = orient, $
               projection=projection, $
               coordinate_system=coordinate_system)
	   endif

	   if weight ne '3' then begin
             title=append_number(title)
             uname = setup_object3d( $
               data=data, $
               badpixval=badpixval, $
               frequency=frequency, $
	       frequnits=frequnits, $
               units=units, $
               faceno = faceno, $
               title=title, $
	       linkweight=wlink, $
               instrume=instrume, $
               version=version, $
               orient = orient, $
               projection=projection, $
               coordinate_system=coordinate_system)
	   Endif

  	   IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) $
              THEN show_object3d, uname

	Endif Else Begin
	     ; handle 2D images
	     wlink = -1

	     if weight Eq 'Y' or weight eq '3' then begin
	       wtitle=append_number(wtitle)
               uname = setup_image( $
                 data=weights, $
		 hidden=1, $
                 badpixval=0, $
                 bandno = bandno, $
                 frequency=freqency, $
                 bandwidth=bandwidth, $
                 units=wunits, $
                 faceno = faceno, $
                 title=wtitle, $
                 instrume=instrume, $
                 version=version, $
                 orient = orient, $
                 projection=projection, $
                 coordinate_system=coordinate_system $
                 )
		 j = execute('wlink = '+uname+'.window')
	      endif

              title=append_number(title)
              uname = setup_image( $
                data=data, $
                badpixval=0, $
                bandno = bandno, $
                frequency=freqency, $
                bandwidth=bandwidth, $
                units=units, $
                faceno = faceno, $
                title=title, $
		linkweight=wlink, $
                instrume=instrume, $
                version=version, $
                orient = orient, $
                projection=projection, $
                coordinate_system=coordinate_system)

  	    	IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) $
                   THEN xdisplay,uname
              EndElse
	END
           
'UIDL': BEGIN
	goto,report
    END

    ; default, oops
    ELSE: BEGIN
        !quiet=0
        Message, /informational, 'Invalid data set type selected.'
        !err = 1
        !error = 1
        !quiet=old_quiet
        goto, report 
        END
    ENDCASE


abort:
report:

if input_type eq 'U' then return
if input_format eq 'CISS' then return
if output_format eq 'UIDL' then return
if input_format eq 'UIDL' then return

;journaling

diojrnl, $
	input_format,   $
        input_type,     $
	instrume, 	$
        dsname,         $
	title, 		$
        faceno,         $
        jdsfield,       $
        subscr,         $
        output_format,  $
        outname

;	delete various local data arrays
	data = 0
	if (weight ne 'N' or weightf ne '') then weights = 0
	frequency = 0

END
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



;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    DATAIN is a routine which reads in FITS and other format files in UIDL
;
;DESCRIPTION:
;    Datain will read FITS, IDL Save Sets, and output FITS or IDL Save Sets.
;    It should be noted that some of the parameters in the list below are
;    instrument specific. This procedure is meant to be used for fields
;    within a FITS file which have been pixelized.
;
;WARNINGS:
;   1) If the input file is an IDL Save Set then the output will be sent
;    directly to an UIMAGE object. If the user wants to get the save
;    set directly to UIDL command line simply use the RESTORE command.
;   2) This routine is to be used for data fields which have been
;      pixelized. Some checking is done for FIRAS PDS's file types
;      to see if none of the fields are pixelized. If this is true
;      then the routine FIRASMOD is called to read in the data. 
;
;      IF YOU TRY AND READ IN A NON-PIXELIZED FIELD FOR A FILE THAT
;      CONTAINS SOME PIXELIZED FIELD THIS ROUTINE WILL BLOW UP DURING
;      THE PIXELIZATION.
;
;CALLING SEQUENCE: datain, dsname, [,format=format] [,intype=intype]
;                    [,face=face] [,dsfield=dsfield] [,subscr=subscr]
;                 [,data=data] [,weights=weights] [,frequency=frequency]
;                 [,instrume=instrume] [,badpixval=badpixval] [,res=res]
;                 [,units=units] [,wunits=wunits] [,rms=rms] [,coord=coord]
;                 [,nu_zero=nu_zero] [,delta_nu=delta_nu] [,num_freq=num_freq]
;                 [,/dbsig] [,/dbqual] [,/simple] [,pixel=pixel]
;                 [,startpix=startpix] [,stoppix=stoppix] [,indfile=indfile]
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    param		I/O   	type		description
;    dsname              I      string		full name of file for input
;    format=format       I      string		FITS or CISS (default to FITS)
;    intype=intype       I      string          DIRBE, FIRAS, DMR, or XAB, not
;                                               needed for most functions
;    face=face           I      integer         face # 0-5, not needed for full
;                                               skymap 
;    dsfield=dsfield     I      string          field to read in (if FITS BE)
;    subscr=subscr       I      string          subscripts to read in if FITS BE
;                                               format '4:7' (note: fields start
; 						at zero (1).
;    dbsig               I      integer         return sigma frm photqual(dirbe)
;    dbqual              I      integer         return flag frm photqual(dirbe)
;    simple              I      integer         set flag if read in simple FITS
;                                               non-skymap, data
;    data=data           O      float           array for data to be returned to
;    weights=weights     O      float           array for weights
;    frequency=frequency O      float           array for frequency info
;    instrume=instrume   O      string          the instrument the data is for
;    badpixval=badpixval O      float		"no data" indicator
;    res=res             O      integer		SKYMAP resolution
;    units=units         O      string		units for data values
;    wunits=wunits       O      string		units for weights
;    rms=rms             O      float		DMR specific
;    coord=coord         O      string          DMR specific
;    nu_zero=nu_zero     O      float           FIRAS specific
;    delta_nu=delta_nu   O      float           FIRAS specific
;    num_freq=num_freq   O      float           FIRAS specific
;    pixel=pixel         O      long            DIRBE specific -output pixel
;    startpix=startpix   I      long            DIRBE specific -start pixel
;    stoppix=stoppix     I      long            DIRBE specific -end pixel-
;                                                 only for DCAF, DSZA 
;    indfile=indfile     I      string          DIRBe specific -Index file
;
;    DMR specific keywords
;    ---------------------
;    rms, coord
;
;    FIRAS specific keywords
;    -----------------------
;    delta_nu, nu_zero, num_freq
;
;    DMR and FIRAS specific keywords
;    -------------------------------
;    weights
;    wunits
;
;    DIRBE and FIRAS specific keywords
;    ---------------------------------
;    frequency, dbsig
;    subscr, dbqual
;
;EXAMPLE: 
;         datain, 'cgis_fits:dmr_skymap_31a.fits',dsfield='SIGNAL',data=data,$
;                 instrume=instrume
;
;         datain, 'cgis_fits:dmr_skymap_53a.fits',dsfield='SIGNAL+WEIGHTS',$
;                 data=data,instrume=instrume
;
;         datain, 'cgis_fits:fip_sky_lhs.fits',dsfield='SIGNAL',data=data,$
;                 frequency=frequency
;
;         datain, 'cgis_fits:fip_sky_lhs.fits',dsfield='SIGNAL',subscr='5:40',$
;                 data=data,frequency=frequency
;
;         datain, '$CGIS_FITS/dirbe_galactic_plane_maps.fits',data=sigma,$
;                 dsfield='photqual',face=0,/dbsig
;#
;COMMON BLOCKS: None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;
;MODIFICATION HISTORY
; SPR 11347 19-Aug-1993   Dalroy Ward,J. Newmark  initial routine
; SPR 11375 13-Oct-1993   Modify for decomposed PhotQual (DIRBE). J. Newmark
; SPR 11759 17-May-1994   Ingest FIRAS PDS's. J. Newmark
; SPR 11905 07-Sep-1994   Ingest DIRBE PDS's. J. Newmark
; SPR 12140 16-Mar-1995   Ingest FIRAS line maps. J. Newmark
;
; 22-Aug-1997 J. Newmark - support for DIRBE DCAF, DSZA, CIO files
; 19-Nov-1997 J. Newmark - ingest FINAL FIRAS products
;
;.TITLE
;Routine DATAIN
;-

	pro datain, dsname, format=format, intype=intype, face=face, $
                    dsfield=dsfield, subscr=subscr, data=data, $
		    weights=weights, frequency=frequency, coord=coord, $
		    instrume=instrume, badpixval=badpixval, res=res, $
		    units=units, wunits=wunits, rms=rms, indfile=indfile,$
                    dbsig=dbsig,dbqual=dbqual,simple=simple, $
                    nu_zero=nu_zero,num_freq=num_freq,delta_nu=delta_nu,$
                    pixel=pixel,startpix=startpix,stoppix=stoppix

	if n_elements(dsname) eq 0 then begin
	   print, 'DATAIN: dsname=filename is a required field'
	   goto, abort
	endif

	if dsname eq '?' then goto, help

	if n_elements(dsfield) eq 0 then dsfield = '?'

	pos = strpos(dsname,'*')
	if pos ge 0 then begin
	   print, 'DATAIN: wildcards not allowed in file specification'
	   goto, abort
	endif

	files = findfile(dsname)
	if files(0) eq '' then begin
	   print, 'DATAIN: '+dsname+' not found'
	   goto, abort
	endif

	if n_elements(format) gt 0 then begin
	   formats = ['FITS','CISS','CSM']
	   pos = where(STRUPCASE(format) eq formats)
	   if pos(0) eq -1 then begin
	      print, 'DATAIN: format is a required field'
	      goto, abort
	   endif
	   inform = formats(pos(0))
	endif else inform = 'FITS'

	if inform eq 'FITS' then begin
	   result = is_fits(dsname,fits_exten)
	   if result(0) ne 1 then begin
	      print, 'DATAIN: you have specified FITS format but this is ',$
                     'not a fits file'
	      goto, abort
	   endif

           if fits_exten ne 'BINTABLE' then dsfield='data'
           if keyword_set(simple) then goto, simple
	   header = headfits(dsname)

          instrume = sxpar(header, 'INSTRUME')
          if strtrim(instrume) eq 'DIRBE' then begin
             product_type = sxpar(header, 'PRODUCT')
             if strmid(product_type,0,4) eq 'DCAF' or $
                strmid(product_type,0,4) eq 'DSZA' then goto, rddcaf
             if strmid(product_type,0,3) eq 'CIO' then goto, rdcio
          endif 
          if strtrim(instrume) eq 'FIRAS' then begin
             product_type = sxpar(header, 'PRODUCT')
;trap for non-skymap data sets. Perform list read only on these.
             if strtrim(product_type) eq 'HLINE-PR' or $
               strtrim(product_type) eq 'LLINE-PR' or $
               strtrim(product_type) eq 'CAL-ERR' or $
               strtrim(product_type) eq 'CCAL-ERR' or $
               strtrim(product_type) eq 'COV-MAT' or $
               strtrim(product_type) eq 'C-VEC' or $
               strtrim(product_type) eq 'A-VEC' or $
               strtrim(product_type) eq 'PHYS-ST' or $
               strtrim(product_type) eq 'ORTHO-ST' or $
               strtrim(product_type) eq 'MOD-ERRS' or $
               strtrim(product_type) eq 'DS-QMAT' or $
               strtrim(product_type) eq 'DS-CSPEC' or $
               strtrim(product_type) eq 'DS-SCMAT' or $
               strtrim(product_type) eq 'CAL-MOD' then goto, rdfiras
          endif

	   weight = sxpar(header, 'WEIGHTS')
	   if !err lt 0 then weight = ''

	   if weight ne '' then begin
	      print, ' '
	      print, 'The file that you have selected has a weighting array'
	      question = 'Would you like to read in the weighting array? '
	      wans = one_line_menu([question,'yes','no'],init=1)
	      if wans eq 1 then weight = 'yes'
	   endif
	endif

	if n_elements(intype) eq 0 then intype = ''

	if n_elements(face) eq 0 then begin
	   face = 7
        endif else begin
	   if face lt 0 or face gt 5 then begin
	      print, 'DATAIN: face must be in the range 0-5 or 7'
	      goto, abort
	   endif
	endelse

        IF (inform EQ 'CISS') THEN dsfield='data'
        dsfield=STRUPCASE(dsfield)
	while strlen(dsfield) lt 8 do dsfield = dsfield + ' '

	if n_elements(subscr) eq 0 then subscr = '' else begin
           Dims   = StrParse( subscr, Delim=':' )
           N_Dims = N_Elements( Dims )
           If (N_Dims gt 2) then begin
             Print, ' '
             Print, 'Error: format of ranges not correct.' 
             Print, ' '
             Print, 'Press any key to continue...'
             Ans = Get_Kbrd( 1 )
             Goto, abort
           EndIf
           For j = 0, N_Dims-1 do begin
             If (Fix(Dims(j)) le  0) then begin
                Print, ' '
                Print, 'Error: subscripts must be greater than 0.'
                Print, ' '
                Print, 'Press any key to continue...'
                Ans = Get_Kbrd( 1 )
                Goto, abort
             EndIf 
           endfor
          endelse
        if keyword_set(dbsig) then dbsig=1 else dbsig=0
        if keyword_set(dbqual) then dbqual=1 else dbqual=0

	output_format = 'UIDL'
	outname = ''

	dconvert, inform, intype, dsname, fits_exten, face, dsfield, subscr, $
                  weight, output_format, outname, data=data, $
                  weights=weights, frequency=frequency, instrume=instrume, $
		  badpixval=badpixval, res=res, units=units, wunits=wunits, $
	          rms=rms, coord=coord, nu_zero=nu_zero, num_freq=num_freq, $
                  delta_nu=delta_nu,dbsig=dbsig,dbqual=dbqual

	return
simple:
        data=readfits(dsname,header)
        question='Do you wish to examine the header? '
	wans = one_line_menu([question,'yes','no'],init=1)
	if wans eq 1 then begin
	    uscroll, header, title='FITS header for'+dsname
        endif
        return
rdfiras:
        firasmod,dsname,dsfield,data=data,units=units
        return
rddcaf:
        rd_dcaf,dsname,product_type,dsfield,data=data,pixel=pixel,$
            startpix=startpix,stoppix=stoppix,indfile=indfile
        return
rdcio:
        if n_elements(stoppix) ne 0 then $
            message,/info,'RD_CIO reads only single pixels, stoppix ignored.'
        rd_cio,dsname,product_type,dsfield,data=data,pixel=pixel,$
            startpix=startpix,indfile=indfile
        return
abort:
	print, 'Usage: datain, dsname, format=format, intype=intype, face=face,$'
        print, '               dsfield=dsfield, subscr=subscr, data=data,$'
	print, '               weights=weights, frequency=frequency, res=res,$'
	print, '               instrume=instrume, badpixval=badpixval, $'
	print, '               units=units, wunits=wunits, coord=coord,$'
	print, '               num_freq=num_freq, delta_nu=delta_nu, rms=rms,$'
        print, '               nu_zero=nu_zero,/dbsig,/dbqual,/simple'

	return

help:
	return
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



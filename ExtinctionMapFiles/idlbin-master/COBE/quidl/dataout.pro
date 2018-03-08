;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    DATAOUT is a routine which writes out FITS files or CISS files
;
;DESCRIPTION:
;	DATAOUT writes out a FITS file or CISS file from UIDL.
;
;WARNINGS:
;       IDL Save Set output is equivalent to using the IDL command SAVE,
;       except a full sky T will be stored in sixpack format.
;
;CALLING SEQUENCE:
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    param		I/O   	type		description
;    datafile            I      string		full name of file for output
;    weightf=weightf     I      string          full name of file for weights
;    format=format       I      string		FITS or CISS (default to FITS)
;    face=face           I      integer         face # 0-5 or 7 for full image
;    data=data           I      float           array for data input
;    weights=weights     I      float           array for weights
;    frequency=frequency I      float           array for frequency info
;    instrume=instrume   I      string          the instrument the data is for
;    badpixval=badpixval I      float		"no data" indicator
;    res=res             I      integer		SKYMAP resolution
;    units=units         I      string		units for data values
;    wunits=wunits       I      string		units for weights
;    title=title	 I      string		title for data
;    wtitle=wtitle       I      string		title for weights
;    rms=rms             I      float           DMR specific
;    coord=coord         I      string          DMR specific
;    nu_zero=nu_zero     I      float           FIRAS specific
;    delta_nu=delta_nu   I      float           FIRAS specific
;    num_freq=num_freq   I      float           FIRAS specific
;
;    DMR specific keywords
;    ---------------------
;    rms, coord
;
;    FIRAS specific keywords
;    -----------------------
;    nu_zero, delta_nu, num_freq
;
;    DMR and FIRAS specific keywords
;    -------------------------------
;    weights, wunits
;
;    DIRBE and FIRAS specific keywords
;    ---------------------------------
;    frequency
;
;EXAMPLE: 
;         dataout, 'dmr_skymap_31a.fits',data=data,instrume=instrume
;
;COMMON BLOCKS: None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;
;MODIFICATION HISTORY
; SPR 11347 19-Aug-1993   Dalroy Ward, J. Newmark initial routine
;
;
;.TITLE
;DATAIN
;-

	pro dataout, dsname, format=format, face=face, $
                    data=data, weights=weights, frequency=frequency, $
 		    instrume=instrume, badpixval=badpixval, coord=coord, $
                    weightf=weightf, title=title, wtitle=wtitle,$
		    res=res, units=units, wunits=wunits, rms=rms, $
                    num_freq=num_freq,delta_nu=delta_nu,nu_zero=nu_zero

	if n_elements(dsname) eq 0 then begin
	   print, 'DATAOUT: datafile is a required field'
	   goto, abort
	endif

	if dsname eq '?' then goto, help

	if n_elements(format) gt 0 then begin
	   formats = ['FITS','CISS']
	   pos = where(STRUPCASE(format) eq formats)
	   if pos(0) eq -1 then begin
	      print, 'DATAOUT: format is a required field'
	      goto, abort
	   endif
	   output_format = formats(pos(0))
	endif else output_format = 'FITS'

	if n_elements(face) eq 0 then begin
           sz=SIZE(data)
	   if sz(1) eq sz(2) then read,$
             'The data appears to be a single face. Please enter face: ',$
             face ELSE face = 7
        endif else begin
	   if face lt 0 or face gt 5 then begin
	      print, 'DATAIN: face must be in the range 0-5 or 7'
	      goto, abort
	   endif
	endelse

	; setup some defaults if the user hasn't given us any information
	if n_elements(badpixval) eq 0 then badpixval = 0
	if n_elements(instrume) eq 0 then instrume = 'UNKNOWN'
	if n_elements(units) eq 0 then units = 'UNKNOWN'
	if n_elements(wunits) eq 0 then wunits = 'UNKNOWN'
	if n_elements(frequency) eq 0 then frequency = 0.
	if n_elements(weightf) eq 0 then weightf = ''
	if n_elements(intype) eq 0 then intype = ''

        inform='UIDL'
	dummy = ''

	dconvert, inform, intype, dummy, fits_exten, face, dummy, subscr, $
                  weightf, output_format, dsname, $
                  data=data, weights=weights, frequency=frequency,$
		  instrume=instrume, badpixval=badpixval, res=res, $
		  units=units, wunits=wunits, rms=rms, coord=coord, $
                  num_freq=num_freq, nu_zero=nu_zero, delta_nu=delta_nu

	return

abort:
	print, 'Usage: dataout, datafile,format=format,intype=intype,$'
        print, '               dsfield=dsfield, data=data,$'
	print, '               weights=weights, frequency=frequency, $
	print, '               instrume=instrume, badpixval=badpixval, $'
	print, '               units=units, wunits=wunits, rms=rms,face=face,$'
        print, '               num_freq=num_freq, coord=coord, $ '
        print, '               nu_zero=nu_zero, delta_nu=delta_nu '

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



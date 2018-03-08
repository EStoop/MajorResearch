Pro csmhdr, 		$
	dsname, $
	name,         $
       	description,  $
       	version,      $
       	units,         $
       	channels,     $
       	index_levels, $
	bandfreq,     $
	bandwidth,    $
	bandunits,    $
       	bandnull

;+
; NAME:
;	CSMHDR
;
; PURPOSE:
;       This procedure reads the header file for the 
;	specified data set 
;	
; CATEGORY:
;	User interface, Menu.
;
; CALLING SEQUENCE:
;	csmhdr, dsname, name, ... 
;	
; INPUTS:
;	dsname	file name of the data set
;
; OUTPUTS:
;	name 		full name of the data base
;	description	text description 
;	version		text version number 
;	units		units descriptor for data, see band units if >1
;	channels	number of channels included, size of S-array
;	index_levels    indicates resolution of map, size of pixel
;	bandfreq	frequency of each band or channel
;	bandwidth	width of each band or channel
;	bandunits	units descriptor for each band or channel
;	bandnull 	signal missing flag
;	
; SIDE EFFECTS:
;	The header data file is opened and read through the data server.
;
; COMMON BLOCKS:
;	None.
;
; RESTRICTIONS:
;	Requires the Data Server software. 
;
; SUBROUTINES/FUNCTIONS CALLED:
;	Read_RMS
;
; MODIFICATION HISTORY:
;       Creation:  Pete Kryszak, GSC, March 1992.
;       SPR 11003 2 Jun 93   Change !version.os to !cgis_os. J Newmark
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;-
;
message,  "this is csmhdr...", /informational
message,  "the name of the data set "+dsname, /informational
;
IF ( !cgis_os EQ 'vms' ) THEN BEGIN; it's vms
   pref = ''
   sep = ':'
   eod = ']'   ;for later
ENDIF
IF (!cgis_os EQ 'windows') THEN BEGIN ; pc idl
    pref = ' '
    sep = '\'
    eod = '\'
ENDIF 
IF (!cgis_os EQ 'unix') THEN BEGIN ;unix compatible
    pref = '$'
    sep = '/'
    eod = '/'
ENDIF


i = STRPOS( dsname, '.')
j = STRPOS( dsname, sep )
k = i - j

header_name = 'ADB_HDR' + STRMID( dsname, j, k) + $
		'*.HDR*;0' 

;print, 'csmhdr.pro ---> ', header_name

header_name = findfile( header_name, count=hfcount)
if hfcount Eq 0 Then Begin
	message, /cont, 'Header file not found.'
	return
	EndIf
header_name = header_name(0)

;print, header_name
;print, i, j, k
;
message,  "the name of the header is: "+header_name, /informational
;
;
fields = 'name,description,version,units,index_levels'
Set_RDF, 'csdr$rdf:xac_header.rdf'
print, 'reading header information'
status = $
 Read_RMS(header_name, fields, name, description, version, units, index_levels)
If status ne 1 then begin
	message,/continue, 'Bad call to Read_RMS.'
	!err = -1
	!error = -1
	return
	EndIf
;
fields = 'channels,band.freq,band.width,band.units,band.null'
Set_RDF, 'csdr$rdf:xac_header.rdf'
status = $
 Read_RMS(header_name, fields, channels, $
		bandfreq, bandwidth, bandunits, bandnull )
If status ne 1 then begin
	message,/continue, 'Bad second call to Read_RMS.'
	!err = -1
	!error = -1
	return
	EndIf
;

;help,index_levels

Return

End
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



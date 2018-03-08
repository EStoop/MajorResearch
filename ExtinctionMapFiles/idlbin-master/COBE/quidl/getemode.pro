PRO getemode, entry_mode,noshow,coord,infmt,re_proj
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    getemode gets the entry mode
;
;DESCRIPTION:
;    This IDL procedure querys the user for entry mode (cursor or
;    long/lat) and if the latter then for coordinate system.
;
;CALLING SEQUENCE:
;    getemode, entry_mode,noshow,coord,infmt,re_proj
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    entry_mode  O     int                 entry mode (1-cur,2-lon/lat)
;    noshow      I     int                 no display flag (1-on)
;    coord       O     str                 Coordinate system
;    infmt       O     str                 Format ('L' if ll entry)
;    re_proj     O     int                 reprojection flag
;
;WARNINGS:
;
;    None.
;
;EXAMPLES: 
;
;    Not user routine.
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    None
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Jun 92
;
; SER 9955    Improved Documentatin/Banner
; 9-SEP-92    J.M. Gales
;
;-
coord_type = 'EGQ'

IF (noshow EQ 1) THEN entry_mode = 2 ELSE BEGIN

	entry_mode = umenu(['Select Entry Mode', $
			    'Cursor', $
			    'Longitude/Latitude', $
			    'HELP'],title=0,init=1)

ENDELSE

CASE entry_mode OF

1 :	IF (re_proj EQ 0) THEN coord = 'E'
2 :	infmt = 'L'

3 :	BEGIN
	;HELP
	END

ENDCASE

; If not cursor entry and quad cube then ask for coord sys
; --------------------------------------------------------
IF ((entry_mode NE 1) AND (re_proj EQ 0)) THEN BEGIN

	entry = umenu(['Coordinate System','Ecliptic', $
                       'Galactic', 'Equatorial'], $
                        title=0,init=1)

	coord = STRMID(coord_type,entry-1,1)

ENDIF

RETURN
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



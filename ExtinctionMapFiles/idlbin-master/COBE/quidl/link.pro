;*******************************************************************************
;
;	This procedure linkimages user-written application utilities.
;
;	AUTHOR:	Song Yom (HUGHES STX) (derived from LINKIMAGES.PRO)
;	DATE:	03/92 (SER# 9525)
;
;	MODIFICATION HISTORY:
;	DATE	SPR/SER		PROGRAMMER	DESCRIPTION
;	----	-------		----------	-----------
;	05/92   9678		Pete Kryszak	read_pixel and read_pixstr
;                                               renamed to DIRBE_...
;       09/92   9823		Pete Kryszak	removed sum.exe
;       05/94   11762                           added btb_moon_angles
;
;*******************************************************************************
SETLOG,'DIRBE_PIXEL_EXE','CGIS$IDL:DIRBE_PIXEL.EXE'
SETLOG,'FMEDIAN_EXE','CGIS$IDL:FMEDIAN.EXE'
SETLOG,'DIRBE_READ_PIXEL_EXE','CGIS$IDL:DIRBE_READ_PIXEL.EXE'
SETLOG,'DIRBE_READ_PIXSTR_EXE','CGIS$IDL:DIRBE_READ_PIXSTR.EXE'
SETLOG,'SETUP_PLT_EXE','CGIS$IDL:SETUP_PLT.EXE'
SETLOG,'BTB_MOON_ANGLES_EXE','CGIS$IDL:BTB_MOON_ANGLES.EXE'
LINKIMAGE,'DIRBE_PIXEL','DIRBE_PIXEL_EXE',1,'DIRBE_PIXEL_WRAPPER'
LINKIMAGE,'DIRBE_READ_PIXEL','DIRBE_READ_PIXEL_EXE',1,'DIRBE_READ_PIXEL_WRAPPER'
LINKIMAGE,'DIRBE_READ_PIXSTR','DIRBE_READ_PIXSTR_EXE',1,'DIRBE_READ_PIXSTR_WRAPPER'
LINKIMAGE,'FMEDIAN','FMEDIAN_EXE',0,'FMEDIAN_WRAPPER'
LINKIMAGE,'SETUP_PLT','SETUP_PLT_EXE',1,'SETUP_PLT_WRAPPER'
LINKIMAGE,'BTB_MOON_ANGLES','BTB_MOON_ANGLES_EXE',1,'BTB_MOON_ANGLES_WRAPPER'
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



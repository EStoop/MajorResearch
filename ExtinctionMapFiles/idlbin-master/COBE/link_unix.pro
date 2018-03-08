;*******************************************************************************
;
;	This procedure linkimages user-written application utilities.
;
;       Written: 11-AUG-1993  J. Newmark
;
;*******************************************************************************
;
; edit PATH to reflect the actual location, e.g. /usr/local/rsi/idl_5/cgis
;
 LINKIMAGE,'cconv','/PATH/cconv.so',1,'cconv' 
 LINKIMAGE,'projtran','/PATH/projtran.so',1,'projtran' 
 LINKIMAGE,'adt2zulu','/PATH/adt2zulu.so',0,'adt2zulu' 
 LINKIMAGE,'rastr','/PATH/rastr.so',0,'rastr' 
 LINKIMAGE,'zulu2adt','/PATH/zulu2adt.so',0,'zulu2adt' 
 LINKIMAGE,'cpusec','/PATH/cpusec.so',1,'cpusec'
 LINKIMAGE,'inpoly','/PATH/inpoly.so',1,'inpoly'
 LINKIMAGE,'pixavg','/PATH/pixavg.so',0,'pixavg'
 LINKIMAGE,'dblsvd','/PATH/dblsvd.so',0,'dblsvd'
 LINKIMAGE,'qrdcmp','/PATH/qrdcmp.so',0,'qrdcmp'

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



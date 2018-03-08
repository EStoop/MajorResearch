;*******************************************************************************
;
;	This procedure linkimages user-written application utilities.
;
;       Written: 11-AUG-1993  J. Newmark
;
;*******************************************************************************
;
; this assumes C executables are in the logical directory #CGIS_IDL#
; edit to reflect actual location
;

 SETLOG,'cconv_exe','#CGIS_IDL#:cconv.exe' 
 SETLOG,'projtran_exe','#CGIS_IDL#:projtran.exe' 
 SETLOG,'adt2zulu_exe','#CGIS_IDL#:adt2zulu.exe' 
 SETLOG,'rastr_exe','#CGIS_IDL#:rastr.exe' 
 SETLOG,'zulu2adt_exe','#CGIS_IDL#:zulu2adt.exe' 
 SETLOG,'cpusec_exe','#CGIS_IDL#:cpusec.exe' 
 SETLOG,'inpoly_exe','#CGIS_IDL#:inpoly.exe'
 SETLOG,'pixavg_exe','#CGIS_IDL#:pixavg.exe'
 SETLOG,'dblsvd_exe','#CGIS_IDL#:dblsvd.exe'
 SETLOG,'qrdcmp_exe','#CGIS_IDL#:qrdcmp.exe'
 LINKIMAGE,'cconv','cconv_exe',1 
 LINKIMAGE,'cpusec','cpusec_exe',1 
 LINKIMAGE,'projtran','projtran_exe',1 
 LINKIMAGE,'adt2zulu','adt2zulu_exe',0 
 LINKIMAGE,'rastr','rastr_exe',0 
 LINKIMAGE,'zulu2adt','zulu2adt_exe',0 
 LINKIMAGE,'inpoly','inpoly_exe',1
 LINKIMAGE,'pixavg','pixavg_exe',0
 LINKIMAGE,'dblsvd','dblsvd_exe',0
 LINKIMAGE,'qrdcmp','qrdcmp_exe',0

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



pro fdecomp,filename,disk,dir,name,qual,version
;+ NAME/ONE LINE DESCRIPTION OF ROUTINE:
;      FDECOMP decomposes a file name into its components.
;
;  DESCRIPTION:
;      IDL procedure to decompose a file name into its
;      constituent components (disk, directory, file name,
;      file name extension, and version number).
;
;  CALLING SEQENCE:
;      FDECOMP,FILENAME,DISK,DIR,NAME,QUAL,VERSION
;
;  ARGUMENTS (I = input, O = output, [] = optional):
;      FILENAME       I      str      File name to be decomposed
;      DISK           O      str      Disk name
;      DIR            O      str      Directory name
;      NAME           O      str      File name
;      QUAL           O      str      File name extension
;      VERSION        O      str      Version number
;
;  WARNINGS:
;      None
;
;  EXAMPLE:
;      To decompose the file name
;      FIRCOADD:[FIRAS.SKYMAP]FCF_SKY_LLSS.ED_8934301_9026410;2:
;
;      filename =
;            'FIRCOADD:[FIRAS.SKYMAP]FCF_SKY_LLSS.ED_8934301_9026410;2'
;      fdecomp,filename,disk,dir,name,qual,ver
;
;      FDECOMP returns these strings for the file name components:
;            disk = 'fircoadd'
;            dir  = '[firas.skymap]'
;            name = 'fcf_sky_llss'
;            qual = 'ed_8934301_9026410'
;            ver  = '2'
;#
;  COMMON BLOCKS:
;      None
;
;  PROCEDURE (AND OTHER PROGRAMMING NOTES):
;      The input and output variables to this procedure are all string
;      arrays.
;
;  PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;      GETTOK
;
;  MODIFICATION HISTORY:
;      version 1  D. Lindler  Oct 1986
;
; SPR 9616
;.TITLE
;Routine FDECOMP
;-
;
st=filename
;
; get disk
;
if strpos(st,':') ge 0 then disk=gettok(st,':')+':' else disk=''
;
; get dir
;
if strpos(st,']') ge 0 then dir=gettok(st,']')+']' else dir=''
;
; get name
;
name=gettok(st,'.')
;
; get qualifier
;
qual=gettok(st,';')
;
; get version
;
version=st
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



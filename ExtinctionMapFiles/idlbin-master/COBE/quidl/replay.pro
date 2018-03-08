;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     REPLAY redisplays an image previously stored using 'capture.pro'
;
;DESCRIPTION:
;     This procedure reads an image from an unformatted file on disk
;     previously saved by 'capture', and redisplays it on an IDL 
;     display window. The display window is automatically sized to
;     the restored image array size.
;
;CALLING SEQUENCE:
;     REPLAY, file, [vimage]
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     FILE       I   str       Name of unformatted file where 'capture' 
;                              saved the IDL image  
;     VIMAGE    [I]  bytarr    Optional name of IDL array into which 
;                              saved file is to be restored. Default 
;                              name is 'image'
;
;WARNINGS:
;     1.  This procedure can only be called from a Workstation 
;         (X-Windows) session.
;     2.  The specified disk file produced by 'Capture' must exist.
;     
;EXAMPLE:
;     To redisplay an image stored in disk file 'output.dat' (as
;     created by capture), and save the image in an array called data. 
;
;     REPLAY, 'output.dat', data
;#
; 
;COMMON BLOCKS:
;     None  
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES) :
;     The image and image dimensions are read from the disk file. The 
;     display window is sized to the dimensions of the image restored.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     Capture
; 
; MODIFICATION HISTORY:
;   Written by Mila Mitra, General Sciences Corporation, March 1991
;   SER 9616
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;.TITLE
;Routine REPLAY
;-
pro replay,file,vimage
; 
openr,1,file+"/unf"
s=lonarr(2)
readu,1,s
if n_elements(vimage) eq 0 then vimage=vimage
vimage=bytarr(s(0),s(1))
readu,1,vimage
close,1  
;
;SETS WINDOW AND DISPLAYS
window,xs=s(0),ys=s(1)
tv,vimage
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



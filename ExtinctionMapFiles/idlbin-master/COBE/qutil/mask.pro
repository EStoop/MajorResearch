;+
; NAME: Mask
;
; PURPOSE: To mask (set equal to 0) selected elements of input array.
;          This is used when input list of elements is a string which
;          can contain spaces, commas or specify a range (x:y) of data.
;
; CATEGORY: Utility
;
; CALLING SEQUENCE: mask,input,array,outarray
;
; INPUTS:
;          input        a string containing a list of numbers. The string
;                       can contain various delimeters such as ' ',',' or
;                       ':'. Use of the latter is inferred as requesting
;                       a range of numbers to be masked.
;
;          array        the input array to be masked
;
; OUTPUTS:
;	   outarray     the masked array
;#
; WARNINGS: None  
;
; COMMON BLOCKS: None
;
; INCLUDE FILES: None
;
; RESTRICTIONS: None
;
; SUBROUTINES CALLED: get_numbers
;
; MODIFICATION HISTORY: Created J. Newmark May 27, 1993
;---------------------------------------------------------------------
;
pro mask,input,array,outarray
len=strlen(input)
outarray=-1
temp=strarr(len)
FOR i=0,len-1 DO temp(i)=strmid(input,i,1)
rpos = where(temp eq ':',cnt)
IF (cnt GT 1) THEN FOR i=0,cnt-2 DO $
   IF (rpos(i)+1 EQ rpos(i+1)) THEN BEGIN
     PRINT,'Entered masked element array with 2 sequential colons. ' + $
           'Please re-enter the array.'
     RETURN
   ENDIF
outarray=array
IF cnt NE 0 THEN BEGIN
   start=0
   total_num=intarr(cnt+1)
   list=intarr(cnt+1)
   srange=intarr(cnt+1)
   erange=intarr(cnt+1)
   FOR i=0,cnt DO BEGIN
     IF (i eq cnt) THEN send=len ELSE send=rpos(i)
     section=STRMID(input,start,send-start)
     get_numbers,section,num,c
     srange(i)=num(c-1)
     erange(i)=num(0)
     outarray(num-1)=0
     IF (i EQ cnt) THEN goto,next ELSE start=rpos(i)+1
   ENDFOR
   next:
   FOR i=0,cnt-1 DO BEGIN
     start=srange(i)+1
     stop=erange(i+1)-1
     outarray(start-1:stop-1)=0
   ENDFOR
ENDIF ELSE BEGIN
  get_numbers,input,num,c
  outarray(num-1)=0
ENDELSE
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



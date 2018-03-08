pro dberror,photqual,sigma,qual_flag,photomet=photomet,badval=badval
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:    
;     DBERROR separates the SIGMA and QUAL_FLAG from PhotQual.
;
;DESCRIPTION: 
;     This routines takes as input the PHOT_QUAL field of the
;     DIRBE GPM FITS file and extracts the data quality flag and
;     sigma. The input arrays should be a 1,2 or 3 dimensional array with
;     the last dimension as wavelength (10 bands). The map must be
;     either in skycube/face or sixpack format. The first (least 
;     significant) byte of the PHOT_QUAL field is ocmposed of a 
;     series of bit flags where a total value of zero indicates the 
;     standard fitting method was used. A non-zero value indicates a 
;     special condition. 
;     The second (most significant) byte gives the standard deviation. 
;     The first 8 bands are represented as RELATIVE SIGMAS. If the PHOTOMET
;     values are input as well than the absolute error is returned.
;     The last two bands (140, 240 um) the standard deviations
;     themselves are given.
;
;CALLING SEQUENCE: 
;     DBERROR,photqual,sigma,qual_flag [,photomet=photomet]
;
;ARGUMENTS (I = input, O = output, [] = optional): 
;     PHOTQUAL   I   intarr  PhotQual field from DIRBE image.
;     PHOTOMET  [I]  fltarr  Photomet field from DIRBE image.
;     BADVAL    [I]  int     Bad Pixel value (def=0).
;     QUAL_FLAG  O   intarr  Data Quality flag
;     SIGMA      O   fltarr  Standard deviations
;
;WARNINGS:  none.
;
;EXAMPLE: 
;     1. To generate the Qual_flag and sigma from the PhotQual using
;        the photomet to produce the absolute error for bands 1-8:
;
;        dberror,photqual,sigma,qual_flag,photomet=photomet
;
;#
;COMMON BLOCKS: none.
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES)
;         Separates the most significant and least significant bytes of
;         the DIRBE PhotQual field.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     none
;     
;MODIFICATION HISTORY:
;     SPR 11375 Written by J. Newmark, 12-Oct-1993.
;
;.TITLE
;Routine DBERROR
;-

relflag=0
IF (KEYWORD_SET(badval) EQ 0) THEN badval=0.
IF (N_ELEMENTS(photomet) NE 0) THEN $
   IF (N_ELEMENTS(photomet) ne N_ELEMENTS(photqual)) THEN BEGIN
     print, 'The number of elements on PHOTOMET and PHOTQUAL must be equal'
     return
   ENDIF ELSE relflag=1
;
; Extract the quality flag (least significant byte).
qual_flag = photqual - photqual/256*256
;
; Extract sigma (most significant byte).
data_list=LONG(photqual)
bad=WHERE(data_list EQ badval)
blow=where(data_list LT 0,cnt)
if cnt gt 0 then data_list(blow)=data_list(blow)+65536
sz=size(data_list)

IF (sz(0) EQ 1) THEN BEGIN
   dbsig=fltarr(sz(1))
   FOR ib=0,9 DO BEGIN
    itemp=data_list(ib)/256
    n=3
    IF (ib GT 7) THEN n=1
    CASE itemp OF
       0: dbsig(ib)=10.0^(-n)
       255: dbsig(ib)=10.0^(-n+4)
       else: dbsig(ib)=10.0^((itemp - 0.5)/63.5 - n)
    ENDCASE
    IF (ib LE 7) AND (relflag EQ 1) THEN dbsig(ib) = $
                           dbsig(ib)*ABS(photomet(ib))
   ENDFOR
ENDIF 
IF (sz(0) EQ 2) THEN BEGIN
   dbsig=fltarr(sz(1),sz(2))
   FOR ib=0,9 DO BEGIN
    itemp=data_list(*,ib)/256
    n=3
    IF (ib GT 7) THEN n=1
    dim=SIZE(itemp)
    FOR k=0l,dim(1)-1 DO $
    CASE itemp(k) OF
       0: dbsig(k,ib)=10.0^(-n)
       255: dbsig(k,ib)=10.0^(-n+4)
       else: dbsig(k,ib)=10.0^((itemp(k) - 0.5)/63.5 - n)
    ENDCASE
    IF (ib LE 7) AND (relflag EQ 1) THEN dbsig(*,ib) = $
                           dbsig(*,ib)*ABS(photomet(*,ib))
   ENDFOR
ENDIF 
IF (sz(0) EQ 3) THEN BEGIN
  dbsig=fltarr(sz(1),sz(2),sz(3))
  FOR ib=0,9 DO BEGIN
    itemp = data_list(*,*,ib)/256
    n = 3
    IF (ib GT 7) THEN n = 1
    dim=SIZE(itemp)
    FOR k=0l,dim(1)-1 DO FOR j=0l,dim(2)-1 DO $ 
        CASE itemp(k,j) OF
           0: dbsig(k,j,ib) = 10.0^(-n) 
           255: dbsig(k,j,ib) = 10.0^(-n+4) 
           else: dbsig(k,j,ib) = $
                10.0^((itemp(k,j) - 0.5)/63.5 - n)
        ENDCASE
    IF (ib LE 7) AND (relflag EQ 1) THEN dbsig(*,*,ib) = $
                           dbsig(*,*,ib)*ABS(photomet(*,*,ib))
  ENDFOR
ENDIF

sigma=dbsig
IF (bad NE -1) THEN sigma(bad)=badval
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



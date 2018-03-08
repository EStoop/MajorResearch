function decompress_mf,indata
;+NAME/ONE LINE DESCRIPTION OF ROUTINE: 
;     DECOMPRESS_MF decompresses dirbe data read from BLI_TOIRSCLD.
;
; DESCRIPTION:
;     IDL function to decompress the science data word which was 
;     compressed by BTL_COMPRESS_SDM_MF.  The algorithm is based on 
;     BTL_DECOMPRESS_SDM_MF.
;
; CALLING SEQUENCE:
;     DATA = DECOMPRESS_MF(indata)
;
; ARGUMENTS:
;     DATA     O   fltarr    An array of floating point (R*4) numbers 
;                            with same dimension as INDATA.
;     INDATA   I   intarr    An array of compressed integers (I*2) of 
;                            arbitrary dimension.
;
; WARNINGS:
;     None.
;
; EXAMPLE:
;     To read DIRBE data from BLI for the 1st hour of 01-January-1990, 
;     and convert the data from a compressed I*2 array to a decompressed 
;     R*4 array:
;
;       data = READ_ARCV('bli_toirscld.dadrbsci2','9000100','9000101')
;       data = DECOMPRESS_MF(data)
;
; COMMON BLOCKS:
;     None.
;
; PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     The algorithm is based on BTL_DECOMPRESS_SDM_MF.FOR.
;
; MODIFICATION HISTORY:
;        Written By: BA Franz, JP Weiland, March 1992 
;        CCR 585
;.TITLE
;Routine DECOMPRESS_MF
;-

outdata = float(indata)

;
; Set minimum data size from BTL_BLI_COMPRESSION_PAR.TXT
;
min_data_size          = 0.5/16./27.12

;
; Divide data into +/- small numbers, large positive and large negative
; numbers, and sentinels.
;
small = where( abs(indata) lt 2048)
lpos  = where((abs(indata) ge 2048) and (indata gt 0))
lneg  = where((abs(indata) ge 2048) and (indata lt 0) and (indata ge -28358))
sent1 = where((indata lt -28358) and (indata gt -32768))
sent2 = where( indata eq -32768)

;
; Decompress small numbers
;
if (small(0) ne -1) then outdata(small) = indata(small)*min_data_size
	
;
; Decompress large positive numbers
;
if (lpos(0) ne -1) then begin
    w11_14 = (indata(lpos) and 30720)/2048
    w0_10  = (indata(lpos) and 2047)
    outdata(lpos) = float(w0_10) * 2.0^w11_14 * min_data_size
endif

;
; Decompress large negative numbers
;
if (lneg(0) ne -1) then begin
    stuff  = -1*indata(lneg)
    w11_14 = (stuff and 30720)/2048
    w0_10  = (stuff and 2047)
    outdata(lneg) = -1.0 * float(w0_10) * 2.0^w11_14 * min_data_size
endif

;
; Form sentinels
;
if (sent1(0) ne -1) then begin
    outdata(sent1) = -16375.0 + (float(indata(sent1)) + 28360.0)
endif

if (sent2(0) ne -1) then begin
    outdata(sent2) = -20783.0
endif

return,outdata
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



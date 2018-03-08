FUNCTION validnum,input
;+
; NAME:
;    VALIDNUM
; PURPOSE:
;    To check the validity of the user input.
; CALLING SEQUENCES:
;	1)  result = validnum('input')
; INPUTS
;
;	input  - character string
; OUTPUTS:
;	result which will have a value 1 if the input is valid
;       and 0 otherwise.
; HISTORY:
;	Created S.R.K. VIDYA SAGAR (ARC)  APR 1992
;-
;------------------------------------------------------------------------
ON_IOERROR, fail
input = strlowcase(input)
y = DOUBLE(input)
;
; Compute the length of input string
;
inp_len = STRLEN(input)
check = BYTE(STRMID(input, inp_len-1, 1))
check = check(0)
;
; If the last character of the string is not a '.' and not between
; 0-9 then it is an invalid number.
;
IF (check LT 48 OR check GT 57 AND (check NE 46)) THEN RETURN,0
;
Label:
start = 0
n_pluses = 0
n_minuses = 0
n_dots = 0
n_exponentials = 0
minus_sign = BYTE(STRMID(input,0,1))
minus_sign1 = minus_sign(0)
minus_sign = BYTE(STRMID(input,1,1))
minus_sign2 = minus_sign(0)
pos_of_exp = 0
pos_of_period = 0
FOR k = 0, inp_len-1 DO BEGIN
    ik = BYTE(STRMID(input, k, 1))
    ik = ik(0)
    IF (ik eq 101) THEN pos_of_exp = k
    IF (ik EQ 46)  THEN pos_of_period = k
ENDFOR
IF (pos_of_period GT pos_of_exp AND pos_of_exp NE 0) THEN RETURN,0
IF (pos_of_exp GE 1) THEN BEGIN
    exp_val = FIX(STRMID(input,pos_of_exp+1,99))
    IF (exp_val GE 38  OR  exp_val LE -38) THEN RETURN,0
ENDIF
FOR K = 0, inp_len-1  DO BEGIN
    ik = STRMID(input,start,1)
    bik = BYTE(ik)
    bik = bik(0)
    IF ((bik EQ 45 OR bik EQ 43) AND start NE 0) THEN BEGIN
       tik = BYTE(STRMID(input,start-1,1))
       tik = tik(0)
       IF (tik NE 101) THEN RETURN,0
    ENDIF
    start = start + 1
    IF (bik EQ 43 OR bik EQ 45 OR bik EQ 46 $
        OR bik EQ 101 OR (bik GE 48 AND bik LE 57)) THEN BEGIN
    ENDIF ELSE BEGIN
        return,0
    ENDELSE
    IF (bik EQ 43) THEN n_pluses = n_pluses + 1
    IF (bik EQ 45) THEN n_minuses = n_minuses + 1
    IF (bik EQ 46) THEN n_dots = n_dots + 1
    IF (bik EQ 101) THEN n_exponentials = n_exponentials + 1
    IF (n_pluses GE 2 OR n_dots GE 2 OR n_exponentials GE 2) THEN BEGIN
        return,0
    ENDIF
    IF ((minus_sign1 EQ 45 AND minus_sign2 EQ 45) OR $
        (minus_sign1 EQ 45 AND n_minuses GE 3) OR $
        (minus_sign1 NE 45 AND n_minuses GE 2)) THEN BEGIN
        return,0
    ENDIF
ENDFOR    
RETURN,1
fail:
RETURN,0
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



FUNCTION weightedsum,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,$
  coef=coef,badval_in=badval_in,outbadval=outbadval
;+NAME/ONE LINE DESCRIPTION:
;    WEIGHTEDSUM adds maps together (properly handling bad pixels).
;
;DESCRIPTION:
;    WEIGHTEDSUM receives up to ten 2-D arrays as input and returns a 2-D
;    array as output.  Any pixel which was designated as "bad" in any
;    of the input arrays will also be designated as bad in the output
;    array.  Pixels which had good values in all input arrays will have,
;    in the output array, a new value which is either the sum or the
;    weighted sum of the values of the corresponding pixels in the input
;    arrays.  (A weighted sum will be performed if an array of coeffients
;    is passed in via the keyword COEF).
;
;CALLING SEQUENCE:
;    result = weightedsum( x1, x2, ... , in, badval_in=[...], coef=[...],$
;               outbadval=[...])
;
;ARGUMENTS (I=input, O=output, []=optional)
;    x1 to x10  I   2-D arr  flt        An array to be included in the
;                                       summation.
;    badval_in [I]  keyword,1-D arr,flt A 1-D array that contains the bad
;                                       pixel values for each input array.
;                                       The size of the array must be
;                                       equal to the number of input
;                                       images.  All zeros will be
;                                       assumed if not supplied.
;    coef      [I]  keyword,1-D arr,flt A 1-D array that contains the
;                                       coefficients for a weighted sum.
;                                       The size of the array must be
;                                       equal to the number of input
;                                       images plus one.  COEF(0) is
;                                       treated as an additive constant.
;                                       [0.,1.,1.,...,1.] is assumed if
;                                       not supplied.  (The number of 1's
;                                       is equal to N_PARAMS(0).)
;    outbadval [O]  keyword  flt        The bad pixel value for result.
;    result     O   2-D arr  flt        The "sum" or "weighted sum" of
;                                       the supplied images.
;
;EXAMPLE:
;    result = weightedsum( x1, x2, badval_in=[0.,0.], coef=[0.,1.,-1.], $
;               outbad=outbad)
;#
;COMMON BLOCKS:  none.
;
;LIBRARY CALLS:  none.
;
;PROCEDURE:
;    If between one and ten 2-D arrays of the same size have not been
;    passed in, then put out a message and exit.  If BADVAL_IN was not
;    supplied, then assume a value of [0.,0.,...,0.] (where the number
;    of 0's is equal to N_PARAMS(0)), and inform the user that such an
;    assumption is being made.  If a value has been assigned by the
;    caller, and if it is not a 1-D float array with the appropriate
;    numer of elements, then put out a message and exit.  If COEF was
;    not supplied, then assume a value of [0.,1.,...,1.] (where the
;    number of 1's is equal to N_PARAMS(0)), and inform the user that
;    such an assumption is being made.  If a value has been assigned by
;    the caller, and if it is not a 1-D float array with the appropriate
;    numer of elements, then put out a message and exit.  Define good
;    pixels to be those which are good for all input images.  If there
;    are no good pixels under this definition, then return an array with
;    all bad pixels.  Try to catch any potential floating point overflows
;    or underflows by looking at results obtained by summing the logs of
;    the extreme good values from the input images with the logs of the
;    extreme coefficient values.  Put out a message and exit if overflow
;    or underflow is determined to be a possibility.  Create result such
;    that all non-good pixels receive a bad value, and all good pixels
;    receive the weighted sum of the corresponding values from the input
;    images.  If any good pixel has received, as a result of the operation,
;    a value which is equal to the bad pixel value, then change the bad
;    pixel value to be 1 less than the minimum of the good pixel values.
;
;REVISION HISTORY:
;    Originally written by Vidia Sagar, ARC, February 1992.
;    Extensively rewritten by John Ewing, ARC, April 1992.
;    Added handling for graphs, John Ewing, ARC, June 1992.
;
;.TITLE
;Routine WEIGHTEDSUM
;-
;
;  Check supplied images.
;  ----------------------
  IF(N_PARAMS() LT 1) THEN BEGIN
    MESSAGE,'Insufficient parameters. (At least one 2-D array is expected).',/CONT
    RETURN,-1
  ENDIF
  sz = SIZE(x1)
  IF(sz(0) NE 2) THEN BEGIN
    MESSAGE,'Arguments must be 2-dimensional.',/CONT
    RETURN,-1
  ENDIF
  dim1x1 = sz(1)
  dim2x1 = sz(2)
  IF(dim2x1 EQ 2) THEN graphs = 1 ELSE graphs = 0
  FOR i=2,N_PARAMS(0) DO BEGIN
    CASE i OF
       2: sz = SIZE(x2)
       3: sz = SIZE(x3)
       4: sz = SIZE(x4)
       5: sz = SIZE(x5)
       6: sz = SIZE(x6)
       7: sz = SIZE(x7)
       8: sz = SIZE(x8)
       9: sz = SIZE(x9)
      10: sz = SIZE(x10)
    ENDCASE
    IF(sz(0) NE 2) THEN BEGIN
      MESSAGE,'Arguments must be 2-dimensional.',/CONT
      RETURN,-1
    ENDIF
    IF(NOT graphs) THEN BEGIN
      IF((sz(1) NE dim1x1) OR (sz(2) NE dim2x1)) THEN BEGIN
        MESSAGE,'Arrays must be of the same size.',/CONT
        RETURN,-1
      ENDIF
    ENDIF
  ENDFOR
;
;  Check what was passed in for BADVAL_IN.
;  ---------------------------------------
  IF (NOT KEYWORD_SET(badval_in)) THEN BEGIN
    badval_in = FLTARR(N_PARAMS(0))
    IF(NOT graphs) THEN BEGIN
      PRINT,'BADVAL_IN was not set.  Assuming a value of 0 for each argument.'
    ENDIF ELSE BEGIN
      temp = N_PARAMS(0)
      badval_in(0) = MIN(x1(*,1)) - 1.
      IF(temp GE 02) THEN badval_in(1) = MIN(x2(*,1)) - 1.
      IF(temp GE 03) THEN badval_in(2) = MIN(x3(*,1)) - 1.
      IF(temp GE 04) THEN badval_in(3) = MIN(x4(*,1)) - 1.
      IF(temp GE 05) THEN badval_in(4) = MIN(x5(*,1)) - 1.
      IF(temp GE 06) THEN badval_in(5) = MIN(x6(*,1)) - 1.
      IF(temp GE 07) THEN badval_in(6) = MIN(x7(*,1)) - 1.
      IF(temp GE 08) THEN badval_in(7) = MIN(x8(*,1)) - 1.
      IF(temp GE 09) THEN badval_in(8) = MIN(x9(*,1)) - 1.
      IF(temp EQ 10) THEN badval_in(9) = MIN(x10(*,1)) - 1.
    ENDELSE
  ENDIF ELSE BEGIN
    sz = SIZE(badval_in)
    IF((sz(0) NE 1) OR (sz(2) NE 4) OR (sz(1) NE N_PARAMS(0))) THEN BEGIN
      MESSAGE, 'An inappropriate BADVAL_IN was supplied.',/CONTINUE
      RETURN,-1
    ENDIF
  ENDELSE
;
;  Check what was passed in for COEF.
;  ----------------------------------
  IF (NOT KEYWORD_SET(coef)) THEN BEGIN
    PRINT,'COEF was not set.  Assuming a value of 1. for each argument and an'
    PRINT,'additive constant of 0.'
    coef = FLTARR(N_PARAMS(0)+1)+1.
    coef(0) = 0.
  ENDIF ELSE BEGIN
    sz = SIZE(coef)
    IF((sz(0) NE 1) OR (sz(2) NE 4) OR (sz(1) NE N_PARAMS(0)+1)) THEN BEGIN
      MESSAGE, 'An inappropriate COEF was supplied.',/CONTINUE
      RETURN,-1
    ENDIF
  ENDELSE
  bv = badval_in
  IF(NOT graphs) THEN BEGIN
;
;  Define GOOD pixels to be those which are good for each input image.
;  -------------------------------------------------------------------
    CASE N_PARAMS(0) OF
       1: good = WHERE(x1 NE bv(0))
       2: good = WHERE((x1 NE bv(0)) AND (x2 NE bv(1)))
       3: good = WHERE((x1 NE bv(0)) AND (x2 NE bv(1)) AND (x3 NE bv(2)))
       4: good = WHERE((x1 NE bv(0)) AND (x2 NE bv(1)) AND (x3 NE bv(2)) $
                   AND (x4 NE bv(3)))
       5: good = WHERE((x1 NE bv(0)) AND (x2 NE bv(1)) AND (x3 NE bv(2)) $
                   AND (x4 NE bv(3)) AND (x5 NE bv(4)))
       6: good = WHERE((x1 NE bv(0)) AND (x2 NE bv(1)) AND (x3 NE bv(2)) $
                   AND (x4 NE bv(3)) AND (x5 NE bv(4)) AND (x6 NE bv(5)))
       7: good = WHERE((x1 NE bv(0)) AND (x2 NE bv(1)) AND (x3 NE bv(2)) $
                   AND (x4 NE bv(3)) AND (x5 NE bv(4)) AND (x6 NE bv(5)) $
                   AND (x7 NE bv(6)))
       8: good = WHERE((x1 NE bv(0)) AND (x2 NE bv(1)) AND (x3 NE bv(2)) $
                   AND (x4 NE bv(3)) AND (x5 NE bv(4)) AND (x6 NE bv(5)) $
                   AND (x7 NE bv(6)) AND (x8 NE bv(7)))
       9: good = WHERE((x1 NE bv(0)) AND (x2 NE bv(1)) AND (x3 NE bv(2)) $
                   AND (x4 NE bv(3)) AND (x5 NE bv(4)) AND (x6 NE bv(5)) $
                   AND (x7 NE bv(6)) AND (x8 NE bv(7)) AND (x9 NE bv(8)))
      10: good = WHERE((x1 NE bv(0)) AND (x2 NE bv(1)) AND (x3 NE bv(2)) $
                   AND (x4 NE bv(3)) AND (x5 NE bv(4)) AND (x6 NE bv(5)) $
                   AND (x7 NE bv(6)) AND (x8 NE bv(7)) AND (x9 NE bv(8)) $
                   AND (x10 NE bv(9)))
    ENDCASE
    IF(good(0) EQ -1) THEN BEGIN
      outbadval = badval_in(0)
      result = x1*0. + outbadval
      RETURN,result
    ENDIF
;
;  Try to trap any potential overflows/underflows.
;  -----------------------------------------------
    data_min = 1.e30
    data_max = 1.e-30
    FOR i=1,N_PARAMS(0) DO BEGIN
      CASE i OF
         1 : tmin=MIN(x1, MAX = tmax)
         2 : tmin=MIN(x2, MAX = tmax)
         3 : tmin=MIN(x3, MAX = tmax)
         4 : tmin=MIN(x4, MAX = tmax)
         5 : tmin=MIN(x5, MAX = tmax)
         6 : tmin=MIN(x6, MAX = tmax)
         7 : tmin=MIN(x7, MAX = tmax)
         8 : tmin=MIN(x8, MAX = tmax)
         9 : tmin=MIN(x9, MAX = tmax)
        10 : tmin=MIN(x10, MAX = tmax)
      ENDCASE
      IF(tmin LT data_min) THEN data_min=tmin
      IF(tmax GT data_max) THEN data_max=tmax
    ENDFOR
    coef_min = MIN(coef(1:*), MAX = coef_max)
    IF(ABS(coef(0)) GT 1.e30) THEN BEGIN
      MESSAGE,'Potential floating point overflow condition: operation not '+$
        'performed.',/CONT
      RETURN,-1
    ENDIF
    IF((ABS(coef_max) GT 0.) AND (ABS(data_max) GT 0.)) THEN BEGIN
      IF(ALOG10(ABS(coef_max))+ALOG10(ABS(data_max)) GT 30) THEN BEGIN
        MESSAGE,'Potential floating point overflow condidtion: operation '+$
          'not performed.',/CONT
        RETURN,-1
      ENDIF
    ENDIF
    IF((ABS(coef_min) GT 0.) AND (ABS(data_min) GT 0.)) THEN BEGIN
      IF(ALOG10(ABS(coef_min))+ALOG10(ABS(data_min)) LT -30) THEN BEGIN
        MESSAGE,'Potential floating point underflow condidtion: operation '+$
          'not performed.',/CONT
        RETURN,-1
      ENDIF
    ENDIF
;
;  Perform the weighted sum.
;  -------------------------
    outbadval = MIN(badval_in)
    result = x1 * 0. + outbadval
    CASE N_PARAMS(0) OF
      1: result(good) = coef(0) + coef(1)*x1(good)
      2: result(good) = coef(0) + coef(1)*x1(good) + coef(2)*x2(good)
      3: result(good) = coef(0) + coef(1)*x1(good) + coef(2)*x2(good) + $
               coef(3)*x3(good)
      4: result(good) = coef(0) + coef(1)*x1(good) + coef(2)*x2(good) + $
               coef(3)*x3(good) + coef(4)*x4(good)
      5: result(good) = coef(0) + coef(1)*x1(good) + coef(2)*x2(good) + $
               coef(3)*x3(good) + coef(4)*x4(good) + coef(5)*x5(good)
      6: result(good) = coef(0) + coef(1)*x1(good) + coef(2)*x2(good) + $
               coef(3)*x3(good) + coef(4)*x4(good) + coef(5)*x5(good) + $
               coef(6)*x6(good)
      7: result(good) = coef(0) + coef(1)*x1(good) + coef(2)*x2(good) + $
               coef(3)*x3(good) + coef(4)*x4(good) + coef(5)*x5(good) + $
               coef(6)*x6(good) + coef(7)*x7(good)
      8: result(good) = coef(0) + coef(1)*x1(good) + coef(2)*x2(good) + $
               coef(3)*x3(good) + coef(4)*x4(good) + coef(5)*x5(good) + $
               coef(6)*x6(good) + coef(7)*x7(good) + coef(8)*x8(good)
      9: result(good) = coef(0) + coef(1)*x1(good) + coef(2)*x2(good) + $
               coef(3)*x3(good) + coef(4)*x4(good) + coef(5)*x5(good) + $
               coef(6)*x6(good) + coef(7)*x7(good) + coef(8)*x8(good) + $
               coef(9)*x9(good)
     10: result(good) = coef(0) + coef(1)*x1(good) + coef(2)*x2(good) + $
               coef(3)*x3(good) + coef(4)*x4(good) + coef(5)*x5(good) + $
               coef(6)*x6(good) + coef(7)*x7(good) + coef(8)*x8(good) + $
               coef(9)*x9(good) + coef(10)*x10(good)
    ENDCASE
;
;  Reset OUTBADVAL if some good pixels ended up having its previous value.
;  ------------------------------------------------------------------------
    w = WHERE(result(good) EQ outbadval)
    IF(w(0) NE -1) THEN BEGIN
      outbadval = MIN(result(good))-1.
      temp = result(good)
      result = result*0. + outbadval
      result(good) = temp
    ENDIF
    RETURN, result
  ENDIF ELSE BEGIN
;
;  Find out where the X arrays overlap precisely.
;  ----------------------------------------------
    xgood = FLTARR(dim1x1)
    ic = 0
    temp = N_PARAMS(0)
    FOR i=0,dim1x1-1 DO BEGIN
      w = WHERE(x1(i,1) NE bv(0))
      IF(w(0) EQ -1) THEN GOTO,badx
      IF(temp GE 2) THEN w = WHERE((x2(*,0) EQ x1(i,0)) AND (x2(*,1) NE bv(1)))
      IF(w(0) EQ -1) THEN GOTO,badx
      IF(temp GE 3) THEN w = WHERE((x3(*,0) EQ x1(i,0)) AND (x3(*,1) NE bv(2)))
      IF(w(0) EQ -1) THEN GOTO,badx
      IF(temp GE 4) THEN w = WHERE((x4(*,0) EQ x1(i,0)) AND (x4(*,1) NE bv(3)))
      IF(w(0) EQ -1) THEN GOTO,badx
      IF(temp GE 5) THEN w = WHERE((x5(*,0) EQ x1(i,0)) AND (x5(*,1) NE bv(4)))
      IF(w(0) EQ -1) THEN GOTO,badx
      IF(temp GE 6) THEN w = WHERE((x6(*,0) EQ x1(i,0)) AND (x6(*,1) NE bv(5)))
      IF(w(0) EQ -1) THEN GOTO,badx
      IF(temp GE 7) THEN w = WHERE((x7(*,0) EQ x1(i,0)) AND (x7(*,1) NE bv(6)))
      IF(w(0) EQ -1) THEN GOTO,badx
      IF(temp GE 8) THEN w = WHERE((x8(*,0) EQ x1(i,0)) AND (x8(*,1) NE bv(7)))
      IF(w(0) EQ -1) THEN GOTO,badx
      IF(temp GE 9) THEN w = WHERE((x9(*,0) EQ x1(i,0)) AND (x9(*,1) NE bv(8)))
      IF(w(0) EQ -1) THEN GOTO,badx
      IF(temp EQ 10) THEN w=WHERE((x10(*,0) EQ x1(i,0)) AND (x10(*,1) NE bv(9)))
      IF(w(0) EQ -1) THEN GOTO,badx
      xgood(ic) = x1(i,0)
      ic = ic + 1
badx:
    ENDFOR
    xgood = xgood(0:ic-1)
    dimout = ic
    IF(dimout EQ 0) THEN BEGIN
      MESSAGE,'No exact overlap in X dimension.',/CONT
      RETURN,-1
    ENDIF
    good1 = WHERE(x1(*,0) EQ xgood)
    IF(temp GE 2)  THEN good2  = WHERE(x2(*,0) EQ xgood)
    IF(temp GE 3)  THEN good3  = WHERE(x3(*,0) EQ xgood)
    IF(temp GE 4)  THEN good4  = WHERE(x4(*,0) EQ xgood)
    IF(temp GE 5)  THEN good5  = WHERE(x5(*,0) EQ xgood)
    IF(temp GE 6)  THEN good6  = WHERE(x6(*,0) EQ xgood)
    IF(temp GE 7)  THEN good7  = WHERE(x7(*,0) EQ xgood)
    IF(temp GE 8)  THEN good8  = WHERE(x8(*,0) EQ xgood)
    IF(temp GE 9)  THEN good9  = WHERE(x9(*,0) EQ xgood)
    IF(temp GE 10) THEN good10 = WHERE(x10(*,0) EQ xgood)
;
;  Try to trap any potential overflows/underflows.
;  -----------------------------------------------
    data_min = 1.e30
    data_max = 1.e-30
    FOR i=1,N_PARAMS(0) DO BEGIN
      CASE i OF
         1 : tempdata = x1(good1,1)
         2 : tempdata = x2(good2,1)
         3 : tempdata = x3(good3,1)
         4 : tempdata = x4(good4,1)
         5 : tempdata = x5(good5,1)
         6 : tempdata = x6(good6,1)
         7 : tempdata = x7(good7,1)
         8 : tempdata = x8(good8,1)
         9 : tempdata = x9(good9,1)
        10 : tempdata = x10(good10,1)
      ENDCASE
      tmin = MIN(tempdata, MAX = tmax)
      IF(tmin LT data_min) THEN data_min=tmin
      IF(tmax GT data_max) THEN data_max=tmax
    ENDFOR
    coef_min = MIN(coef(1:*), MAX = coef_max)
    IF(ABS(coef(0)) GT 1.e30) THEN BEGIN
      MESSAGE,'Potential floating point overflow condition: operation not '+$
        'performed.',/CONT
      RETURN,-1
    ENDIF
    IF((ABS(coef_max) GT 0.) AND (ABS(data_max) GT 0.)) THEN BEGIN
      IF(ALOG10(ABS(coef_max))+ALOG10(ABS(data_max)) GT 30) THEN BEGIN
        MESSAGE,'Potential floating point overflow condidtion: operation '+$
          'not performed.',/CONT
        RETURN,-1
      ENDIF
    ENDIF
    IF((ABS(coef_min) GT 0.) AND (ABS(data_min) GT 0.)) THEN BEGIN
      IF(ALOG10(ABS(coef_min))+ALOG10(ABS(data_min)) LT -30) THEN BEGIN
        MESSAGE,'Potential floating point underflow condidtion: operation '+$
          'not performed.',/CONT
        RETURN,-1
      ENDIF
    ENDIF
;
;  Perform the weighted sum.
;  -------------------------
    outbadval = MIN(badval_in)
    result = FLTARR(dimout,2)
    result(0,0) = xgood
    CASE N_PARAMS(0) OF
      1: result(0,1) = coef(0) + coef(1)*x1(good1,1)
      2: result(0,1) = coef(0) + coef(1)*x1(good1,1) + coef(2)*x2(good2,1)
      3: result(0,1) = coef(0) + coef(1)*x1(good1,1) + coef(2)*x2(good2,1) + $
           coef(3)*x3(good3,1)
      4: result(0,1) = coef(0) + coef(1)*x1(good1,1) + coef(2)*x2(good2,1) + $
           coef(3)*x3(good3,1) + coef(4)*x4(good4,1)
      5: result(0,1) = coef(0) + coef(1)*x1(good1,1) + coef(2)*x2(good2,1) + $
           coef(3)*x3(good3,1) + coef(4)*x4(good4,1) + coef(5)*x5(good5,1)
      6: result(0,1) = coef(0) + coef(1)*x1(good1,1) + coef(2)*x2(good2,1) + $
           coef(3)*x3(good3,1) + coef(4)*x4(good4,1) + coef(5)*x5(good5,1) + $
           coef(6)*x6(good6,1)
      7: result(0,1) = coef(0) + coef(1)*x1(good1,1) + coef(2)*x2(good2,1) + $
           coef(3)*x3(good3,1) + coef(4)*x4(good4,1) + coef(5)*x5(good5,1) + $
           coef(6)*x6(good6,1) + coef(7)*x7(good7,1)
      8: result(0,1) = coef(0) + coef(1)*x1(good1,1) + coef(2)*x2(good2,1) + $
           coef(3)*x3(good3,1) + coef(4)*x4(good4,1) + coef(5)*x5(good5,1) + $
           coef(6)*x6(good6,1) + coef(7)*x7(good7,1) + coef(8)*x8(good8,1)
      9: result(0,1) = coef(0) + coef(1)*x1(good1,1) + coef(2)*x2(good2,1) + $
           coef(3)*x3(good3,1) + coef(4)*x4(good4,1) + coef(5)*x5(good5,1) + $
           coef(6)*x6(good6,1) + coef(7)*x7(good7,1) + coef(8)*x8(good8,1) + $
           coef(9)*x9(good9,1)
     10: result(0,1) = coef(0) + coef(1)*x1(good1,1) + coef(2)*x2(good2,1) + $
           coef(3)*x3(good3,1) + coef(4)*x4(good4,1) + coef(5)*x5(good5,1) + $
           coef(6)*x6(good6,1) + coef(7)*x7(good7,1) + coef(8)*x8(good8,1) + $
           coef(9)*x9(good9,1) + coef(10)*x10(good10,1)
    ENDCASE
    w = WHERE(result(*,1) EQ outbadval)
    IF(w(0) NE -1) THEN outbadval = MIN(result(*,1))-1.
    RETURN, result
  ENDELSE
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



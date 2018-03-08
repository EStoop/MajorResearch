FUNCTION multiply,x1,x2,badval_in=badval_in,outbadval=outbadval
;+NAME/ONE LINE DESCRIPTION:
;    MULTIPLY multiplies two arrays together (properly handling bad pixels).
;
;DESCRIPTION:
;    MULTIPLY receives two 2-D arrays as input and returns a 2-D array as
;    output.  Any pixel which was designated as "bad" in either of the
;    input arrays will also be designated as bad in the output array.
;    Pixels which had good values in both input arrays will have, in
;    the output array, a new value which is the product of the values
;    of the corresponding pixels in the input arrays.  (NOTE - if the
;    second dimension of the arguments is 2, then they are interpreted
;    as graphs).
;
;CALLING SEQUENCE:
;    outarray = multiply( x1, x2, badval_in=[...], outbadval=[...])
;
;ARGUMENTS (I=input, O=output, []=optional)
;    x1         I   2-D arr  flt        An array to be multiplied by x2.
;    x2         I   2-D arr  flt        An array to be multiplied by x1.
;    badval_in [I]  keyword,1-D arr,flt A 1-D array that contains the bad
;                                       pixel values for x1 and x2.
;                                       [0.,0.] will be assumed if not
;                                       supplied.
;    outbadval [O]  keyword  flt        The bad pixel value for outarray.
;    outarray   O   2-D arr  flt        The "product" of x1 and x2.
;
;EXAMPLE:
;    outarray = multiply(x1,x2,badval_in=[0.,0.],outbad=outbad)
;#
;COMMON BLOCKS:  none.
;
;LIBRARY CALLS:  none.
;
;PROCEDURE:
;    If two 2-D arrays of the same size have not been passed in, then
;    put out a message and exit.  If BADVAL_IN was not supplied, then
;    assume a value of [0.,0.], and inform the user that such an
;    assumption is being made.  If a value has been assigned by the
;    caller, and if it is not a 1-D float array with two elements, then
;    put out a message and exit.  If either input image consists of all
;    bad pixels, then return an array that also consists of all bad
;    pixels.  Define good pixels to be those which are good for both
;    input images.  Try to catch any potential floating point overflows
;    or underflows by looking at the sum of the logs of the extreme
;    good values from the input images.  Put out a message and exit if
;    overflow/underflow is determined to be a possibility.  Create
;    outarray such that all non-good pixels receive a bad value, and all
;    good pixels receive the product of the corresponding values from
;    the input images.  If any good pixel has received, as a result of
;    the operation, a value which is equal to the bad pixel value, then
;    change the bad pixel value to be 1 less than the minimum of the
;    good pixel values.
;
;REVISION HISTORY:
;    Originally written by Vidia Sagar, ARC, February 1992.
;    Extensively rewritten by John Ewing, ARC, April 1992.
;    Added handling for graphs, John Ewing, ARC, June 1992.
;
;.TITLE
;Routine MULTIPLY
;-
;
;  Check the validity of the passed parameters.
;  --------------------------------------------
  IF(N_PARAMS() LT 2) THEN BEGIN
    MESSAGE,'Insufficient parameters. (Two 2-D arrays are expected).',/CONT
    RETURN,-1
  ENDIF
  s1 = SIZE(x1)
  s2 = SIZE(x2)
  IF((s1(0) NE 2) OR (s2(0) NE 2)) THEN BEGIN
    MESSAGE,'Arguments must be 2-dimensional.',/CONT
    RETURN,-1
  ENDIF
  IF((s1(2) EQ 2) AND (s2(2) EQ 2)) THEN graphs = 1 ELSE graphs = 0
  IF(NOT graphs) THEN BEGIN
    IF((s1(1) NE s2(1)) OR (s1(2) NE s2(2))) THEN BEGIN
      MESSAGE,'Arrays must be of the same size.',/CONT
      RETURN,-1
    ENDIF
  ENDIF
;
;  Check what was passed in for BADVAL_IN.
;  ---------------------------------------
  IF(NOT KEYWORD_SET(badval_in)) THEN BEGIN
    badval_in = FLTARR(2)
    IF(NOT graphs) THEN BEGIN
      PRINT,'BADVAL_IN was not set.  Assuming [0.,0.].'
    ENDIF ELSE BEGIN
      badval_in(0) = MIN(x1(*,1)) - 1.
      badval_in(1) = MIN(x2(*,1)) - 1.
    ENDELSE
  ENDIF ELSE BEGIN
    sz = SIZE(badval_in)
    IF((sz(0) NE 1) OR (sz(2) NE 4) OR (sz(1) NE 2)) THEN BEGIN
      MESSAGE, 'BADVAL_IN must be a 1-D floating array of 2 elements.',/CONTINUE
      RETURN,-1
    ENDIF
  ENDELSE
  IF(NOT graphs) THEN BEGIN
;
;  Handle the case in which one or more array(s) consist of all bad pixels.
;  ------------------------------------------------------------------------
    good1 = WHERE(x1 NE badval_in(0))
    IF(good1(0) EQ -1) THEN BEGIN
      outbadval = MIN(x1) - 1.
      outarray = x1 * 0 + outbadval
      RETURN,outarray
    ENDIF
    good2 = WHERE(x2 NE badval_in(1))
    IF(good2(0) EQ -1) THEN BEGIN
      outbadval = MIN(x1) - 1.
      outarray = x1 * 0 + outbadval
      RETURN,outarray
    ENDIF
;
;  Catch any potential floating point overflows or underflows.
;  -----------------------------------------------------------
    good = WHERE((x1 NE badval_in(0)) AND (x2 NE badval_in(1)))
    min1 = MIN(ABS(x1(good)))
    min2 = MIN(ABS(x2(good)))
    IF((min1 GT 0) AND (min2 GT 0)) THEN BEGIN
       IF(ALOG10(min1) + ALOG10(min2) LT -35) THEN BEGIN
       MESSAGE,' Floating point underflow anticipated. '+$
                ' Operation not performed.',/continue
        RETURN,-1
     ENDIF 
    ENDIF
    max1 = MAX(ABS(x1))
    max2 = MAX(ABS(x2))
    IF((max1 GT 0) AND (max2 GT 0)) THEN BEGIN
      IF(ALOG10(max1) + ALOG10(max2) GT 35) THEN BEGIN
        MESSAGE, ' Floating point overflow anticipated. '+$
                 ' Operation not performed.',/continue
        RETURN,-1
      ENDIF
    ENDIF
;
;  Perform the multiplication.
;  ---------------------------
    outbadval = MIN(badval_in)
    outarray = x1 * 0. + outbadval
    outarray(good) = x1(good) * x2(good)
;
;  Reset OUTBADVAL if some good pixels ended up having its previous value.
;  -----------------------------------------------------------------------
    w = WHERE(outarray(good) EQ outbadval)
    IF(w(0) NE -1) THEN BEGIN
      outbadval = MIN(outarray(good))-1.
      temp = outarray(good)
      outarray = outarray*0. + outbadval
      outarray(good) = temp
    ENDIF
    RETURN, outarray
  ENDIF ELSE BEGIN
;
;  Find out where the X arrays overlap precisely.
;  ----------------------------------------------
    dim1 = s1(1)
    dim2 = s2(1)
    match = INTARR(dim2)
    FOR i=0,dim2-1 DO match(i) = WHERE(x1(*,0) EQ x2(i,0))
    dimout = 0
    FOR i=0,dim2-1 DO IF(match(i) NE -1) THEN dimout = dimout + 1
    IF(dimout EQ 0) THEN BEGIN
      MESSAGE,'No exact overlap in X dimension.',/CONT
      RETURN,-1
    ENDIF
;
;  Perform the multiplication.
;  ---------------------------
    xout = FLTARR(dimout)
    yout = FLTARR(dimout)
    outbadval = badval_in(0)
    good_flag = INTARR(dimout)
    ic = 0
    FOR i=0,dim2-1 DO BEGIN
      IF(match(i) NE -1) THEN BEGIN
        IF((x1(match(i),1) NE badval_in(0)) AND $
          (x2(i,1) NE badval_in(1))) THEN BEGIN
          temp = alog10(x1(match(i),1)) + alog10(x2(i,1))
          IF((temp LT -30) OR (temp GT 30)) THEN BEGIN
            result = outbadval
            good_flag(ic) = 0
          ENDIF ELSE BEGIN
            result = x1(match(i),1) * x2(i,1)
              good_flag(ic) = 1
          ENDELSE
        ENDIF ELSE BEGIN
          result = outbadval
          good_flag(ic) = 0
        ENDELSE
        xout(ic) = x2(i,0)
        yout(ic) = result
        ic = ic +1
      ENDIF
    ENDFOR
;
;  Redefine OUTBADVAL if any good points got assigned that value.
;  --------------------------------------------------------------
    good = WHERE(good_flag EQ 1)
    IF(good(0) NE -1) THEN BEGIN
      w = WHERE(yout(good) EQ outbadval)
      IF(w(0) NE -1) THEN BEGIN
        outbadval = MIN(yout(good)) - 1.
        bad = WHERE(good_flag EQ 0)
        IF(bad(0) NE -1) THEN yout(bad) = outbadval
      ENDIF
    ENDIF
    RETURN, [[xout],[yout]]
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



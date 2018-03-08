FUNCTION squareroot,input,badval_in=badval_in,outbadval=outbadval
;+NAME/ONE LINE DESCRIPTION:
;    SQUAREROOT takes the square root of a map (properly handling bad pixels).
;
;DESCRIPTION:
;    SQUAREROOT receives a 2-D array as input and returns a 2-D array as
;    output.  Any pixel which was designated as "bad" in the input
;    array will also be designated as bad in the output array.  Also,
;    any good pixel in the input array which has a value less than
;    zero, will be designated as a bad pixel in the output array.
;    All other pixels in the output array will receive a value that
;    is the square root of the value of the corresponding pixel in the
;    input array.  (NOTE - if the second dimension of the argument is 2,
;    then it is interpreted as a graph).
;
;CALLING SEQUENCE:
;    result = squareroot( input, badval_in=[...], outbadval=[...] )
;
;ARGUMENTS (I=input, O=output, []=optional)
;    input      I   2-D arr  flt        Input map.
;    badval_in [I]  keyword  flt        Bad pixel value for input.
;                                       0 will be assumed if not supplied.
;    outbadval [O]  keyword  flt        The bad pixel value for result.
;    result     O   2-D arr  flt        The "square root" of input.
;
;EXAMPLE:
;    result = squareroot(input,badval_in=0.,outbad=outbad)
;#
;COMMON BLOCKS:  none.
;
;LIBRARY CALLS:  none.
;
;PROCEDURE:
;    If a 2-D array has not been passed in, then put out a message and
;    exit.  If BADVAL_IN was not supplied, then assume a value of 0, and
;    inform the user that such an assumption is being made.  If a value
;    has been assigned by the caller, and if it is not a floating point
;    scalar, then put out a message and exit.  If the input image
;    consists of all bad pixels, then return an array that also consists
;    of all bad pixels.  Define good pixels for the output array to be
;    those which are good for the input image and which have a value in
;    the input image which is greater than or equal to 0.  Create
;    result such that all non-good pixels receive a bad value, and all
;    good pixels receive the square root of the corresponding value from
;    the input image.  If any good pixel has received, as a result of
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
;Routine SQUAREROOT
;-
;
;  Check what was passed in for input.
;  -----------------------------------
  IF(N_PARAMS(0) LT 1) THEN BEGIN
    MESSAGE,'Insufficient parameters. (A 2-D array is expected).',/CONT
    RETURN,-1
  ENDIF
  sz = SIZE(input)
  IF(sz(0) NE 2) THEN BEGIN
    MESSAGE,'The argument (input) must be 2-dimensional.',/CONT
    RETURN,-1
  ENDIF
  IF(sz(2) EQ 2) THEN graphs = 1 ELSE graphs = 0
  dim1 = sz(1)
;
;  Check what was passed in for BADVAL_IN.
;  ---------------------------------------
  sz = SIZE(badval_in)
  IF((sz(0)+sz(1)) EQ 0) THEN BEGIN
    IF(NOT graphs) THEN BEGIN
      PRINT,'BADVAL_IN was not set.  Assuming a value of 0.'
      badval_in = 0.
    ENDIF ELSE BEGIN
      badval_in = MIN(input(*,1)) - 1.
    ENDELSE
  ENDIF ELSE BEGIN
    sz = SIZE(badval_in)
    IF((sz(0) NE 0) OR ((sz(1) NE 4) AND (sz(1) NE 2))) THEN BEGIN
      MESSAGE, 'BADVAL_IN must be a floating point scalar.',/CONTINUE
      RETURN,-1
    ENDIF
  ENDELSE
;
;  Create RESULT, in which all good pixels from input which had values
;  greater than or equal to 0 will have, in result, the square root of
;  the value they had in input.  All other pixels will have the value of
;  BADVAL_IN.
;  ---------------------------------------------------------------------
  outbadval = badval_in
  IF(NOT graphs) THEN BEGIN
    good = WHERE((input NE badval_in) AND (input GE 0.))
    result = input * 0. + outbadval
    IF(good(0) NE -1) THEN result(good) = SQRT(input(good)) $
                      ELSE RETURN,result
;
;  Reset OUTBADVAL if some good pixels ended up having its previous value.
;  -----------------------------------------------------------------------
    w = WHERE(result(good) EQ outbadval)
    IF(w(0) NE -1) THEN BEGIN
      outbadval = MIN(result(good))-1.
      temp = result(good)
      result = result*0. + outbadval
      result(good) = temp
    ENDIF
    RETURN,result
  ENDIF ELSE BEGIN
;
;  Handle the special case of a graph.
;  -----------------------------------
    good = WHERE((input(*,1) NE badval_in) AND (input(*,1) GE 0.))
    bad  = WHERE((input(*,1) EQ badval_in) OR  (input(*,1) LT 0.))
    yout = FLTARR(dim1) + outbadval
    IF(good(0) NE -1) THEN yout(good) = SQRT(input(good,1)) $
                      ELSE RETURN,[[input(*,0)],[yout]]
;
;  Redefine OUTBADVAL if any good points got assigned that value.
;  --------------------------------------------------------------
    w = WHERE(yout(good) EQ outbadval)
    IF(w(0) NE -1) THEN BEGIN
      outbadval = MIN(yout(good)) - 1.
      IF(bad(0) NE -1) THEN yout(bad) = outbadval
    ENDIF
    RETURN, [[input(*,0)],[yout]]
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



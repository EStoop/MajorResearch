FUNCTION hiseqmap,inmap,badval=badval,data_min=data_min,data_max=data_max
;+ NAME/ONE LINE DESCRIPTION:
;    HISEQMAP returns a histogram equalized version of an input image.
;
;DESCRIPTION:
;    This routine requires one argument: a 2-D array of numerical data.
;    The keyword BADVAL should be supplied with an appropriate bad pixel
;    value (if known), otherwise a bad pixel value of 0 is assumed.
;    The output image is of the same type as the input image.  Bad pixels
;    in the output image will have the same value which they had in the
;    input image.  Good pixels in the output image will be histogram
;    equalized and will range from 0 to 255.
;
;CALLING SEQUENCE:
;    outmap = hiseqmap(inmap, [badval=...])
;
;ARGUMENTS (I=input, O=output, []=optional):
;    inmap      I   2-D arr, numerical  Input map
;    badval    [I]  keyword  flt        Flag value that identifies bad
;                                       pixels.  If not present, then a
;                                       value of 0. will be assumed.
;    data_min  [I]  keyword  flt        Data minimum.
;    data_max  [I]  keyword  flt        Data maximum.
;    outmap     O   2-D arr             Histogram-equalized map, of the
;                                       same type as inmap.
;
;EXAMPLE:
;    outmap = hiseqmap(inmap, bad=0.)
;#
;PROCEDURE:
;    Check the supplied arguments.  Exit with an error message if a 2-D
;    numerical array was not supplied for inmap.  If a value for badval
;    was not supplied, then assume a value of 0. and print out a statement
;    telling the user that such an assumption is being made.  Identify
;    what are the good pixels in inmap.  If there are none, then print
;    out a message and exit.  Otherwise, next check what is the dynamic
;    range.  If the dynamic range is 0 then print out a message and exit.
;    Otherwise call HIST_EQUAL with a 1-D array of good pixel values as
;    the argument.  Set outmap to be equal to inmap, and then set the
;    good pixels in outmap to be equal to the result returned from
;    HIST_EQUAL.  Then exit, returning the value of outmap.
;
;LIBRARY CALLS:  hist_equal
;
;REVISION HISTORY:
;    Written by John Ewing, ARC            16 April 1992
;
;.TITLE
;Routine HISEQMAP
;-
;
;  Check supplied arguments.
;  -------------------------
  IF(N_PARAMS(0) NE 1) THEN BEGIN
    MESSAGE,'One argument (a 2-D array) must be supplied.',/CONT
    RETURN,-1
  ENDIF
  sz = SIZE(inmap)
  IF(sz(0) NE 2) THEN BEGIN
    MESSAGE,'The argument is not a 2-D array.',/CONT
    RETURN,-1
  ENDIF
  IF((sz(3) LT 1) OR (SZ(3) GT 5)) THEN BEGIN
    MESSAGE,'The argument is of an invalid type.',/CONT
    RETURN,-1
  ENDIF
  sz = SIZE(badval)
  IF((sz(0) + sz(1)) EQ 0) THEN BEGIN
    PRINT,'A value for BADVAL was not supplied.  '+$
      'HISEQMAP will assume a BADVAL of 0.'
    badval=0.
  ENDIF
;
;  Determine what are the good pixels and what's the dynamic range.
;  ----------------------------------------------------------------
  IF((NOT KEYWORD_SET(data_min)) OR (NOT KEYWORD_SET(data_max))) THEN BEGIN
    good = WHERE(inmap NE badval)
    IF(good(0) EQ -1) THEN BEGIN
      MESSAGE,'The supplied image consists of all bad pixels.',/CONT
      RETURN,-1
    ENDIF
    data_min = MIN(inmap(good), MAX = data_max)
  ENDIF ELSE BEGIN
    IF(data_min LT data_max) THEN BEGIN
      MESSAGE,'The supplied image consists of all bad pixels.',/CONT
      RETURN,-1
    ENDIF
  ENDELSE
  IF(data_min EQ data_max) THEN BEGIN
    MESSAGE,'There is no dynamic range in the supplied image.',/CONT
    RETURN,-1
  ENDIF
;
;  Call HIST_EQUAL to histogram equalize the good pixels.
;  ------------------------------------------------------
  binsize = (data_max-data_min)/8192.
  valid = inmap(good)
  heq = hist_equal(valid,minv=data_min,maxv=data_max,bin=binsize)
  result = inmap
  result(good) = heq
  RETURN,result
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



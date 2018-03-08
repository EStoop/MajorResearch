FUNCTION enhanmap,inmap,badval=badval,type=type,data_min=data_min, $
  data_max=data_max
;+ NAME/ONE LINE DESCRIPTION:
;    ENHANMAP returns an edge enhanced version of an input image.
;
;DESCRIPTION:
;    This routine requires one argument: a 2-D array of numerical data.
;    The keyword BADVAL should be supplied with an appropriate bad pixel
;    value (if known), otherwise a bad pixel value of 0 is assumed.
;    The keyword TYPE should be supplied with a string which is either
;    "ROBERTS" or "SOBEL".  If no value is supplied, then a default
;    value of "ROBERTS" is used.  The ouput image is of type integer.
;    Bad pixels in the output image will have a value of 0.  Pixels in
;    the output image which have been contaminated by a neighboring
;    bad pixel will be assigned a value of 0.
;
;CALLING SEQUENCE:
;    outmap = enhanmap(inmap, [badval=...], [type=...])
;
;ARGUMENTS (I=input, O=output, []=optional):
;    inmap      I   2-D arr, numerical  Input map
;    badval    [I]  keyword  flt        Flag value that identifies bad
;                                       pixels.  If not present, then a
;                                       value of 0. will be assumed.
;    type      [I]  keyword  str        Type of edge enhancement. Should
;                                       be either "ROBERTS" or "SOBEL".
;    data_min  [I]  keyword  flt        Data minimum.
;    data_max  [I]  keyword  flt        Data maximum.
;    outmap     O   2-D arr  int        Edge enhanced map.
;
;EXAMPLE:
;    outmap = enhanmap(inmap, bad=0., type='SOBEL')
;#
;PROCEDURE:
;    Check the supplied arguments.  Exit with an error message if a 2-D
;    numerical array was not supplied for inmap.  If a value for BADVAL
;    was not supplied, then assume a value of 0. and print out a statement
;    telling the user that such an assumption is being made.  If a value
;    for TYPE was not supplied, then assume a value of "ROBERTS" and
;    print out a statement telling the user that such an assumption is
;    being made.  If any value was supplied other than "ROBERTS" or
;    "SOBEL", then print out an error message and exit.  Identify what
;    are the good pixels in inmap.  If there are none, then print out a
;    message and exit.  Otherwise, next check what is the dynamic range.
;    If the dynamic range is 0 then print out a message and exit.
;    Otherwise call either ROBERTS or SOBEL to produce an edge enhanced
;    image (outmap).  Set any pixels in outmap which are neighbors to a
;    bad pixel to 0, since their values were messed up by the bad pixel.
;    Then exit, returning the value of outmap.
;
;LIBRARY CALLS:  roberts, sobel
;
;REVISION HISTORY:
;    Written by John Ewing, ARC            16 April 1992
;
;.TITLE
;Routine ENHANMAP
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
  dim1 = sz(1)
  dim2 = sz(2)
  IF((sz(3) LT 1) OR (SZ(3) GT 5)) THEN BEGIN
    MESSAGE,'The argument is of an invalid type.',/CONT
    RETURN,-1
  ENDIF
  sz = SIZE(badval)
  IF((sz(0) + sz(1)) EQ 0) THEN BEGIN
    PRINT,'A value for BADVAL was not supplied.  '+$
      'ENHANMAP will assume a BADVAL of 0.'
    badval=0.
  ENDIF
  sz = SIZE(type)
  IF((sz(0) + sz(1)) EQ 0) THEN BEGIN
    PRINT,'A value for TYPE was not supplied.  '+$
      'ENHANMAP will assume a TYPE of "ROBERTS".'
    type='ROBERTS'
  ENDIF
  sz = SIZE(type)
  IF((sz(0) NE 0) OR (sz(1) NE 7)) THEN BEGIN
    MESSAGE,'TYPE must be a string (either "ROBERTS" or "SOBEL").',/CONT
    RETURN,-1
  ENDIF
  type=STRUPCASE(type)
  IF((type NE 'SOBEL') AND (type NE 'ROBERTS')) THEN BEGIN
    MESSAGE,'"ROBERTS" or "SOBEL" are the only allowed values for TYPE.',/CONT
    RETURN,-1
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
    IF(data_min GT data_max) THEN BEGIN
      MESSAGE,'The supplied image consists of all bad pixels.',/CONT
      RETURN,-1
    ENDIF
  ENDELSE
  IF(data_min EQ data_max) THEN BEGIN
    MESSAGE,'There is no dynamic range in the supplied image.',/CONT
    RETURN,-1
  ENDIF
;
;  Call either ROBERTS or SOBEL to perform the edge enhancement.
;  -------------------------------------------------------------
  IF(type EQ 'ROBERTS') THEN outmap = ROBERTS(inmap*1000./data_max) $
                        ELSE outmap = SOBEL(inmap*1000./data_max)
;
;  Zero out any pixels which were contaminated by a neighboring bad pixel.
;  -----------------------------------------------------------------------
  bad = WHERE(inmap EQ badval)
  IF(bad(0) NE -1) THEN BEGIN
    temp = inmap
    signal = MIN(outmap)-10
    temp(bad) = signal
    FOR j=1,dim2-2 DO FOR i=1,dim1-2 DO $
      IF(MIN(temp(i-1:i+1,j-1:j+1)) EQ signal) THEN outmap(i,j) = 0
    FOR i=0,dim1-2 DO BEGIN
      IF(temp(i+1,0) EQ signal) THEN outmap(i,0) = 0
      IF(temp(i+1,dim2-1) EQ signal) THEN outmap(i,dim2-1) = 0
    ENDFOR
    FOR i=1,dim1-1 DO BEGIN
      IF(temp(i-1,0) EQ signal) THEN outmap(i,0) = 0
      IF(temp(i-1,dim2-1) EQ signal) THEN outmap(i,dim2-1) = 0
    ENDFOR
    FOR j=0,dim2-2 DO BEGIN
      IF(temp(0,j+1) EQ signal) THEN outmap(0,j) = 0
      IF(temp(dim1-1,j+1) EQ signal) THEN outmap(dim1-1,j) = 0
    ENDFOR
    FOR j=1,dim2-1 DO BEGIN
      IF(temp(0,j-1) EQ signal) THEN outmap(0,j) = 0
      IF(temp(dim1-1,j-1) EQ signal) THEN outmap(dim1-1,j) = 0
    ENDFOR
  ENDIF
  RETURN,outmap
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



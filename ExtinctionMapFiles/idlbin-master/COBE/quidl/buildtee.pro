        FUNCTION buildtee, f0, f1, f2, f3, f4, f5, facenums=facenums, $
                                 badval=badval, sixpack=sixpack
;
;+NAME/ONE-LINE DESCRIPTION:
;     BUILDTEE takes 1-6 cube faces and places them into an unfolded T
;
; DESCRIPTION:
;    This function takes up to six square arrays representing skycube
;    faces and places them into the appropriate position on the unfolded
;    full sky "right T", thus enabling the user to build a full sky map  
;    out of single faces.  The faces can be supplied in any order in the
;    calling sequence, provided that you also supply an optional
;    FACENUMS array telling which face is which.  The FACENUMS array is
;    mandatory if fewer than all six faces are being supplied, in order
;    to know which faces they are.  If all six faces are supplied and
;    no FACENUMS array is specified, the input arrays are assumed to be
;    in the order 0,1,2,3,4,5.
;
;    Any numerical data type can be accommodated; the output will always
;    be a full sky right T with the same data type as the input.  The 
;    input faces may also be spectral cubes (e.g. FIRAS data).  The user
;    may also supply an optional BADVAL bad pixel value that will fill 
;    all unoccupied positions of the output array. 
;
;    If the /SIXPACK qualifier is used then the routine builds a 3 x 2
;    "sixpack" instead of a right T.
;
; CALLING SEQUENCE:
;    outmap = BUILDTEE (face0, [face1,...,face5], [facenums=facenums], $
;                                     [badval=badval], [/sixpack])
;
; ARGUMENTS: (I=input, O=output, []=optional)
;
;    OUTMAP     O   any type   Output unfolded right-T (all six faces)
;                                 with same data type as input faces
;    FACE0      I   any type   Input square array (possibly with a 3rd, 
;                                 spectrum dimension), holding one face 
;                                 of sky cube.  Need NOT be face 0.
;   [FACE1..5] [I]  any type   Other cube faces, in any order
;   [FACENUMS] [I]  int [arr]  Vector containing face numbers of 
;                                 corresponding input cube faces
;   [BADVAL]   [I]  int        Bad or unoccupied pixel value, default=0 
;   [/SIXPACK] [I]  keyword    Flag to pack output into sixpack instead
;                                 of right T
;
;
; WARNINGS:
;    It is up to the user to make sure that his input faces are indeed
;    square in the x and y dimensions.
;
; EXAMPLES:
;    If you have all six faces in order f0...f5 and these contain zeros  
;    as their flag value, then the call
;
;    UIDL> SKYMAP = BUILDTEE (F0,F1,F2,F3,F4,F5)
;
;    ...will build a complete skymap.  You could (perversely) shuffle the
;    order of the input arrays and still build a correct skymap by
;    providing the order in the FACENUMS array, e.g.:
;
;    UIDL> SKYMAP = BUILDTEE (F1,F4,F2,F5,F3,F0,FACENUMS=[1,4,2,5,3,0])
;
;    Now suppose you had only two faces, called GALCEN for the galactic
;    center (face 4) and ECLPOLE for the North Ecliptic Pole (face 0), 
;    your sentinal value is -99.  You can build a full skymap with
;    data in the appropriate two faces by calling:
;
;    UIDL> SKY = BUILDTEE (GALCEN, ECLPOLE, FACENUMS=[4,0], BADVAL=-99)
;
;    Most of the map would then be filled with the sentinel value.
;#
;  COMMON BLOCKS:  None
;
;  PROCEDURE AND PROGRAMMING NOTES:
;      The routine makes extensive use of CASE statements in order to
;      set the output data type to be the same as the input, and in order
;      to be able to step through an arbitrary number of parameters in
;      arbitrary order.  Note also that all work is done on 3-D arrays
;      (spectral cubes).  2-D data is first REFORMed into 3-D arrays with
;      a degenerate 3rd dimension (1), the work is done, and both the
;      input and output arrays are REFORMED back to 2-D.
;
;  LIBRARY CALLS:  None
;
;  MODIFICATION HISTORY:
;     Written by Rich Isaacman, General Sciences Corp.  8 May 1992 
;     Modified by RBI to include /SIXPACK switch       15 May 1992
;     SPR 9705, CCR 596
;
;.TITLE
;Routine BUILDTEE
;-
;
;  Establish error condition to return to calling routine in case of
;  trouble.  First check that there is some input!
;
        ON_ERROR, 2
        IF (N_PARAMS() EQ 0) THEN MESSAGE, "No input faces provided!"
;
;  Default number of faces to 6 if facenums array is not specified,
;  then check to see if numbern of input faces is same as number of
;  elements in facenums.
;
        nfaces = N_ELEMENTS(facenums)
        IF (nfaces EQ 0) THEN BEGIN
           nfaces = 6
           facenums = [0, 1, 2, 3, 4, 5]
        ENDIF
        IF (nfaces NE N_PARAMS()) THEN MESSAGE, $
                   "Input data does not match expected number of faces"
 ;
 ;  Convert 2-D arrays to 3-D temporarily
 ;
        dims = SIZE(f0)
        IF ((dims(0) lt 2) or (dims(0) gt 3)) THEN MESSAGE, $
                                   "Input faces must be 2-D or 3-D"
        
        xlen = dims(1)                   ; the size of a face
        IF (dims(0) eq 2) THEN BEGIN
                                 f0 = REFORM (f0, xlen, xlen, 1) 
           IF (nfaces ge 2) THEN f1 = REFORM (f1, xlen, xlen, 1) 
           IF (nfaces ge 3) THEN f2 = REFORM (f2, xlen, xlen, 1)
           IF (nfaces ge 4) THEN f3 = REFORM (f3, xlen, xlen, 1) 
           IF (nfaces ge 5) THEN f4 = REFORM (f4, xlen, xlen, 1)
           IF (nfaces ge 6) THEN f5 = REFORM (f5, xlen, xlen, 1) 
        ENDIF
        dims = SIZE(f0)
        IF KEYWORD_SET(sixpack) THEN BEGIN
            xfaces = 3
            yfaces = 2
        ENDIF ELSE BEGIN
            xfaces = 4
            yfaces = 3
        ENDELSE
        xout = xfaces * xlen
        yout = yfaces * xlen
        IF (N_ELEMENTS(badval) EQ 0) THEN badval = 0.
;
;  Create the output array to be same data type as input faces
;
        CASE dims(4) OF 1:  outmap = bytarr(xout, yout, dims(3)) + badval  
                        2:  outmap = intarr(xout, yout, dims(3)) + badval
                        3:  outmap = lonarr(xout, yout, dims(3)) + badval
                        4:  outmap = fltarr(xout, yout, dims(3)) + badval
                        5:  outmap = dblarr(xout, yout, dims(3)) + badval
                        6:  outmap = complexarr(xout, yout, dims(3)) + badval
                     ELSE:  MESSAGE, "Non-numeric input data type!"
        ENDCASE
;
;  Set up the coordinates of each of the six face positions within the
;  unfolded right T.  "F0XLO" means "Face 0, low X coord", i.e. the x 
;  coord of the left edge.  Similarly, "F4YHI" is the Y coordinate of the 
;  top of face 4, etc.  Each face is XLEN on a side.
;
        IF KEYWORD_SET(sixpack) THEN BEGIN
            f0xlo = 2 * xlen          &   f0ylo = xlen  
            f0xhi = f0xlo + xlen - 1  &   f0yhi = f0ylo + xlen - 1  
            f1xlo = f0xlo             &   f1ylo = 0  
            f1xhi = f0xhi             &   f1yhi = xlen - 1      
            f2xlo = xlen              &   f2ylo = 0
            f2xhi = f1xhi - xlen      &   f2yhi = f1yhi
            f3xlo = 0                 &   f3ylo = 0
            f3xhi = xlen - 1          &   f3yhi = f1yhi
            f4xlo = 0                 &   f4ylo = xlen 
            f4xhi = xlen - 1          &   f4yhi = f0yhi
            f5xlo = xlen              &   f5ylo = xlen
            f5xhi = f2xhi             &   f5yhi = f0yhi

        ENDIF ELSE BEGIN
            f0xlo = 3 * xlen          &   f0ylo = 2 * xlen  
            f0xhi = f0xlo + xlen - 1  &   f0yhi = f0ylo + xlen - 1  
            f1xlo = f0xlo             &   f1ylo = xlen  
            f1xhi = f0xhi             &   f1yhi = f1ylo + xlen - 1      
            f2xlo = f1xlo - xlen      &   f2ylo = f1ylo
            f2xhi = f1xhi - xlen      &   f2yhi = f1yhi
            f3xlo = f2xlo - xlen      &   f3ylo = f1ylo
            f3xhi = f2xhi - xlen      &   f3yhi = f1yhi
            f4xlo = 0                 &   f4ylo = f1ylo 
            f4xhi = xlen - 1          &   f4yhi = f1yhi
            f5xlo = f0xlo             &   f5ylo = 0
            f5xhi = f0xhi             &   f5yhi = xlen - 1
        ENDELSE
;
;  Now do the work.  Step through the FACENUMS array to see which face
;  the a given parameter corresponds to, and stuff it into the output.
;
        FOR j=0, nfaces-1 DO BEGIN
           status = EXECUTE("fj = f" + STRTRIM(STRING(j),1))
           CASE facenums(j) OF 0:  outmap(f0xlo:f0xhi, f0ylo:f0yhi, *) = fj
                               1:  outmap(f1xlo:f1xhi, f1ylo:f1yhi, *) = fj
                               2:  outmap(f2xlo:f2xhi, f2ylo:f2yhi, *) = fj
                               3:  outmap(f3xlo:f3xhi, f3ylo:f3yhi, *) = fj
                               4:  outmap(f4xlo:f4xhi, f4ylo:f4yhi, *) = fj
                               5:  outmap(f5xlo:f5xhi, f5ylo:f5yhi, *) = fj
                            ELSE:  MESSAGE, "Illegal face number (must be 0-5)"
           ENDCASE
        ENDFOR
;
;  Go straight home if it was a 3-D data cube all along; otherwise,
;  REFORM it back to 2-D before returning.
;
        IF (dims(0) EQ 3) THEN RETURN, outmap 
        
                              f0 = REFORM (f0, xlen, xlen) 
        IF (nfaces ge 2) THEN f1 = REFORM (f1, xlen, xlen) 
        IF (nfaces ge 3) THEN f2 = REFORM (f2, xlen, xlen)
        IF (nfaces ge 4) THEN f3 = REFORM (f3, xlen, xlen) 
        IF (nfaces ge 5) THEN f4 = REFORM (f4, xlen, xlen)
        IF (nfaces ge 6) THEN f5 = REFORM (f5, xlen, xlen) 
        RETURN, REFORM(outmap, 4*xlen, 3*xlen)
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



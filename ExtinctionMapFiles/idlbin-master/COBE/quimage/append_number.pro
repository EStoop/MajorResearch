FUNCTION append_number, intitle
;+
;  APPEND_NUMBER - a UIMAGE-specific routine.  This function takes a
;  string as input and returns a string as output.
;     Let's imagine that the input string is "Histogram plot".  Then if
;  there are no existing UIMAGE data-objects whose title begins with that
;  string, then this routine will return "Histogram plot [1]".  If there
;  are some other data-objects with titles "Histogram plot [1]" and
;  "Histogram plot [2]", then this routine will return "Histogram plot [3]".
;     This routine always returns a string which is equal to the input string
;  appended by characters representing a space, a left bracket, a number, and
;  a right bracket.  The value of the number is set to be 1 if there are no
;  data-objects whose title begins with the input string, or else set to
;  be one plus the largest number contained within any data-object's title,
;  in which the first part of that title matched the input string.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 06 93  Add more comments.  J Ewing
;  SPR 10456  Jan 13 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 11101  Jun 28 93  Add ability to handle 3D objects. J. Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;----------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
;
;  Search through the UIMAGE data environment to see what is the title of
;  each existing data-object (exlude any appendages to the titles of the
;  form [N]).  Whenever a data-object has a title (after the [N] tail has
;  been removed) which is equal to INTITLE, then set NUM to be the number
;  that is bounded by the brackets for that title string (use a value of
;  zero if there is no bracketed number).  Set MAXNUM to be the highest
;  value of NUM that is encountered.
;
;  Search through MAPs and FACEs.
;  ------------------------------
  maxnum = 0
  type = ['map', 'face']
  FOR it = 0, 1 DO BEGIN
    FOR res = 6, 9 DO BEGIN
      arr_struc = type(it) + string(res, '(i1)')
      i = EXECUTE('sz = SIZE(' + arr_struc + ')')
      IF(sz(2) EQ 8) THEN BEGIN
        num_entries = sz(3)
        FOR i = 0, num_entries - 1 DO BEGIN
          arr_elem = arr_struc + '(' + STRTRIM(STRING(i), 2) + ')'
          j = EXECUTE('title = ' + arr_elem + '.title')
          toptitle = STRMID(title, 0, STRLEN(intitle))
          tail = STRMID(title, STRLEN(intitle), STRLEN(title))
          IF(toptitle EQ intitle) THEN BEGIN
            ic1 = STRPOS(tail, '[')
            ic2 = STRPOS(tail, ']')
            IF((ic1 NE -1) AND (ic2 NE -1)) THEN BEGIN
              num = FIX(STRMID(tail, ic1 + 1, ic2 - ic1 - 1))
              IF(num GT maxnum) THEN maxnum = num
            ENDIF
          ENDIF
        ENDFOR
      ENDIF
    ENDFOR
  ENDFOR
;
;  Search through PROJ_MAPs, GRAPHs, and ZOOMEDs.
;  ----------------------------------------------
  type = ['proj_map', 'proj2_map','graph', 'zoomed']
  FOR it = 0, 3 DO BEGIN
    arr_struc = type(it)
    i = EXECUTE('sz = SIZE(' + arr_struc + ')')
    IF(sz(2) EQ 8) THEN BEGIN
      num_entries = sz(3)
      FOR i = 0, num_entries - 1 DO BEGIN
        arr_elem = arr_struc + '(' + STRTRIM(STRING(i), 2) + ')'
        j = EXECUTE('title = ' + arr_elem + '.title')
        toptitle = STRMID(title, 0, STRLEN(intitle))
        tail = STRMID(title, STRLEN(intitle), STRLEN(title))
        IF(toptitle EQ intitle) THEN BEGIN
          ic1 = STRPOS(tail, '[')
          ic2 = STRPOS(tail, ']')
          IF((ic1 NE -1) AND (ic2 NE -1)) THEN BEGIN
            num = FIX(STRMID(tail, ic1 + 1, ic2 - ic1 - 1))
            IF(num GT maxnum) THEN maxnum = num
          ENDIF
        ENDIF
      ENDFOR
    ENDIF
  ENDFOR
;
; Search through 3D objects
    siz = SIZE(object3d)
    index3d = 0
    IF (siz(0) NE 0) THEN BEGIN
      WHILE(object3d(index3d).inuse EQ 1) DO BEGIN
        title=object3d(index3d).title
        toptitle = STRMID(title, 0, STRLEN(intitle))
        tail = STRMID(title, STRLEN(intitle), STRLEN(title))
        IF(toptitle EQ intitle) THEN BEGIN
          ic1 = STRPOS(tail, '[')
          ic2 = STRPOS(tail, ']')
          IF((ic1 NE -1) AND (ic2 NE -1)) THEN BEGIN
            num = FIX(STRMID(tail, ic1 + 1, ic2 - ic1 - 1))
            IF(num GT maxnum) THEN maxnum = num
          ENDIF
        ENDIF
        index3d = index3d + 1
        IF (index3d EQ 3) THEN GOTO, done 
      ENDWHILE
    ENDIF
;
;  Return INTITLE with an appendage which is a string that consists
;  of MAXNUM+1 bounded by brackets.
;  ----------------------------------------------------------------
done:
  RETURN, intitle + ' [' + STRTRIM(STRING(maxnum + 1), 2) + ']'
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



PRO trace_3d_links, index3d, nlinks, names, unlink = unlink
;+
;  TRACE_3D_LINKS - a UIMAGE-specific routine.  This routine identifies
;  all the 2-D data-objects which are linked to a specified 3-D data-
;  object.  INDEX3D is an input argument, NLINKS and NAMES are output
;  arguments.
;     INDEX3D is a number, either 0, 1, or 2, which identifies the 3-D
;  data-object.  This routine will set NLINKS to the number of 2-D
;  objects which are linked to the specified 3-D object (i.e., the value
;  of their .LINK3D field is INDEX3D+1).  NAMES will be set to a 1-D
;  array of strings of size NLINKS.  Each element in NAMES will be something
;  like "MAP6(0)", identifying the name of a 2-D object that is linked
;  to the 3-D object.
;     The UNLINK keyword severs the links between the children 2-D objects
;  and the parent 3-D object.  This is desirable for the situation in
;  which the parent 3-D object is being deleted, in which it would be
;  undesirable for the children 2-D objects to report a link to a non-
;  existant parent.  Those links are severed by setting the .LINK3D
;  field to 0 in each child structure.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10740  Mar 25 93  .LINK3D = -1 means no linkage.  J Ewing
;  SPR 10829  Apr  1993  Change info_3d common block
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  nlinks = 0
;
;  Determine the type of the 3-D object.
;  -------------------------------------
  dim1 = object3d(index3d).dim1
  dim2 = object3d(index3d).dim2
  IF(dim1 EQ dim2) THEN BEGIN
    CASE dim1 OF
       32 : type = 'FACE6'
       64 : type = 'FACE7'
      128 : type = 'FACE8'
      256 : type = 'FACE9'
      ELSE: BEGIN
              MESSAGE, 'Unsupported resolution.', /CONT
              RETURN
            END
    ENDCASE
    GOTO, identify_links
  ENDIF
  IF((dim2/3)*4 EQ dim1) THEN BEGIN
    CASE dim1 OF
       128 : type = 'MAP6'
       256 : type = 'MAP7'
       512 : type = 'MAP8'
      1024 : type = 'MAP9'
      ELSE : BEGIN
              MESSAGE, 'Unsupported resolution.', /CONT
              RETURN
            END
    ENDCASE
    GOTO, identify_links
  ENDIF
  IF((dim1 EQ 512) AND (dim2 EQ 256)) THEN BEGIN
    type = 'PROJ_MAP'
    GOTO, identify_links
  ENDIF
  IF((dim1 EQ 1024) AND (dim2 EQ 512)) THEN BEGIN
    type = 'PROJ2_MAP'
    GOTO, identify_links
  ENDIF
  MESSAGE, 'Unsupported 2-D data array size.', /CONT
  RETURN
;
;  See what, if any, 2-D objects are linked to the 3-D object.
;  -----------------------------------------------------------
identify_links:
  j = EXECUTE('siz = SIZE(' + type + ')')
  IF(siz(2) EQ 8) THEN BEGIN
    num_objects = siz(3)
    FOR i = 0, num_objects - 1 DO BEGIN
      j = EXECUTE('link3d = ' + type + '(i).link3d')
      IF(link3d EQ index3d) THEN BEGIN
        name = type + '(' + STRTRIM(STRING(i), 2) + ')'
        IF(nlinks EQ 0) THEN names = [name] ELSE names = [names, name]
        nlinks = nlinks + 1
        IF(KEYWORD_SET(unlink)) THEN j = EXECUTE(type + '(i).link3d = -1')
      ENDIF
    ENDFOR
  ENDIF
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



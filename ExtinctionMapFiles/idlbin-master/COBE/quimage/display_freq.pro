PRO display_freq
;+
;  DIPLAY_FREQ - a UIMAGE-specific routine.  This routine allows the
;  user to select a 3-D object whose associated table of indexes and
;  frequencies will then be printed out.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 06 93  Add more comments.  J Ewing
;  SPR 10456  Jan 13 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 08 93  Change journaling.  J Ewing
;  SPR 10829  Apr  1993  Change info_3d common block
;--------------------------------------------------------------------------
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
   freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON journal,journal_on,luj
  menu_title = 'Display frequency/wavelength table'
;
;  Have the user identify the 3-D object he wants to work with.
;  ------------------------------------------------------------
  select_object3d, index3d
  IF(index3d EQ -1) THEN RETURN
;
;  Print out a table of index's and frequencies.
;  ---------------------------------------------
  freqname = 'FREQ3D_' + STRING(index3d, '(i1)')
  PRINT, ' '
  PRINT, 'Title of 3-D object:  ' + bold(object3d(index3d).title)
  PRINT, underline('    Index    Frequency/Wavelength   ')
  j = EXECUTE('FOR i = 0, object3d(index3d).dim3 - 1 DO PRINT, i + 1, ' + $
    freqname + '(i)')
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand (3-D):  ' + object3d(index3d).title
    PRINTF, luj, '    Index    Frequency/Wavelength   '
    j = EXECUTE('FOR i=0, object3d(index3d).dim3-1 DO PRINTF, luj, ' + $
      'i + 1,' + freqname + '(i)')
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
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



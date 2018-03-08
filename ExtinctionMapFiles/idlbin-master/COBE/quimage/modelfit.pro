PRO modelfit
;+
;  MODELFIT - a UIMAGE-specific routine.  This routine puts up the
;  menu labeled 'Modeling and Fitting', and takes appropriate
;  action depending upon selections made by the user.
;#
;  Written by John Ewing.
;  SPR 10332  Jan 05 93  Put "init = sel" into call to UMENU.  J Ewing
;  SPR 10383  Jan 06 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10870  Apr 29 93  Add call to FITSPEC driver routine. J Newmark
;  SPR 10880  Apr 30 93  Add call to FITMU driver routine. J Newmark
;  SPR 10890  May 04 93  Add call to MKBGRMOD driver. J Newmark
;  SPR 10939  May 11 93  Add call to SVDFIT driver. J Newmark
;  SPR 10941  May 14 93  Add call to MULTIPOLE driver. J Newmark
;--------------------------------------------------------------------------
;
;  Define routines which can be called from here, as well as the
;  phrases, to appear in a menu, which identify those routines.
;  [The VT field below tells whether that routine should be accessible
;  from a VT200 (non-X-window) terminal.]
;  -------------------------------------------------------------------
  mtitle = 'Modeling and Fitting'
  imenu = [mtitle, 'Subtract a Dipole', 'Fit Spectra with a 3 B-Body Model',$
    'Fit Spectra with a Bose-Ein. plus BB Model',$
    '2-D Background Modeling Tool','1-D Polynomial Fitting', $
    'Subtract a Dipole + Quadrupole','1-D Gaussian Fitting', $
    ' ','HELP', 'Return to MAIN MENU']
show_menu:
  sel = UMENU(imenu, title = 0, init = sel)
  CASE sel OF
    1 : udipole
    2 : ubbfit
    3 : ubefit
    4 : ubackmod
    5 : usvdfit
    6 : umulti
    7 : ugfit
    8 : GOTO, show_menu
    9 : BEGIN
          uimage_help, mtitle
          GOTO, show_menu
        END
    10 : RETURN
  ENDCASE
  GOTO, show_menu
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



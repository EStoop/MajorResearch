;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    CLEANPLOT resets IDL system plotting variables to defaults.
;
;DESCRIPTION:
;    CLEANPLOT reinitializes NEARLY ALL IDL system plotting
;    variables back to the default values that existed before
;    the user mucked with them.  See WARNINGS section for a
;    list of system plotting variables which are not reset.
;
;CALLING SEQUENCE:
;    CLEANPLOT
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    There are NO arguments.
;
;WARNINGS:
;  1.  This does NOT reset the plotting device.
;  2.  This does not change any system variables that don't control 
;       plotting.
;  3.  THE FOLLOWING PLOTTING SYSTEM VARIABLES ARE NOT RESET:
;          !P BLOCK:   !p.channel   
;                      !p.clip
;                      !p.t, !p.t3d
;          !X BLOCK:   !x.crange
;                      !x.s
;                      !x.window
;                      !x.tickv, !x.tickname
;          !Y BLOCK:  (see !X block)
;          !Z BLOCK:  (see !X block)
;
;EXAMPLE:
;    To clear all existing plotting labels, axis minima and maxima 
;    limits, tick labelling, etc:
;
;     CLEANPLOT
;
;#
;COMMON BLOCKS:
;     None.
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     Explicitly resets the !C, !P, !X, !Y and !Z
;     system variables to default values.  A small number of
;     these were considered impractical to reset -- see
;     WARNINGS section.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS,ETC.:
;     None.
;
;MODIFICATION HISTORY:
;     Written by JP Weiland, General Sciences Corp. Dec 1991
;         Idea from obsolete Version 1 code from Landsman.    
;
; SPR 9616
;.TITLE
; Routine CLEANPLOT
;-
Pro  CleanPlot, Dummy	   ;Set System Plot Variables to Default Values
;
;;;;;;;;;;;;;;;;;;;;;;;;;;; !C BLOCK ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   !C = 0     ;cursor variable
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;; !P BLOCK ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   !p.background = 0    ; background color
   !p.color = 255
   !p.font = 0
   !p.ticklen = 0     ; can be changed to get a grid
   !p.noclip = 0
   !p.linestyle = 0
   !p.title = ''
   !p.NoErase = 0
   !p.nsum = 0
   !p.psym = 0
   !p.charsize = 0
   !p.charthick = 0
   !p.subtitle = ''
   !p.thick = 0
   !p.position(0:3) = 0
   !p.multi(0:4) = 0
   !p.region(0:3) = 0
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; !X BLOCK ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   !x.range(0:1) = 0
   !X.ticks = 0
   !x.title = ''
   !x.Type = 0
   !x.style = 0
   !x.thick = 0
   !x.charsize = 0
   !x.minor = 0
   !x.margin(0) = 10
   !x.margin(1) = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; !Y BLOCK ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   !y.range(0:1) = 0
   !Y.ticks = 0
   !y.title = ''
   !y.Type = 0
   !y.style = 0
   !y.thick = 0
   !y.charsize = 0
   !y.minor = 0
   !y.margin(0) = 4
   !y.margin(1) = 2
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; !Z BLOCK ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   !z.range(0:1) = 0
   !z.ticks = 0
   !z.title = ''
   !z.Type = 0
   !z.style = 0
   !z.thick = 0
   !z.charsize = 0
   !z.minor = 0
   !z.margin(0) = 0
   !z.margin(1) = 0
;
;
Return					;Completed
End
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



;+
; NAME:
;  EXPORT 
; PURPOSE:
;       This function gets all the output data set information from
;	the user and then does the data set conversion.
; CATEGORY:
;  User interface, Menu, Data I/O
; CALLING SEQUENCE:
;  export
; INPUTS:
;   none
;
; OUTPUTS:
;   none
;
; COMMON BLOCKS:
;  None.
; RESTRICTIONS:
;  None.
; MODIFICATION HISTORY:
;       Creation:  Pete Kryszak, GSC, March 1992.
;
;	16-Dec-1992	10330	Kryszak - Quit functionality
;	23-Dec-1992	10387	Kryszak - added help 
;       31-Mar-1993     10703   Ward - removed /zoomed parameter from
;                                       call to select_object
;       13-Apr-1993     10819   Ward - changed call to dconvert 
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark

;-
;
Pro Export

COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
COMMON uimage_data2,proj_map,proj2_map
COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
  freq3D_2,wt3d_0,wt3d_1,wt3d_2

	menu_title='Select Image to Export'

SelObj:
	name = select_object(/map6, /map7, /map8, /map9,  $
		        /face6,/face7,/face8,/face9, $
		        /object3d, /help,/exit,$
			title=menu_title)

	If name Eq 'NO_OBJECTS' Then Begin
	   message, 'No objects in memory to export.', /continue
	   !err = -1
	   !error = -1
	   Return
	EndIf

	If strupcase(name) Eq 'EXIT' Then goto, back

	If strupcase(name) Eq 'HELP' Then Begin
 	   action, 'export.hlp'
	   Goto, SelObj 
	EndIf

	; so dsname (source) and format ('UIMAGE')
	dsname = name
	inform = 'UIMAGE'

	SelOutpt, output_format, outname
	If !ERR Eq 99 Then goto, back

	If ((output_format eq '') or (outname eq '')) then goto, back

	weights = ''
	if output_format eq 'FITS' then begin
	   weight=0
	   first_letter=STRMID(name,0,1)
	   j=EXECUTE('linkweight= ' +name+ '.linkweight')

	   ; Handle Case where selected object is a 3-dim one.
	   IF (first_letter EQ 'O') THEN BEGIN
	      str_index3d=STRMID(name,9,1)
	      index3d=FIX(str_index3d)
	      j=EXECUTE('havewt= wt3d_' + str_index3d)
	      sz=SIZE(havewt)
	      IF (sz(0) EQ 3) THEN weights = 'yes'
	   ENDIF ELSE IF (linkweight NE -1) THEN weights = 'yes'
	endif

	if weights eq 'yes' then begin
	   print, ' '
	   print, 'The object that you have selected has a weighting array '+$
		  'associated with it.'
	   question = 'Would you like to write it out as well?'
	   ans = one_line_menu([question,'yes','no'],init=1)

nameagain:
	   print, 'Enter the name of the fits file to write the weighting array'
	   print, 'to.  Enter back to return to the main menu.'
	   wname = ''
	   if ans eq 1 then read,'Dataset name? ', wname

	   if strupcase(wname) eq 'BACK' then goto, back

	   if strupcase(wname) eq strupcase(outname) then begin
	      print, 'You must give the file for the weights array a different'
	      print, 'name than you gave for your data arrays'
	      goto, nameagain
	   endif

	   ; set things up for dconvert
	   weights = wname
	endif

	; set those things not needed 
	; for exporting UIMAGE data set
	; dconvert needs the place holders
	;
	intype = ''
	fits_exten = ''
	fnum = ''
	field = ''
	subscr = ''

	dconvert, inform, intype, dsname, fits_exten, fnum, field, subscr, $
          	  weights, output_format, outname
	return

back:
	!err = 1
	!error = 1
	return
	end
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



PRO ucontour
;+
; NAME:	
;	UCONTOUR
; PURPOSE:
;       To plot contours for a given object.
; CATEGORY:
;	
; CALLING SEQUENCE:
;       UCONTOUR
; INPUTS:
;	None.
; OUTPUTS:
;       None.
; COMMON BLOCKS:
;	Uimage_data
; RESTRICTIONS:
;       None.
; PROCEDURE:
;       The SELECT_OBJECT function is called to allow the user to select
;       the desired 'object' for plotting. The user is given the option of
;	entering the levels. The final plot(only contours, no image) maybe
;	sent to a PostScript file.
;
; SUBROUTINES CALLED:
;       Select_object
;	Add_object
; MODIFICATION HISTORY:
;       Creation:  Sarada Chintala, January 1992.
;  SPR 10293  Dec 03,92  Changed to work on zoomed images.  S Chintala
;  SPR 10453  Feb 09,93  Multi-coloured contours can be plotted. S.Chintala
;  SPR 11003  2 Jun 93   Change !version.os to !cgis_os. J Newmark
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;-
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON nonxwindow,term_type
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  COMMON zback,zbgr,zdsrcmap,zbsub

IF(NOT DEFINED(uimage_version)) THEN uimage_version = 2
IF(NOT DEFINED(journal_on)) THEN journal_on = 0
first=1
begin_contour:
ans=' '
res_ans=' '
ps_ans=0
ps_print = 'false'
equally_spaced='false'
;
; Get terminal type for non-Xwindow terminals
;
IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    sz = SIZE(term_type)
    IF(NOT DEFINED(term_type)) THEN BEGIN
       supported_terms=['Terminal type','regis','tek','Other']
       sel=UMENU(supported_terms,title=0)
       term_type=STRUPCASE(supported_terms(sel))
    ENDIF
ENDIF ELSE BEGIN
   zoom_factor = 2^zoom_index
ENDELSE
;
; Select object for generation Contours
;
get_argument:
  menu_title='Contour plots'
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
	/face8,/face9,/proj_map,/proj2_map,/zoomed,/help,$
        /exit,title=menu_title)
  IF(name EQ 'EXIT') THEN BEGIN
    IF((journal_on EQ 1) AND (first NE 1)) THEN $
      PRINTF,luj,'----------------------------------------'+$
                 '--------------------------------------'
    RETURN
  ENDIF ELSE BEGIN
    IF((name NE 'HELP') AND (first NE 1)) THEN BEGIN
       IF(journal_on EQ 1) THEN $
         PRINTF,luj,'----------------------------------------'+$
                 '--------------------------------------'
    ENDIF
    IF(first EQ 1) THEN first = 0
  ENDELSE
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help,menu_title
    GOTO,get_argument
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE,'There are no objects available.',/CONT
    RETURN
  ENDIF
;
; For zoomed objects, extract relevant fields/values from the parent window.
;
  IF(STRMID(name,0,6) EQ 'ZOOMED') THEN BEGIN
    j = EXECUTE('win_orig = '+name+'.win_orig')
    j = EXECUTE('start_x = '+name+'.start_x')
    j = EXECUTE('stop_x = '+name+'.stop_x')
    j = EXECUTE('start_y = '+name+'.start_y')
    j = EXECUTE('stop_y = '+name+'.stop_y')
    bflag=0
    IF (win_orig LT 0) THEN BEGIN
      win_orig=ABS(win_orig)
      bflag=1
    ENDIF
    j=EXECUTE('zoomflag = ' + name + '.zoomflag')
    IF (zoomflag NE 0) THEN bflag=1
    name_orig=get_name(win_orig)
    j = EXECUTE('sz = SIZE('+name_orig+'.data)')
    dim1_orig = sz(1)
    dim2_orig = sz(2)
    x0 = start_x > 0
    x1 = stop_x < (dim1_orig-1)
    y0 = start_y > 0
    y1 = stop_y < (dim2_orig-1)
    j = EXECUTE('vimage = '+name_orig+'.data(x0:x1,y0:y1)')
    istatus = EXECUTE('minimum = '+name_orig+'.data_min')
    istatus = EXECUTE('maximum = '+name_orig+'.data_max')
    IF (bflag EQ 1) THEN BEGIN
      IF (zoomflag EQ 1) THEN vimage=zbgr
      IF (zoomflag EQ 2) THEN vimage=zdsrcmap
      IF (zoomflag EQ 3) THEN vimage=zbsub
      minimum=MIN(vimage)
      maximum=MAX(vimage)
    ENDIF
    j = EXECUTE('specific_zoom = '+name+'.specific_zoom')
    zoom_factor=(2^zoom_index)
    total_zoom = zoom_factor * specific_zoom
    sz = SIZE(vimage)
    IF(start_x GE 0) THEN sx1 = 0 ELSE sx1 = -1*start_x*total_zoom
    sx1 = sx1 + total_zoom / 2
    sx2 = sx1 + (sz(1) - 1) * total_zoom
    IF(start_y GE 0) THEN sy1 = 0 ELSE sy1 = -1*start_y*total_zoom
    sy1 = sy1 + total_zoom / 2
    sy2 = sy1 + (sz(2) - 1) * total_zoom
    istatus = EXECUTE('badpixval='+ name_orig + '.badpixval')
    istatus = EXECUTE('projection='+ name_orig + '.projection')
    j = EXECUTE('window = '+name+'.window')
    WSET,window
  ENDIF ELSE BEGIN
    istatus = EXECUTE('badpixval='+ name + '.badpixval')
    istatus = EXECUTE('vimage='+ name + '.data')
    istatus = EXECUTE('projection='+ name + '.projection')
    istatus = EXECUTE('minimum = '+name+'.data_min')
    istatus = EXECUTE('maximum = '+name+'.data_max')
  ENDELSE
;
; Remove the bad pixel values from the data. The user is given the
; option of re-entering the maximum and minimum values.
;
  z = WHERE(vimage NE 0)
  IF(z(0) EQ -1) THEN BEGIN
     MESSAGE,'This image has no data values.',/CONT
     GOTO, begin_contour
  ENDIF
  IF(minimum EQ maximum) THEN BEGIN
    PRINT,'Image has no dynamic range.'
    IF(uimage_version EQ 2) THEN GOTO,begin_contour ELSE RETURN
  ENDIF
  IF(minimum EQ maximum) THEN BEGIN
    PRINT,'Image consists of all bad pixels.'
    IF(uimage_version EQ 2) THEN GOTO,begin_contour ELSE RETURN
  ENDIF
  IF( badpixval eq 9999999) THEN badpixval = data_min
  w = WHERE(vimage eq badpixval)
  IF( w(0) ne -1) THEN vimage(w) = 1.0E6
  PRINT,'The minimum and maximum values of the image are :  '$
	,strtrim(string(minimum),2),'  ',strtrim(string(maximum),2)
  IF(MIN(vimage) EQ MAX(vimage)) THEN BEGIN
     MESSAGE,'Cannot generate contours for an image in which the maximum'+$
          'and minimum values are the same.',/CONT  
     GOTO, begin_contour
  ENDIF
  IF((STRMID(name,0,4) EQ 'MAP8') AND (projection NE 'SKY-CUBE')) THEN BEGIN
     vimage=vimage(*,0:255)
  ENDIF
; 
; Get the number of contour levels from the user
;
  xlevels = 0
  index = 0
  str=' '
  q = ['Enter the Contour levels?','Yes','No']
  i1 = one_line_menu(q,init=2)
  IF( i1 eq 1) THEN BEGIN
     val=0
     WHILE(val EQ 0) DO BEGIN
        READ,'Enter number of contour levels:  ',str
        val=validnum(str)
	IF(val EQ 1) THEN xlevels=FIX(str) $
		ELSE PRINT,'Invalid number, please re-enter.'
     ENDWHILE
  ENDIF
  IF ( xlevels GT 0) THEN BEGIN
     q = ['Do you want the levels to be equally spaced?','Yes','No']
     i2 = one_line_menu(q,init=1)
     IF(i2 eq 1) THEN BEGIN
        equally_spaced='true'
     ENDIF 
     lstr = FLTARR(xlevels)
  ENDIF
  WHILE((equally_spaced ne 'true') and (xlevels GT index)) DO BEGIN
rd_lev:
     tstr = STRING(index+1,'(i3)')
     val = 0
     WHILE(val EQ 0) DO BEGIN
        READ,'Enter level '+ STRTRIM(tstr,2) +':  ',str
        val = validnum(str)
        IF(val EQ 1) THEN l=FLOAT(str) $
                   ELSE PRINT,'Invalid number, please re-enter.'
     ENDWHILE
     IF(index EQ 0) THEN BEGIN
	lstr(0) = l
     ENDIF ELSE BEGIN
        k = 0
	WHILE ( k LT index) DO BEGIN
  	   IF(l EQ lstr(k)) THEN BEGIN
	      PRINT,'This value has already been entered. Please renter.'
	      GOTO,rd_lev
	   ENDIF
	   k = k+1
	ENDWHILE
	lstr(index) = l
     ENDELSE
     IF((l lt minimum) or (l gt maximum)) THEN BEGIN
        PRINT,'The level should be within the minimum and maximum.',$
		 ' Please reenter'
        GOTO,rd_lev
     ENDIF 
     index = index + 1
  ENDWHILE
  IF((equally_spaced ne 'true')and (xlevels GT 0)) THEN BEGIN
     lstr = lstr(SORT(lstr))
  ENDIF
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
     q=['Would you like multi-coloured contours?','Yes','No']
     ic = one_line_menu(q,init=2)
  ENDIF ELSE BEGIN
     ic=1
  ENDELSE
  IF(ic EQ 1) THEN BEGIN
     colr=[2,3,4,5,6,7,8,9,1]
  ENDIF ELSE BEGIN
     colr=[1]
  ENDELSE
;
; Call the Contour function with the given Max,Min 
;
plot_contour:
  IF(STRMID(name,0,6) EQ 'ZOOMED') THEN i3=2
  IF(((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) and (ps_print ne 'true')$
	AND (STRMID(name,0,6) NE 'ZOOMED')) THEN BEGIN
     q = ['Overlay the contours over another image?',$
		'Yes','No']
     i3 = one_line_menu(q,init=2)
     IF( i3 eq 1) THEN BEGIN
get_argument1:
        menu_title='Select object for overlaying Contours'
        types = strmid(name,0,strpos(name,')')-2)
	j=EXECUTE('overname=select_object(/'+types+',/help,/exit,title=menu_title)')
        IF(overname EQ 'EXIT') THEN RETURN
        IF(overname EQ 'HELP') THEN BEGIN
           uimage_help,menu_title
           GOTO,get_argument1
        ENDIF
        IF(overname EQ 'NO_OBJECTS') THEN BEGIN
           MESSAGE,'There are no objects available.',/CONT
           RETURN
        ENDIF
       istatus = EXECUTE('window='+ overname + '.window')
       WSET,window
       istatus = EXECUTE('xpos='+ overname + '.pos_x')
       istatus = EXECUTE('ypos='+ overname + '.pos_y')
     ENDIF ELSE BEGIN
       istatus = EXECUTE('window='+ name + '.window')
       WSET,window
       istatus = EXECUTE('xpos='+ name + '.pos_x')
       istatus = EXECUTE('ypos='+ name + '.pos_y')
     ENDELSE
     sz = SIZE(vimage)
     sx1 = xpos
     sy1 = ypos
     sx2 = xpos + sz(1)*zoom_factor - 1
     sy2 = ypos + sz(2)*zoom_factor - 1
  ENDIF
  pj_sx1 = 0
  pj_xlevels = 0
  pj_lstr = 0
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
     IF (xlevels GT 0) THEN BEGIN
        IF(equally_spaced eq 'true') THEN BEGIN
             exstr='CONTOUR,vimage,POSITION=[sx1,sy1,sx2,sy2],MAX_VAL=maximum,'+$
		  'NLEVELS=xlevels,/NOERASE,/DEV,XSTYLE=5,YSTYLE=5,'+$
		  '/FOLLOW,C_CHARSIZE=1.25,C_COLOR=colr'
             pj_sx1 = 1
             pj_xlevels = 1
        ENDIF ELSE BEGIN
             exstr='CONTOUR,vimage,POSITION=[sx1,sy1,sx2,sy2],MAX_VAL=maximum,'+$
		  'LEVELS=lstr,/NOERASE,/DEV,XSTYLE=5,YSTYLE=5,'+$
		  '/FOLLOW,C_CHARSIZE=1.25,C_COLOR=colr'
             pj_sx1 = 1
             pj_lstr = 1
	ENDELSE
     ENDIF ELSE BEGIN
        exstr='CONTOUR,vimage,POSITION=[sx1,sy1,sx2,sy2],MAX_VAL=maximum,'+$
		'/NOERASE,/DEV,XSTYLE=5,YSTYLE=5,/FOLLOW,C_CHARSIZE=1.25,'+$
		'C_COLOR=colr'
        pj_sx1 = 1
     ENDELSE
  ENDIF ELSE BEGIN
     IF(xlevels GT 0) THEN BEGIN
        IF(equally_spaced eq 'true') THEN BEGIN
           exstr='CONTOUR,vimage,MAX_VAL=maximum,NLEVELS=xlevels,'+$
			'/FOLLOW,C_CHARSIZE=1.25'
           pj_xlevels = 1
	ENDIF ELSE BEGIN
           exstr='CONTOUR,vimage,MAX_VAL=maximum,LEVELS=lstr,/FOLLOW,'+$
			'C_CHARSIZE=1.25'
           pj_lstr = 1
	ENDELSE
     ENDIF ELSE BEGIN
        exstr='CONTOUR,vimage,MAX_VAL=maximum,/FOLLOW,C_CHARSIZE=1.25'
     ENDELSE
  ENDELSE

  j=EXECUTE(exstr)

  IF((journal_on EQ 1) and (ps_print NE 'true')) THEN BEGIN
    PRINTF,luj,menu_title
    istatus = EXECUTE('title='+ name + '.title')
    PRINTF,luj,'  operand:  '+title
    PRINTF,luj,'  the following statement was executed:'
    WHILE(STRLEN(exstr) GT 76) DO BEGIN
       jstr=STRMID(exstr,0,76)
       PRINTF,luj,'  '+jstr
       exstr=STRMID(exstr,76,STRLEN(exstr))
    ENDWHILE
    PRINTF,luj,'  '+exstr
    PRINTF,luj,'  where:'
    PRINTF,luj,'    image --> '+title
    PRINTF,luj,'    maximum = '+STRTRIM(STRING(maximum),2)
    IF(pj_sx1) THEN BEGIN
      PRINTF,luj,'    sx1 = ' + STRTRIM(STRING(sx1),2)
      PRINTF,luj,'    sx2 = ' + STRTRIM(STRING(sx2),2)
      PRINTF,luj,'    sy1 = ' + STRTRIM(STRING(sy1),2)
      PRINTF,luj,'    sy2 = ' + STRTRIM(STRING(sy2),2)
    ENDIF
    IF(pj_xlevels) THEN PRINTF,luj,'    xlevels = ' + STRTRIM(STRING(xlevels),2)
    IF(pj_lstr) THEN PRINTF,luj,'    lstr = ' , STRTRIM(STRING(lstr), 2)
    IF(i1 EQ 1) THEN BEGIN
      PRINTF,luj,'  the number of contours selected is :  '+$
	   STRTRIM(STRING(xlevels),2)
      IF(i2 eq 1) THEN BEGIN
        PRINTF,luj,'  the contours selected are equally spaced'
      ENDIF ELSE BEGIN
         PRINTF,luj,'  the contours selected are not equally spaced'
      ENDELSE
    ENDIF ELSE BEGIN
      PRINTF,luj,'  the default number of contours has been selected'
    ENDELSE
    IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
      IF(i3 EQ 1) THEN BEGIN
         istatus = EXECUTE('overtitle='+ overname + '.title')
         PRINTF,luj,'  the contours have been overlayed over image '+overtitle
      ENDIF ELSE BEGIN
         PRINTF,luj,'  the contours have been overlayed on the same image 
      ENDELSE
    ENDIF
  ENDIF

;
; After printing the PostScript plot, call set_plot with what it
; was previously set to before calling set_plot,'ps'. 
;
  IF ( ps_print EQ 'true') THEN BEGIN
     PRINT,'The PS output was sent to '+filename+'.'
     IF(journal_on EQ 1) THEN BEGIN
        PRINTF,luj,'  the PostScript output was sent to '+filename
         PRINTF,luj,'----------------------------------------'+$
                 '--------------------------------------'
     ENDIF
     DEVICE,/CLOSE
     SET_PLOT,term
     set_viewport,0.150005,0.955005,0.110005,0.945005
     RETURN
  ENDIF
;
; In case of non-X terminals, if the user hits any key after the plot is
; drawn, then he/she goes back to the screen where the IDL is running.
; Then the user is asked the question related to the printing of PostScript 
; output. This following code is to prevent overwriting the plot with question.
;
  q = ['Do you want a PostScript output?','Yes','No']
  IF ((!D.NAME ne 'X') AND (!D.NAME NE 'WIN')) then BEGIN
    CASE term_type OF
       'TEK'     : BEGIN
		          i=get_kbrd(1)&text
		       END      
       'REGIS' : BEGIN
		          i=get_kbrd(1)&erase
		      END
    ENDCASE
    print,' '
  ENDIF
  ps_ans=one_line_menu(q,init=2)
  IF(journal_on EQ 1) THEN BEGIN
     IF(ps_ans EQ 1) THEN BEGIN
        PRINTF,luj,'  the PostScript option has been selected for contours'
     ENDIF ELSE BEGIN
        PRINTF,luj,'  no PostScript file was produced' 
     ENDELSE
  ENDIF
;
; If PostScript option has been specified, then save the previous device
; name used in SET_PLOT in variable term. 
;
  IF ( ps_ans eq 1) THEN BEGIN
     print,'The PostScript file will ONLY contain the contours (not the image also).'  
     IF(!cgis_os EQ 'vms') THEN BEGIN
        print,'The PostScript output will by default be sent to '+bold('SYS$LOGIN:IDL.PS')
     ENDIF ELSE BEGIN 
           IF (!cgis_os EQ 'unix') THEN $
        print,'The PostScript output will by default be sent to $HOME/idl.ps.'$
                ELSE IF (!cgis_os EQ 'windows') THEN $
          print,'The PostScript output will by default be sent to HOME\idl.ps.'
         ENDELSE
     menu = ['Output file name:','Use default','Specify a different name']
     isel = one_line_menu(menu,init=1)
     IF(isel EQ 2) THEN BEGIN
        read_file=0
        file=''
        dr=''
        WHILE(read_file EQ 0) DO BEGIN
           read,'Enter directory name:  ',dr
           read,'Enter filename:  ',file
	   filename=dr+file
	   OPENW,unit,filename,/GET_LUN,/DELETE,error=err
	   IF(err NE 0) THEN PRINT,!ERR_STRING ELSE read_file=1
        ENDWHILE
        CLOSE,unit
     ENDIF ELSE BEGIN
        IF(!cgis_os EQ 'vms') THEN BEGIN
	   filename='SYS$LOGIN:IDL.PS'
        ENDIF ELSE BEGIN
           IF(!cgis_os EQ 'unix') THEN filename='$HOME/idl.ps' $
           ELSE IF (!cgis_os EQ 'windows') THEN filename='HOME\idl.ps'
        ENDELSE
     ENDELSE
     term=!D.name
     SET_PLOT,'PS'
     SET_VIEWPORT,0.1,0.9,0.1,0.8
     IF(ic EQ 1) THEN BEGIN
        DEVICE,/color,filename=filename
     ENDIF ELSE BEGIN
        DEVICE,filename=filename
     ENDELSE
     ypos = 0.93
     ps_print='true'
     GOTO,plot_contour
  ENDIF
;
; Before returning, call set_viewport with system default sizes for window views.
;
  set_viewport,0.150005,0.955005,0.110005,0.945005
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
     menu = ['Remove contours from image?','Yes','No']
     sel = one_line_menu(menu,init=2)
     IF(i3 NE 1) THEN BEGIN
       IF(sel EQ 1) THEN IF(STRMID(name,0,6) NE 'ZOOMED') $
          THEN xdisplay,name,/rescale ELSE show_zoomed,name
     ENDIF ELSE BEGIN
       IF(sel EQ 1) THEN IF(STRMID(overname,0,6) NE 'ZOOMED') $
          THEN xdisplay,overname,/rescale ELSE show_zoomed,overname
     ENDELSE
  ENDIF

  IF(uimage_version EQ 2) THEN GOTO, begin_contour ELSE RETURN
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



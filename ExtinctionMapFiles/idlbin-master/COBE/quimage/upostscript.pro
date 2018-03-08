]PRO upostscript
;+
; NAME:	
;	UPOSTSCRIPT
; PURPOSE:
;       To create a PostScript for a given object.
; CATEGORY:
;	
; CALLING SEQUENCE:
;       UPOSTSCRIPT
; INPUTS:
;	None.
; OUTPUTS:
;       None.
; COMMON BLOCKS:
;	Uimage_data
;	uimage_data2
;	xwindow
;	journal
;	history
;	colors
;	color_value
;	pscolor
; RESTRICTIONS:
;       None.
; PROCEDURE:
;	The UPOSTSCRIPT routine prompts the user for the following options:
;	Color or No color?
;	Screen dump? (with overlays)
;	Portrait or Landscape
;	At the end, the user is given the option to specify where the output
;	file should be written. 
;
; SUBROUTINES CALLED:
;       Select_object
; MODIFICATION HISTORY:
;       Creation:  Sarada Chintala, January 1992.
;  SPR 10087  Dec 10,92  Change to add title for Landscape.  S Chintala
;  SPR 10625  Mar 30,93  Change to switch from B&W to color. S.Chintala
;  SPR 10648  Mar 30,93  Change to user interface. S.Chintala
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
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  COMMON colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
  COMMON color_values,c_badpix,c_draw,c_scalemin
  COMMON zback,zbgr,zdsrcmap,zbsub

begin_postscript:
  IF(NOT DEFINED(uimage_version)) THEN uimage_version = 2
  IF(NOT DEFINED(journal_on)) THEN journal_on = 0
  first=1
get_argument:
  menu_title='Create a PostScript file'
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
	/face8,/face9,/proj_map,/proj2_map,/graph,/zoomed,$
        /help,/exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help,menu_title
    GOTO,get_argument
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE,'There are no objects available.',/CONT
    RETURN
  ENDIF
  IF(journal_on EQ 1) THEN BEGIN
    PRINTF,luj,menu_title
    istatus = EXECUTE('xtitle='+ name + '.title')
    PRINTF,luj,'  operand:  '+xtitle
    first = 0
  ENDIF
  IF(!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN') THEN BEGIN
     zoom_factor = 2^zoom_index
  ENDIF
  first_letter = STRMID(name,0,1)
  term = !D.NAME
;;;
;;; Extract various fields for zoomed objects
;;;
  IF(first_letter EQ 'Z') THEN BEGIN
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
     j = EXECUTE('data = '+name_orig+'.data(x0:x1,y0:y1)')
     IF (bflag EQ 1) THEN BEGIN
       IF (zoomflag EQ 1) THEN data=zbgr
       IF (zoomflag EQ 2) THEN data=zdsrcmap
       IF (zoomflag EQ 3) THEN data=zbsub
     ENDIF
     j = EXECUTE('badpixval='+ name_orig + '.badpixval')
     j = EXECUTE('projection='+ name_orig + '.projection')
     j = EXECUTE('scale_min='+ name + '.scale_min')
     j = EXECUTE('scale_max='+ name + '.scale_max')
     j = EXECUTE('title_orig= '+name_orig+'.title')
     j = EXECUTE('specific_zoom = '+name+'.specific_zoom')
     zoom_factor=(2^zoom_index)*specific_zoom
     specific_zoom = STRTRIM(STRING(specific_zoom),2)
     start_x   = STRTRIM(STRING(start_x),2)
     start_y   = STRTRIM(STRING(start_y),2)
     stop_x    = STRTRIM(STRING(stop_x),2)
     stop_y    = STRTRIM(STRING(stop_y),2)
     vimage = BYTSCL(data,scale_min,scale_max)
  ENDIF ELSE BEGIN
;;;
;;; Extract fields for non-Graph objects like Map,Face and Projections
;;;
     IF(first_letter NE 'G') THEN BEGIN
    	istatus = EXECUTE('data='+ name + '.data')
        j = EXECUTE('scale_min = '+name+'.scale_min')
        j = EXECUTE('scale_max = '+name+'.scale_max')
        j = EXECUTE('badpixval = '+name+'.badpixval')
        types = strmid(name,0,strpos(name,')')-2)
        CASE types OF
           'MAP6'     : BEGIN
			   data(0:95,0:31)=scale_max
			   data(0:95,64:95)=scale_max
		        END      
           'MAP7'     : BEGIN
			   data(0:191,0:63)=scale_max
			   data(0:191,128:191)=scale_max
		        END
           'MAP8'     : BEGIN
			   data(0:383,0:127)=scale_max
			   data(0:383,256:383)=scale_max
		        END
           'MAP9'     : BEGIN
			   data(0:767,0:255)=scale_max
			   data(0:767,512:767)=scale_max
		        END
	    ELSE      :
        ENDCASE
    ENDIF
  ENDELSE
;;;
;;; Handle overlays ( like contours)
;;; If the contours are multicoloured, the colors donot show on the PS output.
;;; The color table on the workstation vary from 0 to 244, whereas the colors
;;; on the color printer vary from 0 to 20. An internal Congrid is applied to
;;; rescale the colors - therefore they donot appear very distinct.
;;;
  overlay='false'
  IF(((term EQ 'X') OR (term EQ 'WIN')) AND ( first_letter NE 'G')) THEN BEGIN
     print,' '
     print,'Select one of the two display methods shown below:'
     menu = ['','Data alone(no overlays)','Data+overlays(screen capture)']
     sel = one_line_menu(menu,init=1)
     IF(journal_on EQ 1) THEN BEGIN
       IF(sel EQ 1) THEN BEGIN
          PRINTF,luj,'  Data alone(no overlays)'
       ENDIF ELSE BEGIN
          PRINTF,luj,'  A PostScript output of a screen dump of the image has been selected'
       ENDELSE
     ENDIF
     IF(sel EQ 2) THEN BEGIN
     	istatus = EXECUTE('window='+ name + '.window')
        WSET,window
	vimage=TVRD(0,0,!D.X_SIZE,!D.Y_SIZE)
	overlay='true'
     ENDIF 
     tvlct,rcolor,gcolor,bcolor,/get
  ENDIF ELSE BEGIN
     IF(journal_on EQ 1) THEN BEGIN
        PRINTF,luj,'  Data alone(no overlays)'
     ENDIF
  ENDELSE

create_ps:

;
;  Extract fields for an element of MAP* or FACE*
;  ------------------------------------------------
  IF((first_letter EQ 'M')$
     OR (first_letter EQ 'F') OR (first_letter EQ 'P'))THEN BEGIN
     j = EXECUTE('siz = SIZE('+name+'.data)')
     dim1 = siz(1)
     dim2 = siz(2)
     j = EXECUTE('bandno    = '+name+'.bandno')
     j = EXECUTE('bandwidth = '+name+'.bandwidth')
     j = EXECUTE('coordinate_system = '+name+'.coordinate_system')
     j = EXECUTE('faceno    = '+name+'.faceno')
     j = EXECUTE('frequency = '+name+'.frequency')
     j = EXECUTE('instrume  = '+name+'.instrume')
     j = EXECUTE('orient    = '+name+'.orient')
     j = EXECUTE('projection= '+name+'.projection')
     j = EXECUTE('units     = '+name+'.units')
     j = EXECUTE('version   = '+name+'.version')
     j = EXECUTE('good = WHERE('+name+'.data NE badpixval)')
     j = EXECUTE('data_min = '+name+'.data_min')
     j = EXECUTE('data_max = '+name+'.data_max')
     IF(data_min EQ data_max) THEN BEGIN
        PRINT,'Image has no dynamic range.'
        IF(uimage_version EQ 2) THEN GOTO,begin_postscript ELSE RETURN
     ENDIF
     IF(data_min GT data_max) THEN BEGIN
        PRINT,'Image consists of all bad pixels.'
        IF(uimage_version EQ 2) THEN GOTO,begin_postscript ELSE RETURN
     ENDIF
     CASE first_letter OF
      	'M' : res = STRMID(name,3,1)
      	'F' : res = STRMID(name,4,1)
        ELSE:
     ENDCASE
     dim1 = STRTRIM(STRING(dim1),2)
     dim2 = STRTRIM(STRING(dim2),2)
     badpixval = STRTRIM(STRING(badpixval),2)
     bandno =    STRTRIM(STRING(bandno),2)
     bandwidth = STRTRIM(STRING(bandwidth),2)
     data_min  = STRTRIM(STRING(data_min),2)
     data_max  = STRTRIM(STRING(data_max),2)
     faceno    = STRTRIM(STRING(faceno),2)
     frequency = STRTRIM(STRING(frequency),2)
     scale_min = STRTRIM(STRING(scale_min),2)
     scale_max = STRTRIM(STRING(scale_max),2)
  ENDIF
;
; User chooses B/W or color
;
  PRINT,'You will be prompted below to supply some information regarding the'
  PRINT,'format of the PostScript file, as well as where it should be stored.'
  menu = ['Color or Black&White?','Color','Black&White']
  IF(NOT DEFINED(csel)) THEN csel = 2
  csel = one_line_menu(menu,init=2)
  IF(journal_on EQ 1) THEN BEGIN
     IF( csel EQ 1) THEN BEGIN
        PRINTF,luj,'  Color'
     ENDIF ELSE BEGIN
        PRINTF,luj,'  Black & White'
     ENDELSE
  ENDIF
;
; Portrait or Landscape?
;
  menu = ['Portrait or Landscape?','Portrait','Landscape']
  sel = one_line_menu(menu,init=1)
  IF(journal_on EQ 1) THEN BEGIN
     IF(sel EQ 1) THEN BEGIN
        PRINTF,luj,'  Portrait'
     ENDIF ELSE BEGIN
        PRINTF,luj,'  Landscape'
     ENDELSE
  ENDIF
;
; The PS file will be written be to appropriate directories on VMS and ULTRIX
;
  IF(!cgis_os EQ 'vms') THEN BEGIN
     print,'The PostScript output will by default be sent to '+$
	       bold('SYS$LOGIN:IDL.PS')
  ENDIF ELSE BEGIN
     IF(!cgis_os EQ 'unix') THEN $
        print,'The PostScript output will by default be sent to $HOME/IDL.PS' $
     ELSE IF (!cgis_os EQ 'windows') THEN $
        print,'The PostScript output will by default be sent to HOME\IDL.PS' 
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
        ELSE IF (!cgis_os EQ 'windows') THEN BEGIN
                home=getenv('HOME')
                filename=home+'idl.ps'
              ENDIF
    ENDELSE
  ENDELSE
  IF(journal_on EQ 1) THEN BEGIN
     PRINTF,luj,'  the PostScript output was sent to '+filename
  ENDIF
;
; Save the current graphic output device; Set the output device to PS using
; the SET_PLOT procedure.
;
  j = EXECUTE('title = '+name+'.title')
  old_dname=!D.NAME
  SET_PLOT,'PS'
;
; For a non-Xwindow user, a color table is provided for choice of color.
;
  IF((csel EQ 1) AND (old_dname NE 'X') AND (old_dname NE 'WIN')) THEN BEGIN    
     menu_title='Standard IDL color table'
     menu_items = [menu_title,$
                'B-W Linear',           $  ; 0
                'Blue/White',           $  ; 1
                'Grn-Red-Blu-Wht',      $  ; 2
                'Red Temperature',      $  ; 3
                'Blue/Green/Red/Yellow',$  ; 4
                'Std Gamma-II',         $  ; 5
                'Prism',                $  ; 6
                'Red-Purple',           $  ; 7
                'Green/White Linear',   $  ; 8
                'Green/Wht Exponential',$  ; 9
                'Green-Pink',           $  ;10
                'Blue-Red',             $  ;11
                '16 Level',             $  ;12
                '16 Level II',          $  ;13
                'Steps',                $  ;14
                'Stern Special',        $  ;15
                '',$
                'HELP',   $
                'Return to previous menu']
show_menu:
     menu_select = umenu(menu_items,TITLE=0,INIT=1)
     IF(menu_select EQ 17) THEN GOTO,show_menu
     IF(menu_select EQ 18) THEN BEGIN
       uimage_help,menu_title
       GOTO,show_menu
     ENDIF
     IF(menu_select EQ 19) THEN BEGIN
       IF((journal_on EQ 1) AND (first NE 1)) THEN $
         PRINTF,luj,'----------------------------------------'+$
                 '--------------------------------------'
       RETURN
     ENDIF
     LOADCT,menu_select-1
  ENDIF ELSE IF (csel EQ 1) THEN tvlct,rcolor,gcolor,bcolor
;
; The image is scaled using BYTSCL
;
  IF((first_letter NE 'G') AND (overlay NE 'true')) THEN BEGIN
     bad = WHERE(data EQ badpixval)
     IF((old_dname EQ 'X') OR (old_dname EQ 'WIN')) THEN BEGIN
	 c_scalemax=N_ELEMENTS(r)-1
         vimage = BYTSCL(data,scale_min,scale_max,TOP=(c_scalemax-c_scalemin))+$
		  c_scalemin
     ENDIF ELSE BEGIN
         vimage = BYTSCL(data,scale_min,scale_max)
     ENDELSE
     IF(bad(0) NE -1) THEN vimage(bad)=0
  ENDIF
;
; Text describing the image is written onto the PS file for Portrait plots.
; The Landscape plots have the title and the color table.
;
  CASE sel OF
      1 : BEGIN

	      IF(first_letter EQ 'P') THEN BEGIN
		 IF(csel EQ 1) THEN BEGIN
	            DEVICE,/portrait,/color,xoff=2.54,yoff=15,xsiz=16,ysiz=8,$
			   filename=filename
		 ENDIF ELSE BEGIN
	            DEVICE,/portrait,xoff=2.54,yoff=15,xsiz=16,ysiz=8,$
			filename=filename,color=0
		 ENDELSE
	      ENDIF ELSE BEGIN
	         IF(csel EQ 1) THEN BEGIN
	            DEVICE,/portrait,/color,xoff=2.54,xsiz=16.51,ysiz=12.3825,$
			filename=filename
	         ENDIF ELSE BEGIN
		    IF(first_letter NE 'G') THEN BEGIN
	               DEVICE,/portrait,xoff=2.54,xsiz=16.51,ysiz=12.3825,$
			  filename=filename,color=0
		    ENDIF ELSE BEGIN
	               DEVICE,/portrait,xoff=2.54,yoff=10,xsiz=16.51,$
			  ysiz=12.3825,filename=filename,color=0
		    ENDELSE
	         ENDELSE
	      ENDELSE
	      xpos=0
	      IF(first_letter NE 'G') THEN BEGIN
		 IF(first_letter EQ 'P') THEN BEGIN
	            ypos=-2000
	            xpos1=4500
                    XYOUTS,xpos,9000,title,charthick=2,font=6,size=1.5,/device
		 ENDIF ELSE BEGIN
	            ypos=-1500
	            xpos1=4500
                    XYOUTS,xpos,13500,title,charthick=2,font=6,size=1.5,/device
		 ENDELSE
	      ENDIF
  	      IF(first_letter EQ 'Z') THEN BEGIN
		 ypos=ypos-1500
		 xpos1=xpos1-1500
		 gap=700
                 XYOUTS,xpos,ypos,'!6Data origin',size=1.1,/device
		 XYOUTS,xpos1,ypos,'= '+title_orig,size=1.1,/device
		 ypos=ypos-gap
                 XYOUTS,xpos,ypos,'!6X range',size=1.1,/device
		 XYOUTS,xpos1,ypos,'= '+start_x+' to '+stop_x,size=1.1,/device
		 ypos=ypos-gap
                 XYOUTS,xpos,ypos,'!6Y range',size=1.1,/device
		 XYOUTS,xpos1,ypos,'= '+start_y+' to '+stop_y,size=1.1,/device
		 ypos=ypos-gap
		 XYOUTS,xpos,ypos,'!6Zoom factor',size=1.1,/device
		 XYOUTS,xpos1,ypos,'= '+specific_zoom,size=1.1,/device
  	      ENDIF
  	      IF((first_letter EQ 'M')$
	 	 OR (first_letter EQ 'F') OR (first_letter EQ 'P'))THEN BEGIN
;
;  Print out fields for an element of MAP* or FACE*
;  ------------------------------------------------
		  IF(first_letter EQ 'P') THEN gap=700 ELSE gap=600
                  XYOUTS,xpos,ypos,'!6Array dimensions',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+dim1+' X '+dim2,size=1.1,/device
		  ypos=ypos-gap
                  XYOUTS,xpos,ypos,'!6Bad pixel value',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+badpixval,size=1.1,/device
		  ypos=ypos-gap
                  XYOUTS,xpos,ypos,'!6Band number',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+bandno,size=1.1,/device
		  ypos=ypos-gap
                  XYOUTS,xpos,ypos,'!6Bandwidth',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+bandwidth,size=1.1,/device
		  ypos=ypos-gap
                  XYOUTS,xpos,ypos,'!6Coordinate-system',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+coordinate_system,size=1.1,/device
		  ypos=ypos-gap
                  XYOUTS,xpos,ypos,'!6Data minimum',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+data_min,size=1.1,/device
		  ypos=ypos-gap
                  XYOUTS,xpos,ypos,'!6Data maximum',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+data_max,size=1.1,/device
		  ypos=ypos-gap
    		  IF(first_letter EQ 'F') THEN BEGIN
                     XYOUTS,xpos,ypos,'!6Face number',size=1.1,/device
		     XYOUTS,xpos1,ypos,'= '+faceno,size=1.1,/device
		     ypos=ypos-gap
		  ENDIF
                  XYOUTS,xpos,ypos,'!6Frequency',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+frequency,size=1.1,/device
		  ypos=ypos-gap
                  XYOUTS,xpos,ypos,'!6Instrument',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+instrume,size=1.1,/device
		  ypos=ypos-gap
                  XYOUTS,xpos,ypos,'!6Orientation',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+orient,size=1.1,/device
		  ypos=ypos-gap
                  XYOUTS,xpos,ypos,'!6Projection',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+projection,size=1.1,/device
		  ypos=ypos-gap
    		  IF(first_letter NE 'P') THEN BEGIN
                     XYOUTS,xpos,ypos,'!6Resolution',size=1.1,/device
		     XYOUTS,xpos1,ypos,'= '+res,size=1.1,/device
		     ypos=ypos-gap
		  ENDIF
                  XYOUTS,xpos,ypos,'!6Scaling minimum',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+scale_min,size=1.1,/device
		  ypos=ypos-gap
                  XYOUTS,xpos,ypos,'!6Scaling maximum',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+scale_max,size=1.1,/device
		  ypos=ypos-gap
                  XYOUTS,xpos,ypos,'!6Units',size=1.1,/device
		  XYOUTS,xpos1,ypos,'= '+units,size=1.1,/device
  	      ENDIF
  	      IF(first_letter NE 'G') THEN BEGIN
		 TV,vimage
	      ENDIF ELSE BEGIN
 		plot_graph,name
	        xpos=3700
                XYOUTS,xpos,14000,title,charthick=1.25,font=6,size=1.5,/device
	      ENDELSE
	   END
       2 : BEGIN
	      IF(first_letter EQ 'P') THEN BEGIN
		 ypos=13000
		 xpos=12000
		 IF(csel EQ 1) THEN BEGIN
	            DEVICE,/landscape,/color,xoff=5,yoff=25,xsiz=22,ysiz=11,$
			filename=filename
		 ENDIF ELSE BEGIN
	            DEVICE,/landscape,xoff=5,yoff=25,xsiz=22,ysiz=11,$
			filename=filename,color=0
		 ENDELSE
	      ENDIF ELSE BEGIN
		 IF(csel EQ 1) THEN BEGIN
	            DEVICE,/landscape,/color,filename=filename,xoff=3,$
                        yoff=25,xsiz=22,ysize=16.5
		 ENDIF ELSE BEGIN
	            DEVICE,/landscape,filename=filename,xoff=3,yoff=25,$
			xsize=22,ysize=16.5,color=0
		 ENDELSE
	      ENDELSE
  	      IF(first_letter NE 'G') THEN BEGIN
		 IF((first_letter EQ 'M') AND (csel EQ 1) AND $
                    (overlay NE 'true')) THEN BEGIN
		    temp=CONGRID(data,1024,768)
        	    bad = WHERE(temp EQ badpixval)
		    vimage=BYTSCL(temp,scale_min,scale_max)
        	    IF(bad(0) NE -1) THEN vimage(bad)=0
		    a=indgen(256)
		    cb=intarr(260,40)
		    for i=2,37 do cb(2,i)=a
		    vimage(210,70)=cb
		    TV,vimage
		    str=' '
		    len=0
		    WHILE(len EQ 0) DO BEGIN 
		       READ,'Enter text for color table title :  ',str
		       len=STRLEN(str)
		    ENDWHILE
		    siz=30./len
		    IF(len GT 80) THEN BEGIN
		       str=STRMID(str,0,80)
		       len=80
		       siz=0.9
		    ENDIF ELSE BEGIN
		    	siz=(siz > 1) < 1.3
		    ENDELSE
         	    XYOUTS,8000,3000,str,charthick=1.25,font=6,size=siz,$
			  /device,alignment=0.5
		 ENDIF ELSE BEGIN
		    TV,vimage
		 ENDELSE
		 len=STRLEN(title)
		 IF (first_letter EQ 'P') THEN ypos=14000 ELSE $
                   IF (first_letter EQ 'Z') THEN ypos=18100 $
                   ELSE ypos=17000
		 IF(len GT 60) THEN BEGIN
		    xpos=0
		    IF(len GT 80) THEN BEGIN
		       title=STRMID(title,0,80)
		       len=80
		    ENDIF
		    siz=60./len
		    siz=0.9
		 ENDIF ELSE BEGIN
		    xpos=900
		    siz=60./len
		    siz=(siz > 0.9) < 2.4
		 ENDELSE
         	 XYOUTS,xpos,ypos,title,charthick=1.25,font=6,size=siz,$
			/device
	      ENDIF ELSE BEGIN
 		 plot_graph,name
	      ENDELSE
	   END
   ENDCASE
   IF(isel EQ 2) THEN print,'The PostScript output has been sent to '+bold(filename)
   DEVICE,/close	
   set_plot,term
   set_viewport,0.150005,0.955005,0.110005,0.945005
   IF(journal_on EQ 1) THEN $
      PRINTF,luj,'----------------------------------------'+$
                 '--------------------------------------'
  IF(uimage_version EQ 2) THEN GOTO, get_argument ELSE RETURN
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



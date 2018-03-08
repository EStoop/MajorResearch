PRO uhistogram,name
;+
; NAME:	
;	UHISTOGRAM
; PURPOSE:
;       To plot a histogram for a given object.
; CATEGORY:
;	
; CALLING SEQUENCE:
;       UHISTOGRAM
; INPUTS:
;	None.
; OUTPUTS:
;       None.
; COMMON BLOCKS:
;	Uimage_data
; 	uimage_data2
;	nonxwindow
;	journal
;	history
;
; RESTRICTIONS:
;       None.
; PROCEDURE:
;       The SELECT_OBJECT function is called to allow the user to select
;       the desired 'object' for plotting. The user has to enter the
; 	following: 1) bin size 2) Linear/Logarithmic plot. He/she also
;	has the option of altering the Minimum and Maximum for the plot.
;	At the end the plot may be resized if the user selects that option.
;
; SUBROUTINES CALLED:
;       Select_object
;	Add_object
; MODIFICATION HISTORY:
;       Creation:  Sarada Chintala, January 1992.
; spr 10389, called append_number to add a number to the title so that it
;	     may be unique.
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;  SPR 11174  Jun 14 94  Upgrade for use of PXINPOLY, J. Newmark
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;-

COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
COMMON uimage_data2,proj_map,proj2_map
COMMON nonxwindow,term_type
COMMON journal,journal_on,luj
COMMON history,uimage_version
COMMON zback,zbgr,zdsrcmap,zbsub

IF(NOT DEFINED(journal_on)) THEN journal_on = 0
sz = SIZE(uimage_version)
IF(NOT DEFINED(uimage_version)) THEN uimage_version = 2
first_resize = 1
nblocks_x=5
nblocks_y=5
begin_histogram:
ans=' '
res_ans=' '
str=' '
units=' '
exstr=' '
area_selected='false'
;
; Get terminal type for non-Xwindow terminals
;
IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    IF(NOT DEFINED(term_type)) THEN BEGIN
       supported_terms=['Terminal type','Regis','Tek','Other']
       sel=UMENU(supported_terms,title=0)
       term_type=STRUPCASE(supported_terms(sel))
    ENDIF
ENDIF
;
; Select object for plotting Histogram
;
get_argument:
  menu_title='Statistics & Histogram'
  IF(N_PARAMS(0) EQ 1) THEN GOTO,have_name
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
	/face8,/face9,/proj_map,/proj2_map,/zoomed,/help,$
        /exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help,menu_title
    GOTO,get_argument
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE,'There are no objects available.',/CONT
    RETURN
  ENDIF
have_name:
;
; For zoomed objects, the fields/values are extracted from the parent image
;
  istatus = EXECUTE('title='+ name + '.title')
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
    istatus = EXECUTE('badpixval='+ name_orig + '.badpixval')
  ENDIF ELSE BEGIN
    istatus = EXECUTE('badpixval='+ name + '.badpixval')
  ENDELSE
;
; Remove the bad pixel values from the data. The user is given the
; option of re-entering the maximum and minimum values.
;
  PRINT,'Do you want to work with the whole image or a section of the image?'
  q = ['Select Mode:','Whole image','A section']
  sel = one_line_menu(q,init=1)
  IF ( sel eq 2) THEN BEGIN
again:
    area_selected = 'true'
    area=define_area(name)
    good = WHERE(area NE badpixval)
    IF(good(0) EQ -1) THEN BEGIN
       MESSAGE,'The selected area consists of all bad pixels.',/CONT
       IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
          menu = ['Remove area boundary from image?','Yes','No']
          sel = one_line_menu(menu,init=1)
          IF(sel EQ 1) THEN BEGIN
	     show_zoomed,name
	     area_selected = 'false'
	  ENDIF
       ENDIF
       menu = ['Select another area on this image?','Yes','No']
       sel = one_line_menu(menu,init=2)
       IF(sel EQ 1) THEN GOTO,again ELSE return
    ENDIF ELSE BEGIN
       vimage = area
       title= append_number('Histogram of user defined area')
    ENDELSE
  ENDIF ELSE BEGIN
     title = 'Histogram of '+'"'+title+'"'
     IF (strlen(title) gt 40) THEN title = append_number('Histogram plot')
     IF(STRMID(name,0,6) EQ 'ZOOMED') THEN BEGIN
        j = EXECUTE('sz = SIZE('+name_orig+'.data)')
        dim1_orig = sz(1)
        dim2_orig = sz(2)
        x0 = start_x > 0
        x1 = stop_x < (dim1_orig-1)
        y0 = start_y > 0
        y1 = stop_y < (dim2_orig-1)
        j = EXECUTE('vimage = '+name_orig+'.data(x0:x1,y0:y1)')
        IF (bflag EQ 1) THEN BEGIN
           IF (zoomflag EQ 1) THEN vimage=zbgr
           IF (zoomflag EQ 2) THEN vimage=zdsrcmap
           IF (zoomflag EQ 3) THEN vimage=zbsub
        ENDIF
     ENDIF ELSE BEGIN
        istatus = EXECUTE('vimage='+ name + '.data')
     ENDELSE
     IF(STRMID(name,0,4) EQ 'PROJ') THEN BEGIN
       j = EXECUTE('win_orig = '+name+'.win_orig')
       name_orig=get_name(win_orig)
       istatus = EXECUTE('vimage='+ name_orig + '.data')
     ENDIF
  ENDELSE
  IF((STRMID(name,0,6) EQ 'ZOOMED') OR (area_selected EQ 'true') $
     OR (STRMID(name,0,4) EQ 'PROJ')) THEN BEGIN
     minimum=MIN(vimage,MAX=maximum)
  ENDIF ELSE BEGIN
     j = EXECUTE('minimum = '+name+'.data_min')
     j = EXECUTE('maximum = '+name+'.data_max')
  ENDELSE
  IF(minimum GT maximum) THEN BEGIN
    PRINT,'Image consists of all bad pixels.'
    IF(uimage_version EQ 2) THEN GOTO,begin_histogram ELSE RETURN
  ENDIF
  IF( badpixval eq 9999999) THEN badpixval = minimum
  gooddata = WHERE( vimage ne badpixval)
  IF(gooddata(0) EQ -1) THEN BEGIN
     PRINT,'The image does not have any good pixels.'
     GOTO, begin_histogram
  ENDIF
  IF((STRMID(name,0,6) EQ 'ZOOMED') OR (area_selected EQ 'true') $
     OR (STRMID(name,0,4) EQ 'PROJ')) THEN BEGIN
     minimum=MIN(vimage(gooddata))
     maximum=MAX(vimage(gooddata))
  ENDIF
  med = MEDIAN(vimage(gooddata))
  sz = SIZE(gooddata)
  dim1 = sz(1)
  IF(sz(1) EQ 1) THEN BEGIN
    sd = 0.
    mean = MAX(area(good))
  ENDIF ELSE BEGIN
         sd = STDEV(vimage(gooddata))
         mean = TOTAL(vimage(gooddata),/double)/N_ELEMENTS(vimage(gooddata))
        ENDELSE
  PRINT,underline('Image/Area Profile:')
  PRINT,' Minimum     =',STRING(minimum,'(g12.5)'),' ',units
  PRINT,' Maximum     =',STRING(maximum,'(g12.5)'),' ',units
  PRINT,' Mean        =',STRING(mean,'(g12.5)'),' ',units
  PRINT,' Median      =',STRING(med,'(g12.5)'),' ',units
  PRINT,' Stan Dev    =',STRING(sd,'(g12.5)')
  PRINT,' # of pixels = ',STRTRIM(STRING(dim1),2)
  IF(minimum EQ maximum) THEN BEGIN
     PRINT,'Image has no dynamic range.'  
     IF(uimage_version EQ 2) THEN GOTO, begin_histogram ELSE RETURN
  ENDIF
  q = ['Do you want to make a Histogram?','Yes','No']
  i = one_line_menu(q,init=1)
  IF( i EQ 2) THEN GOTO, print_jou
  q = ['Define a special MIN & MAX for the Histogram?','Yes','No']
  i = one_line_menu(q,init=2)
  IF ( i eq 1) THEN BEGIN
     err=1
     WHILE(err NE 0) DO BEGIN
        val = 0
        WHILE(val EQ 0) DO BEGIN
           READ,'Enter MIN :  ',str
           val = validnum(str)
           IF(val EQ 1) THEN minimum=FLOAT(str) $
                   ELSE PRINT,'Invalid number, please re-enter.'
        ENDWHILE
        IF(minimum LT maximum) THEN err=0 ELSE print,'The minimum is >= maximum' 
     ENDWHILE 
     err=1
     WHILE(err NE 0) DO BEGIN
	val=0
        WHILE(val EQ 0) DO BEGIN
           READ,'Enter MAX :  ',str
           val = validnum(str)
           IF(val EQ 1) THEN maximum=FLOAT(str) $
                   ELSE PRINT,'Invalid number, please re-enter.'
        ENDWHILE
        IF(maximum GT minimum) THEN err=0 ELSE print,'The maximum is <= minimum'
     ENDWHILE 
  ENDIF
; 
; Get the number of Bins from the user
;
  err = 1
  range = maximum - minimum
  WHILE ( err NE 0) DO BEGIN
     val=0
     WHILE(val EQ 0) DO BEGIN
        READ,'Enter number of Bins :  ',str
        val = validnum(str)
        IF(val EQ 1) THEN nbin=FIX(str) $
                   ELSE PRINT,'Invalid number, please re-enter.'
     ENDWHILE
     IF( (nbin ge 2) and (nbin le 1024)) THEN err = 0 ELSE BEGIN
        PRINT,'Please enter the No. of bins such that: 2 <= No of bins <=1024'
     ENDELSE
  ENDWHILE
  binsz = range/nbin
;
; Call the Histogram function with the given Max,Min and number of bins
;
  exstr = 'HISTOGRAM(vimage(gooddata), '+$
		'binsize=binsz, min=minimum, max=maximum)'
  j = EXECUTE('hist = '+exstr)
  istatus = EXECUTE('wno='+ name + '.window')
;
; Setup menu for user to select Linear or Logarithmic plotting
;
  list = STRARR(2)
  q = ['Choose type of plot','Linear','Logarithmic']
  flag = one_line_menu(q,init=1) 
  lflag = 0
  IF(flag eq 2) THEN BEGIN
     lflag = 1
     hist = hist > 10 
  ENDIF
  sz = SIZE(hist)
  dim1 = nbin
  x =  findgen(dim1)
  ypos = 0.88
;
; Plot Linear or Logarithmic based on the selection
;
  name_graph = setup_graph(x,hist,dim1,logflag=lflag,psymbol=10,$
  		title = title,x_title='Bin number',y_title='Number of pixels',$
		label_title='Physical units',topmin=minimum,topmax=maximum,$
		label_size=1.2,nblocks_x=nblocks_x,nblocks_y=nblocks_y)
  plot_graph,name_graph
;
; Resize graph if necessary 
;
  IF((first_resize EQ 1) AND ((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN'))) $
    THEN BEGIN
     sel=UMENU(['Resize graph?',$
       'Yes / Resize the graph','No / Exit this menu'],title=0)
     IF(sel EQ 1) THEN BEGIN
       resize_graph,name_graph
       j = EXECUTE('nblocks_x = '+name_graph+'.nblocks_x')
       j = EXECUTE('nblocks_y = '+name_graph+'.nblocks_y')
     ENDIF
     first_resize=0
   ENDIF
;
; print journaling information
;
print_jou:
  IF((area_selected eq 'true') AND ((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN'))) $
    THEN BEGIN
     menu = ['Remove area boundary from image?','Yes','No']
     sel = one_line_menu(menu,init=1)
     IF(sel EQ 1) THEN IF(STRMID(name,0,6) NE 'ZOOMED') $
        THEN xdisplay,name,/rescale ELSE show_zoomed,name
  ENDIF
  IF(journal_on EQ 1) THEN BEGIN
    PRINTF,luj,menu_title
    istatus = EXECUTE('xtitle='+ name + '.title')
    PRINTF,luj,'  operand:  '+xtitle
    IF(area_selected EQ 'true') THEN BEGIN
       PRINTF,luj,'  an area of the image was selected'
    ENDIF ELSE BEGIN
       PRINTF,luj,'  the whole image was selected'
    ENDELSE
    PRINTF,luj,'  minimum = ',STRTRIM(minimum,2)
    PRINTF,luj,'  maximum = ',STRTRIM(maximum,2)
    PRINTF,luj,'  mean    = ',STRTRIM(mean,2)
    PRINTF,luj,'  median  = ',STRTRIM(med,2)
    PRINTF,luj,'  standard deviation = ',STRTRIM(sd,2)
    PRINTF,luj,'  # of pixels = ',STRTRIM(dim1,2)
    IF(exstr NE ' ') THEN BEGIN
       PRINTF,luj,'  the following command was executed:'
       PRINTF,luj,'  result ='+exstr
       PRINTF,luj,'  where:'
       PRINTF,luj,"    image is the data in the operand's structure"
       PRINTF,luj,'    gooddata is the good pixels'
       PRINTF,luj,'    binsz = ',STRTRIM(binsz,2)
       PRINTF,luj,'    hist --> '+title
    ENDIF
    PRINTF,luj,'----------------------------------------'+$
                 '--------------------------------------'
  ENDIF
  IF(uimage_version EQ 2) THEN GOTO, begin_histogram ELSE RETURN
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



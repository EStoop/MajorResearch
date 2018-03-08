PRO usvdfit, name
;+
;  USVDFIT - a UIMAGE-specific routine.  This is UIMAGE's driver routine
;  for the SVDFIT procedure. This fits a data set (e.g. a spectrum) with 
;  a polynomial. 
;#
;  Written by J Newmark
;  SPR 10939  May 11 93 Creation.
;  SPR 11106  Jun 28 93 Clean up I/O and add overplotting of fits. J Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON journal,journal_on,luj
  COMMON xscreen,magnify,block_usage,scrsiz
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = '1-D Polynomial Fitting'
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Put up a menu of potential operands.
;  ------------------------------------
show_menu:
  name = select_object(/graph,/help,/exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no appropriate objects available.', /CONT
    RETURN
  ENDIF
;
;  An operand has been selected - make sure it is valid.
;  -----------------------------------------------------
have_name:
parent = name
j = EXECUTE('title = ' + name + '.title')
j = EXECUTE('badpixval = ' + name + '.badpixval')
j = EXECUTE('linkweight= ' + name + '.linkweight')
j = EXECUTE('scatter = ' + name + '.scatter')
j=EXECUTE('num = ' +name+ '.num')
j=EXECUTE('swindow = ' + name + '.window')
IF (scatter EQ 1) THEN BEGIN
  j=EXECUTE('win_orig1 = ' + name + '.win_orig1')
  j=EXECUTE('win_orig2 = ' + name + '.win_orig2')
  name1=get_name(win_orig1)
  name2=get_name(win_orig2)
  j = EXECUTE('badpixval1 = ' + name1 + '.badpixval')
  j = EXECUTE('badpixval2 = ' + name2 + '.badpixval')
  j = EXECUTE('good = WHERE((' + name1 + '.data NE badpixval1) AND (' + $
      name2 + '.data NE badpixval2))')
  j=EXECUTE('freq = ' + name1 +'.data(good)')
  j=EXECUTE('spec = ' + name2 +'.data(good)')
ENDIF ELSE BEGIN
 j=EXECUTE('data = ' + name + '.data')
 freq=data(0:num-1,0)
 spec=data(0:num-1,1)
ENDELSE
j=EXECUTE('funits = ' + name + '.x_title')
j=EXECUTE('fsunit = ' + name + '.y_title')
freq=DOUBLE(freq)
spec=DOUBLE(spec)
scale=MAX(freq)
freq=freq/scale
wsel=0
IF( linkweight NE -1) THEN BEGIN
  wsel=1
  namew=get_name(linkweight)
  j = EXECUTE('titlew = ' + namew + '.title')
  PRINT, 'Weights will be taken from ' + bold(titlew)
  j = EXECUTE('badvalw = ' + namew + '.badpixval')
  j = EXECUTE('weight = ' + namew + '.data(0:num-1,1)')
  goodw = WHERE(weight NE badvalw)
  badw = WHERE(weight EQ badvalw)
  IF(goodw(0) NE -1) THEN BEGIN
    minw = MIN(weight(goodw))
    IF(badw(0) NE -1) THEN weight(badw) = minw
  ENDIF                              
ENDIF ELSE BEGIN
  PRINT,'There are no weights associated with the chosen object.'
  READ,'Enter an estimate for the Sigma: ',sig
    weight=spec*0. + 1./sig^2
ENDELSE
bad=WHERE(spec EQ badpixval)
IF (bad(0) NE -1) THEN weight(bad)=0.
weight=DOUBLE(weight)
;        
; Define various input parameters.
;
pix_parm=BYTARR(7)
plot_res=0
plot_fit=0
get_parms:
parmenu = ['Define Input/Output Parameters ',$
          'Degree of polynomial(mandatory)', $
          'Overlay fit (mandatory for scatter plots)?',$
          'Output residual of fit (not appropriate for scatter plots)?',$
          'Output chisq?', 'Output Covariance?','Output Sigma?',$
          'Output singular points?','Mask any frequencies?',$
          'HELP', 'DONE ']
inparm = UMENU(parmenu, title=0, init=inparm)
CASE inparm OF
  1: BEGIN
    READ,'Enter the degree of polynomial to use in fit: ',M
    M = M + 1
    pix_parm(0)=1
    GOTO, get_parms
    END
  2: BEGIN
    Print,'The fit will be overlayed on the graph'
    yfit=1
    plot_fit=1
    pix_parm(1)=1
    GOTO, get_parms
    END
  3: BEGIN
    Print,'A graph of the residual will be created'
    yfit=1
    pix_parm(1)=1
    plot_res=1
    GOTO, get_parms
    END
  4: BEGIN
    Print,'Chi-square will be listed.'
    chisq=1
    pix_parm(2)=1
    GOTO, get_parms
    END
  5: BEGIN
    Print,'Covariance matrix will be listed.'
    covar=1
    pix_parm(3)=1
    GOTO, get_parms
    END
  6: BEGIN
    Print,"Sigma's will be listed."
    variance=1
    pix_parm(4)=1
    GOTO, get_parms
    END
  7: BEGIN
    Print,'Singular points will be listed.'
    singular=1
    pix_parm(5)=1
    GOTO, get_parms
    END
  8: BEGIN
    index=''
    PRINT,'The range of the frequency/wavelength array is 1 -',FIX(num)
    READ,'Enter the masked element array values',index
    mask,index,weight,mask
    weight=mask
    pix_parm(6)=1
    GOTO, get_parms
    END
  9: BEGIN
      uimage_help, menu_title
      GOTO, get_parms
    END
  10: PRINT,'Parameters are defined, proceeding..'
 ENDCASE
; 
;set up proper strings for call
;
exstr= 'fit= SVDFIT(freq,spec,M,weight=weight'
; Input values
IF (pix_parm(1)) THEN exstr=exstr+',yfit=yfit'
IF (pix_parm(2)) THEN exstr=exstr+',chisq=chisq'
IF (pix_parm(3)) THEN exstr=exstr+',covar=covar'
IF (pix_parm(4)) THEN exstr=exstr+',variance=variance'
IF (pix_parm(5)) THEN exstr=exstr+',singular=singular'
exstr=exstr+')'
;
;  Call SVDFIT
;  ------------
 j = EXECUTE(exstr)
 freq=freq*scale
;
;  Store the result in a UIMAGE data-object.
;  -----------------------------------------
;
; Output for graphs
;
IF (plot_fit EQ 1) THEN BEGIN
  IF ((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
   WSET,swindow
   WSHOW,swindow
  ENDIF
  IF (scatter EQ 1) THEN BEGIN
      plot,freq,spec,psym=1
      oplot,freq(sort(freq)),yfit(sort(yfit))
  ENDIF ELSE BEGIN
      plot,freq,spec
      oplot,freq,yfit,linestyle=1
  ENDELSE
;  For non-X-window terminals, wait for a key-press before continuing.
;  -------------------------------------------------------------------
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    CASE !D.NAME OF
       'TEK'   : BEGIN
		   i = GET_KBRD(1)
                   text
                   PRINT, ' '
		 END      
       'REGIS' : BEGIN
		   i = GET_KBRD(1)
                   ERASE
                   PRINT, ' '
		 END
       ELSE    :
    ENDCASE
  ENDIF
ENDIF
IF ((plot_res EQ 1) AND (scatter NE 1)) THEN BEGIN
  residual=spec-yfit
  ntitle='Residual = ' + title +'- Fit'
  ntitle=append_number(ntitle)
  gname= setup_graph(freq,residual,num,title=ntitle,x_title=funits,$
         y_title=sunit)
  plot_graph,gname
ENDIF
PRINT,'Fit coefficients (A0 + A1*X + A2*X^2 + ... + AN*X^N)'
for i=0,m-1 do BEGIN
  var=STRTRIM('  A',0)+STRTRIM(STRING(i),2) + ' = ' +string(fit(i))
  IF (pix_parm(4) EQ 1) THEN var=var + ' +/- ' + STRING(SQRT(variance(i)))
  PRINT,var
endfor
IF (pix_parm(2) EQ 1) THEN PRINT,'Chisq = '+STRING(chisq)
IF (pix_parm(5) EQ 1) THEN PRINT,'Singular points = '+STRING(singular)
IF (pix_parm(3) EQ 1) THEN BEGIN
   PRINT,'Covar matrix = '
   FOR i=0,m-1 DO PRINT,STRING(covar(i*m:(i+1)*m-1))
ENDIF
;
;  If journaling is enabled, then write info to the journal file.
;  --------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title
    PRINTF, luj, '  the following statement was executed:'
    PRINTF, luj, '  ' + STRMID(exstr, 0, 76)
    exstr = STRMID(exstr, 76, STRLEN(exstr))
    WHILE(STRLEN(exstr) GE 1) DO BEGIN
      PRINTF, luj, '    ' + STRMID(exstr, 0, 74)
      exstr = STRMID(exstr, 74, STRLEN(exstr))
    ENDWHILE
    IF(wsel EQ 1) THEN BEGIN
      PRINTF, luj, '    weights came from:  ' + titlew
    ENDIF
    PRINTF,luj,'Fit coefficients (A0 + A1*X + A2*X^2 + ... + AN*X^N)'
    for i=0,m-1 do BEGIN
      var=STRTRIM('  A',0)+STRTRIM(STRING(i),2) + ' = ' +string(fit(i))
      IF (pix_parm(4) EQ 1) THEN var=var + ' +/- ' + STRING(SQRT(variance(i)))
      PRINTF,luj,var
   endfor
   IF (pix_parm(2) EQ 1) THEN PRINTF,luj,'Chisq = '+STRING(chisq)
   IF (pix_parm(5) EQ 1) THEN PRINTF,luj,'Singular points = '+STRING(singular)
   IF (pix_parm(3) EQ 1) THEN BEGIN
      PRINTF,luj,'Covar matrix = '
      FOR i=0,m-1 DO PRINTF,luj,STRING(covar(i*m:(i+1)*m-1))
    ENDIF
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
  ENDIF
  IF(uimage_version EQ 2) THEN GOTO, show_menu
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


